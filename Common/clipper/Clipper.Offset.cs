/*******************************************************************************
* Author    :  Angus Johnson                                                   *
* Date      :  11 October 2025                                                 *
* Website   :  https://www.angusj.com                                          *
* Copyright :  Angus Johnson 2010-2025                                         *
* Purpose   :  Path Offset (Inflate/Shrink)                                    *
* License   :  https://www.boost.org/LICENSE_1_0.txt                           *
*******************************************************************************/

#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;

namespace Clipper2Lib
{
  public enum JoinType
  {
    Miter,
    Square,
    Bevel,
    Round
  }

  public enum EndType
  {
    Polygon,
    Joined,
    Butt,
    Square,
    Round
  }

  public class ClipperOffset
  {

    private class Group
    {
      internal Paths64 inPaths;
      internal JoinType joinType;
      internal EndType endType;
      internal bool pathsReversed;
      internal int outerSign, outerIdx;
      internal readonly List<(int sign, UInt128 magnitude)> areas = new();
      internal bool[] isStray = Array.Empty<bool>();
      internal bool IsStray(int i) => isStray.Length != 0 && isStray[i];

      public Group(Paths64 paths, JoinType joinType, EndType endType = EndType.Polygon)
      {
        this.joinType = joinType;
        this.endType = endType;
        bool joined = endType == EndType.Polygon || endType == EndType.Joined;
        inPaths = new Paths64(paths.Count);
        foreach (Path64 path in paths)
        {
          Path64 clean = Clipper.StripDuplicates(path, joined);
          if (clean.Count != 0) inPaths.Add(clean);
        }
        if (endType != EndType.Polygon) return;
        UInt128 largest = 0;
        for (int k = 0; k < inPaths.Count; k++)
        {
          Path64 path = inPaths[k];
          Int128 area = 0;
          // Match upstream's exact, wrapping 128-bit shoelace sum.
          unchecked
          {
            for (int i = 0, j = path.Count - 1; i < path.Count; j = i++)
              area += (Int128)path[j].X * path[i].Y - (Int128)path[i].X * path[j].Y;
          }
          int sign = area.CompareTo(Int128.Zero);
          UInt128 magnitude = unchecked((UInt128)(area < 0 ? -area : area));
          areas.Add((sign, magnitude));
          if (magnitude <= largest) continue;
          largest = magnitude;
          outerSign = sign;
          outerIdx = k;
        }
      }
    }

    private const double Tolerance = 1.0E-12;

    // Clipper2 approximates arcs by using series of relatively short straight
    //line segments. And logically, shorter line segments will produce better arc
    // approximations. But very short segments can degrade performance, usually
    // with little or no discernable improvement in curve quality. Very short
    // segments can even detract from curve quality, due to the effects of integer
    // rounding. Since there isn't an optimal number of line segments for any given
    // arc radius (that perfectly balances curve approximation with performance),
    // arc tolerance is user defined. Nevertheless, when the user doesn't define
    // an arc tolerance (ie leaves alone the 0 default value), the calculated
    // default arc tolerance (offset_radius / 500) generally produces good (smooth)
    // arc approximations without producing excessively small segment lengths.
    // See also: https://github.com/PolyClipping/Clipper2/wiki/Trigonometry
    private const double arc_const = 0.002; // <-- 1/500

    private readonly List<Group> _groupList = new List<Group>();
    private Path64 pathOut = new Path64();
    private readonly PathD _normals = new PathD();
    private Paths64 _solution = new Paths64();
    private PolyTree64? _solutionTree;

    private double _groupDelta; //*0.5 for open paths; *-1.0 for negative areas
    private double _delta;
    private double _mitLimSqr;
    private double _stepsPerRad;
    private double _stepSin;
    private double _stepCos;
    private JoinType _joinType;
    public double ArcTolerance { get; set; }
    public bool MergeGroups { get; set; }
    public double MiterLimit { get; set; }
    public bool PreserveCollinear { get; set; }
    public bool ReverseSolution { get; set; }

    public delegate double DeltaCallback64(Path64 path,
      PathD path_norms, int currPt, int prevPt);
    public DeltaCallback64? DeltaCallback { get; set; }

#if USINGZ
    internal void ZCB(Point64 bot1, Point64 top1,
        Point64 bot2, Point64 top2, ref Point64 ip)
    {
      if (bot1.Z != 0 &&
        ((bot1.Z == bot2.Z) || (bot1.Z == top2.Z))) ip.Z = bot1.Z;
      else if (bot2.Z != 0 && bot2.Z == top1.Z) ip.Z = bot2.Z;
      else if (top1.Z != 0 && top1.Z == top2.Z) ip.Z = top1.Z;
      else ZCallback?.Invoke(bot1, top1, bot2, top2, ref ip);
    }
    public ClipperBase.ZCallback64? ZCallback { get; set; }
#endif
    public ClipperOffset(double miterLimit = 2.0,
      double arcTolerance = 0.0, bool
      preserveCollinear = false, bool reverseSolution = false)
    {
      MiterLimit = miterLimit;
      ArcTolerance = arcTolerance;
      MergeGroups = true;
      PreserveCollinear = preserveCollinear;
      ReverseSolution = reverseSolution;
#if USINGZ
      ZCallback = null;
#endif
    }
    public void Clear()
    {
      _groupList.Clear();
    }

    public void AddPath(Path64 path, JoinType joinType, EndType endType)
    {
      int cnt = path.Count;
      if (cnt == 0) return;
      Paths64 pp = new Paths64(1) { path };
      AddPaths(pp, joinType, endType);
    }

    public void AddPaths(Paths64 paths, JoinType joinType, EndType endType)
    {
      int cnt = paths.Count;
      if (cnt == 0) return;
      _groupList.Add(new Group(paths, joinType, endType));
    }

    private int CalcSolutionCapacity()
    {
      int result = 0;
      foreach (Group g in _groupList)
        result += (g.endType == EndType.Joined) ? g.inPaths.Count * 2 : g.inPaths.Count;
      return result;
    }

    private static bool IsInside(Path64 path, Func<Point64, PointInPolygonResult> locate)
    {
      foreach (Point64 pt in path)
      {
        var pip = locate(pt);
        if (pip != PointInPolygonResult.IsOn) return pip == PointInPolygonResult.IsInside;
      }
      return true;
    }

    private static bool HasVertexInside(Path64 path, Func<Point64, PointInPolygonResult> locate)
    {
      bool allOn = true;
      foreach (Point64 pt in path)
      {
        var pip = locate(pt);
        if (pip == PointInPolygonResult.IsInside) return true;
        if (pip == PointInPolygonResult.IsOutside) allOn = false;
      }
      return allOn;
    }

    private void SetGroupOrientations()
    {
      var paths = new List<(UInt128 size, int group, int path)>();
      var bounds = new List<Rect64>();
      int[] ends = new int[_groupList.Count], outers = new int[_groupList.Count];
      bool ccw = false, cw = false, holes = false;
      for (int i = 0; i < _groupList.Count; i++)
      {
        Group g = _groupList[i];
        g.pathsReversed = g.outerSign < 0;
        g.isStray = Array.Empty<bool>();
        ccw |= g.outerSign > 0;
        cw |= g.outerSign < 0;
        for (int j = 0; j < g.areas.Count; j++)
        {
          if (g.areas[j].sign == 0) continue;
          holes |= g.areas[j].sign != g.outerSign;
          if (j == g.outerIdx) outers[i] = paths.Count;
          paths.Add((g.areas[j].magnitude, i, j));
          bounds.Add(InternalClipper.GetBounds(g.inPaths[j]));
        }
        ends[i] = paths.Count;
      }
      bool mixed = ccw && cw;
      if (!mixed && !holes) return;
      IEnumerable<int> order = Enumerable.Range(0, paths.Count);
      if (mixed) order = order.OrderByDescending(i => paths[i].size)
        .ThenByDescending(i => paths[i].path == _groupList[paths[i].group].outerIdx);
      Path64 PathOf(int i) => _groupList[paths[i].group].inPaths[paths[i].path];
      int SignOf(int i) => _groupList[paths[i].group].areas[paths[i].path].sign;
      bool IsOuter(int i) => (SignOf(i) < 0) == _groupList[paths[i].group].pathsReversed;
      BoundsTree? tree = null;
      BoundsTree GetTree()
      {
        if (tree != null) return tree;
        tree = new BoundsTree(bounds.ToArray());
        for (int i = 0; i < paths.Count; i++) tree.SetActive(i, IsOuter(i));
        return tree;
      }
      var indexes = new PolygonIndex?[paths.Count];
      PointInPolygonResult Locate(int i, Point64 pt)
      {
        Path64 polygon = PathOf(i);
        if (polygon.Count < 64) return Clipper.PointInPolygon(pt, polygon);
        return (indexes[i] ??= new PolygonIndex(polygon)).Locate(pt);
      }
      foreach (int k in order)
      {
        int groupIdx = paths[k].group;
        Group g = _groupList[groupIdx];
        Path64 path = PathOf(k);
        Rect64 pathBounds = bounds[k];
        if (paths[k].path == g.outerIdx)
        {
          if (!mixed) continue;
          var enclosers = new List<int>();
          GetTree().Query(r => r.Contains(pathBounds), i =>
          {
            if (paths[i].group != groupIdx && paths[i].size > paths[k].size &&
              !_groupList[paths[i].group].IsStray(paths[i].path)) enclosers.Add(i);
          });
          UInt128? innermost = null;
          bool isHole = false;
          foreach (int i in enclosers.OrderBy(i => paths[i].size))
          {
            if (innermost.HasValue && paths[i].size > innermost.Value) break;
            if (!IsInside(path, pt => Locate(i, pt))) continue;
            innermost = paths[i].size;
            if (IsOuter(i) && SignOf(i) != g.outerSign) isHole = true;
          }
          if (!isHole) continue;
          g.pathsReversed = !g.pathsReversed;
          for (int i = groupIdx == 0 ? 0 : ends[groupIdx - 1]; i < ends[groupIdx]; i++)
            GetTree().SetActive(i, IsOuter(i));
          continue;
        }
        if (IsOuter(k)) continue;
        int ownOuter = outers[groupIdx];
        bool enclosed = IsOuter(ownOuter) && HasVertexInside(path, pt => Locate(ownOuter, pt));
        if (!enclosed)
          GetTree().Query(r => !enclosed && r.Intersects(pathBounds), i =>
          {
            if (paths[i].size >= paths[k].size) enclosed = HasVertexInside(path, pt => Locate(i, pt));
          }, true);
        if (enclosed) continue;
        if (g.isStray.Length == 0) g.isStray = new bool[g.inPaths.Count];
        g.isStray[paths[k].path] = true;
      }
    }

    private Clipper64 InitClipper()
    {
      var c = new Clipper64 { PreserveCollinear = PreserveCollinear, Compatibility = Compatibility };
#if USINGZ
      c.ZCallback = ZCB;
#endif
      return c;
    }

    private void KeepGrownRegions()
    {
      var c = InitClipper();
      c.AddSubject(_solution[^1]);
      _solution.RemoveAt(_solution.Count - 1);
      var grown = new Paths64();
      c.Execute(ClipType.Union, FillRule.Positive, grown);
      _solution.AddRange(grown);
    }

    /// <summary>Retain the first polygon group's input winding for existing consumers.
    /// Set false for upstream C++ positive outer winding. ReverseSolution still applies.</summary>
    private bool? _preserveInputOrientation;
    public bool PreserveInputOrientation
    {
      get => _preserveInputOrientation ?? Compatibility == ClipperCompatibility.Legacy;
      set => _preserveInputOrientation = value;
    }

    /// <summary>Legacy retains established square joins and clipping rounding.
    /// Upstream enables the current C++ coordinate calculations.</summary>
    public ClipperCompatibility Compatibility { get; set; } = ClipperCompatibility.Legacy;

    private void ExecuteInternal(double delta)
    {
      if (_groupList.Count == 0) return;
      _solution.EnsureCapacity(CalcSolutionCapacity());

      SetGroupOrientations();
      if (Math.Abs(delta) < 0.5)
      {
        foreach (Group group in _groupList)
        {
          if (group.endType != EndType.Polygon) continue;
          for (int i = 0; i < group.inPaths.Count; i++)
          {
            if (group.IsStray(i)) continue;
            Path64 path = new Path64(group.inPaths[i]);
            if (group.pathsReversed) path.Reverse();
            _solution.Add(path);
          }
        }
      }
      else
      {
        _delta = delta;
        _mitLimSqr = MiterLimit <= 1 ? 2.0 : 2.0 / Clipper.Sqr(MiterLimit);
        foreach (Group group in _groupList) DoGroupOffset(group);
      }
      if (_solution.Count == 0) return;
      Clipper64 c = InitClipper();
      bool reversed = PreserveInputOrientation &&
        (_groupList.FirstOrDefault(g => g.endType == EndType.Polygon)?.outerSign < 0);
      c.ReverseSolution = ReverseSolution != reversed;
      c.AddSubject(_solution);
      if (_solutionTree != null) c.Execute(ClipType.Union, FillRule.Positive, _solutionTree);
      else c.Execute(ClipType.Union, FillRule.Positive, _solution);
    }

    public void Execute(double delta, Paths64 solution)
    {
      solution.Clear();
      _solutionTree = null;
      _solution = solution;
      ExecuteInternal(delta);
    }

    public void Execute(double delta, PolyTree64 solutionTree)
    {
      solutionTree.Clear();
      _solutionTree = solutionTree;
      _solution = new Paths64();
      try { ExecuteInternal(delta); }
      finally { _solutionTree = null; }
    }


    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal static PointD GetUnitNormal(Point64 pt1, Point64 pt2)
    {
      double dx = (pt2.X - pt1.X);
      double dy = (pt2.Y - pt1.Y);
      if ((dx == 0) && (dy == 0)) return new PointD();

      double f = 1.0 / Math.Sqrt(dx * dx + dy * dy);
      dx *= f;
      dy *= f;

      return new PointD(dy, -dx);
    }

    public void Execute(DeltaCallback64 deltaCallback, Paths64 solution)
    {
      DeltaCallback = deltaCallback;
      Execute(1.0, solution);
    }    
    
  [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static PointD TranslatePoint(PointD pt, double dx, double dy)
    {
#if USINGZ
      return new PointD(pt.x + dx, pt.y + dy, pt.z);
#else
      return new PointD(pt.x + dx, pt.y + dy);
#endif
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static PointD ReflectPoint(PointD pt, PointD pivot)
    {
#if USINGZ
      return new PointD(pivot.x + (pivot.x - pt.x), pivot.y + (pivot.y - pt.y), pt.z);
#else
      return new PointD(pivot.x + (pivot.x - pt.x), pivot.y + (pivot.y - pt.y));
#endif
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static bool AlmostZero(double value, double epsilon = 0.001)
    {
      return Math.Abs(value) < epsilon;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double Hypotenuse(double x, double y)
    {
      return Math.Sqrt(Math.Pow(x, 2) + Math.Pow(y, 2));
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static PointD NormalizeVector(PointD vec)
    {
	    double h = Hypotenuse(vec.x, vec.y);
	    if (AlmostZero(h)) return new PointD(0,0);
        double inverseHypot = 1 / h;
	    return new PointD(vec.x* inverseHypot, vec.y* inverseHypot);
    }


    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static PointD GetAvgUnitVector(PointD vec1, PointD vec2)
    {
	    return NormalizeVector(new PointD(vec1.x + vec2.x, vec1.y + vec2.y));
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private Point64 GetPerpendic(Point64 pt, PointD norm)
    {
#if USINGZ
      return new Point64(pt.X + norm.x * _groupDelta,
        pt.Y + norm.y * _groupDelta, pt.Z);
#else
      return new Point64(pt.X + norm.x * _groupDelta,
        pt.Y + norm.y * _groupDelta);
#endif
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private PointD GetPerpendicD(Point64 pt, PointD norm)
    {
#if USINGZ
      return new PointD(pt.X + norm.x * _groupDelta,
        pt.Y + norm.y * _groupDelta, pt.Z);
#else
      return new PointD(pt.X + norm.x * _groupDelta,
        pt.Y + norm.y * _groupDelta);
#endif
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void DoBevel(Path64 path, int j, int k)
    {
      Point64 pt1, pt2;
      if (j == k)
      {
        double absDelta = Math.Abs(_groupDelta);
#if USINGZ
        pt1 = new Point64(
          path[j].X - absDelta * _normals[j].x, 
          path[j].Y - absDelta * _normals[j].y, path[j].Z);
        pt2 = new Point64(
          path[j].X + absDelta * _normals[j].x, 
          path[j].Y + absDelta * _normals[j].y, path[j].Z);
#else
        pt1 = new Point64(
          path[j].X - absDelta * _normals[j].x,
          path[j].Y - absDelta * _normals[j].y);
        pt2 = new Point64(
          path[j].X + absDelta * _normals[j].x,
          path[j].Y + absDelta * _normals[j].y);
#endif
      }
      else
      {
#if USINGZ
        pt1 = new Point64(
          path[j].X + _groupDelta * _normals[k].x,
          path[j].Y + _groupDelta * _normals[k].y, path[j].Z);
        pt2 = new Point64(
          path[j].X + _groupDelta * _normals[j].x,
          path[j].Y + _groupDelta * _normals[j].y, path[j].Z);
#else
        pt1 = new Point64(
          path[j].X + _groupDelta * _normals[k].x,
          path[j].Y + _groupDelta * _normals[k].y);
        pt2 = new Point64(
          path[j].X + _groupDelta * _normals[j].x,
          path[j].Y + _groupDelta * _normals[j].y);
#endif
      }
      pathOut.Add(pt1);
      pathOut.Add(pt2);
    }

    private void SquareIntersect(PointD a, PointD b, PointD c, PointD d, out PointD pt)
    {
      if (Compatibility == ClipperCompatibility.Upstream)
      {
        InternalClipper.GetLineIntersectPt(a, b, c, d, out pt);
        return;
      }
      // Preserve the established infinite-line square join for existing designs.
      if (InternalClipper.IsAlmostZero(a.x - b.x))
      {
        if (InternalClipper.IsAlmostZero(c.x - d.x)) { pt = new PointD(); return; }
        double m = (d.y - c.y) / (d.x - c.x);
        pt = new PointD(a.x, m * a.x + (c.y - m * c.x));
      }
      else if (InternalClipper.IsAlmostZero(c.x - d.x))
      {
        double m = (b.y - a.y) / (b.x - a.x);
        pt = new PointD(c.x, m * c.x + (a.y - m * a.x));
      }
      else
      {
        double m1 = (b.y - a.y) / (b.x - a.x), m2 = (d.y - c.y) / (d.x - c.x);
        if (InternalClipper.IsAlmostZero(m1 - m2)) { pt = new PointD(); return; }
        double b1 = a.y - m1 * a.x, b2 = c.y - m2 * c.x;
        double x = (b2 - b1) / (m1 - m2);
        pt = new PointD(x, m1 * x + b1);
      }
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void DoSquare(Path64 path, int j, int k)
    {
      PointD vec;
      if (j == k)
      {
        vec = new PointD(_normals[j].y, -_normals[j].x);
      }
      else
      {
        vec = GetAvgUnitVector(
          new PointD(-_normals[k].y, _normals[k].x),
          new PointD(_normals[j].y, -_normals[j].x));
      }

      double absDelta = Math.Abs(_groupDelta);
      // now offset the original vertex delta units along unit vector
      PointD ptQ = new PointD(path[j]);
      ptQ = TranslatePoint(ptQ, absDelta * vec.x, absDelta * vec.y);

      // get perpendicular vertices
      PointD pt1 = TranslatePoint(ptQ, _groupDelta * vec.y, _groupDelta * -vec.x);
      PointD pt2 = TranslatePoint(ptQ, _groupDelta * -vec.y, _groupDelta * vec.x);
      // get 2 vertices along one edge offset
      PointD pt3 = GetPerpendicD(path[k], _normals[k]);

      if (j == k)
      {
        PointD pt4 = new PointD(
          pt3.x + vec.x * _groupDelta,
          pt3.y + vec.y * _groupDelta);
        SquareIntersect(pt1, pt2, pt3, pt4, out PointD pt);
#if USINGZ
        pt.z = ptQ.z;
#endif    
        //get the second intersect point through reflecion
        pathOut.Add(new Point64(ReflectPoint(pt, ptQ)));
        pathOut.Add(new Point64(pt));
      }
      else
      {
        PointD pt4 = GetPerpendicD(path[j], _normals[k]);
        SquareIntersect(pt1, pt2, pt3, pt4, out PointD pt);
#if USINGZ
        pt.z = ptQ.z;
#endif
        pathOut.Add(new Point64(pt));
        //get the second intersect point through reflecion
        pathOut.Add(new Point64(ReflectPoint(pt, ptQ)));
      }
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void DoMiter(Path64 path, int j, int k, double cosA)
    {
      double q = _groupDelta / (cosA + 1);
#if USINGZ
      pathOut.Add(new Point64(
          path[j].X + (_normals[k].x + _normals[j].x) * q,
          path[j].Y + (_normals[k].y + _normals[j].y) * q,
          path[j].Z));
#else
      pathOut.Add(new Point64(
          path[j].X + (_normals[k].x + _normals[j].x) * q,
          path[j].Y + (_normals[k].y + _normals[j].y) * q));
#endif
    }

    private void SetArcSteps(double delta)
    {
      double absDelta = Math.Abs(delta);
      double arcTol = ArcTolerance > Tolerance ? Math.Min(absDelta, ArcTolerance) : absDelta * arc_const;
      double steps = Math.Min(Math.PI / Math.Acos(1 - arcTol / absDelta), absDelta * Math.PI);
      _stepSin = Math.Sin(2 * Math.PI / steps);
      _stepCos = Math.Cos(2 * Math.PI / steps);
      if (delta < 0) _stepSin = -_stepSin;
      _stepsPerRad = steps / (2 * Math.PI);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void DoRound(Path64 path, int j, int k, double angle)
    {
      if (DeltaCallback != null) SetArcSteps(_groupDelta);

      Point64 pt = path[j];
      PointD offsetVec = new PointD(_normals[k].x * _groupDelta, _normals[k].y * _groupDelta);
      if (j == k) offsetVec.Negate();
#if USINGZ
      pathOut.Add(new Point64(pt.X + offsetVec.x, pt.Y + offsetVec.y, pt.Z));
#else
      pathOut.Add(new Point64(pt.X + offsetVec.x, pt.Y + offsetVec.y));
#endif
      int steps = (int) Math.Ceiling(_stepsPerRad * Math.Abs(angle));
      for (int i = 1; i < steps; i++) // ie 1 less than steps
      {
        offsetVec = new PointD(offsetVec.x * _stepCos - _stepSin * offsetVec.y,
            offsetVec.x * _stepSin + offsetVec.y * _stepCos);
#if USINGZ
        pathOut.Add(new Point64(pt.X + offsetVec.x, pt.Y + offsetVec.y, pt.Z));
#else
        pathOut.Add(new Point64(pt.X + offsetVec.x, pt.Y + offsetVec.y));
#endif
      }
      pathOut.Add(GetPerpendic(pt, _normals[j]));
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void BuildNormals(Path64 path)
    {
      int cnt = path.Count;
      _normals.Clear();
      if (cnt == 0) return;
      _normals.EnsureCapacity(cnt);
      for (int i = 0; i < cnt - 1; i++)
        _normals.Add(GetUnitNormal(path[i], path[i + 1]));
      _normals.Add(GetUnitNormal(path[cnt - 1], path[0]));
    }

    private void OffsetPoint(Group group, Path64 path, int j, ref int k)
    {
      if (path[j] == path[k]) { k = j; return; }

      // Let A = change in angle where edges join
      // A == 0: ie no change in angle (flat join)
      // A == PI: edges 'spike'
      // sin(A) < 0: right turning
      // cos(A) < 0: change in angle is more than 90 degree
      double sinA = InternalClipper.CrossProduct(_normals[j], _normals[k]);
      double cosA = InternalClipper.DotProduct(_normals[j], _normals[k]);
      if (sinA > 1.0) sinA = 1.0;
      else if (sinA < -1.0) sinA = -1.0;

      if (DeltaCallback != null)
      { 
        _groupDelta = DeltaCallback(path, _normals, j, k);
        if (group.pathsReversed) _groupDelta = -_groupDelta;
      }
      if (Math.Abs(_groupDelta) < Tolerance)
      {
        pathOut.Add(path[j]);
        return;
      }

      if (cosA > -0.999 && (sinA * _groupDelta < 0)) // test for concavity first (#593)
      {
        // is concave
        // by far the simplest way to construct concave joins, especially those joining very 
        // short segments, is to insert 3 points that produce negative regions. These regions 
        // will be removed later by the finishing union operation. This is also the best way 
        // to ensure that path reversals (ie over-shrunk paths) are removed.
        pathOut.Add(GetPerpendic(path[j], _normals[k]));
        pathOut.Add(path[j]); // (#405, #873, #916)
        pathOut.Add(GetPerpendic(path[j], _normals[j]));
      }
      else if ((cosA > 0.999) && (_joinType != JoinType.Round))
      {
        // almost straight - less than 2.5 degree (#424, #482, #526 & #724) 
        DoMiter(path, j, k, cosA);
      }
      else switch (_joinType)
      {
        // miter unless the angle is sufficiently acute to exceed ML
        case JoinType.Miter when cosA > _mitLimSqr - 1:
          DoMiter(path, j, k, cosA);
          break;
        case JoinType.Miter:
          DoSquare(path, j, k);
          break;
        case JoinType.Round:
          DoRound(path, j, k, Math.Atan2(sinA, cosA));
          break;
        case JoinType.Bevel:
          DoBevel(path, j, k);
          break;
        default:
          DoSquare(path, j, k);
          break;
      }

      k = j;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void OffsetPolygon(Group group, Path64 path)
    {
      pathOut = new Path64();
      int cnt = path.Count, prev = cnt - 1;
      for (int i = 0; i < cnt; i++)
        OffsetPoint(group, path, i, ref prev);
      if (group.pathsReversed) pathOut.Reverse();
      _solution.Add(pathOut);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void OffsetOpenJoined(Group group, Path64 path)
    {
      OffsetPolygon(group, path);
      path = Clipper.ReversePath(path);
      BuildNormals(path);
      OffsetPolygon(group, path);
    }

    private void OffsetOpenPath(Group group, Path64 path, EndType endType)
    {
      pathOut = new Path64();
      int highI = path.Count - 1;

      if (DeltaCallback != null) 
        _groupDelta = DeltaCallback(path, _normals, 0, 0);

      // do the line start cap
      if (Math.Abs(_groupDelta) < Tolerance)
        pathOut.Add(path[0]);
      else
        switch (endType)
        {
          case EndType.Butt:
            DoBevel(path, 0, 0);
            break;
          case EndType.Round:
            DoRound(path, 0, 0, Math.PI);
            break;
          default:
            DoSquare(path, 0, 0);
            break;
        }

      // offset the left side going forward
      for (int i = 1, k = 0; i < highI; i++)
        OffsetPoint(group, path, i, ref k);

      // reverse normals ...
      for (int i = highI; i > 0; i--)
        _normals[i] = new PointD(-_normals[i - 1].x, -_normals[i - 1].y);
      _normals[0] = _normals[highI];

      if (DeltaCallback != null)
        _groupDelta = DeltaCallback(path, _normals, highI, highI);
      // do the line end cap
      if (Math.Abs(_groupDelta) < Tolerance)
        pathOut.Add(path[highI]);
      else
        switch (endType)
        {
          case EndType.Butt:
            DoBevel(path, highI, highI);
            break;
          case EndType.Round:
            DoRound(path, highI, highI, Math.PI);
            break;
          default:
            DoSquare(path, highI, highI);
            break;
        }

      // offset the left side going back
      for (int i = highI -1, k = highI; i > 0; i--)
        OffsetPoint(group, path, i, ref k);

      _solution.Add(pathOut);
    }

    private static bool IsFlat(Path64 path)
    {
      int i = 1;
      while (i < path.Count && path[i] == path[0]) i++;
      if (i == path.Count) return true;
      Point64 pt1 = path[i++];
      for (; i < path.Count; i++)
        if (InternalClipper.CrossProductSign(path[0], pt1, path[i]) != 0) return false;
      return true;
    }

    private void DoGroupOffset(Group group)
    {
      if (group.endType == EndType.Polygon)
      {
        double delta = group.outerSign != 0 ? _delta : Math.Abs(_delta);
        _groupDelta = group.pathsReversed ? -delta : delta;
      }
      else _groupDelta = Math.Abs(_delta);
      double absDelta = Math.Abs(_groupDelta);
      _joinType = group.joinType;
      if (group.joinType == JoinType.Round || group.endType == EndType.Round) SetArcSteps(_groupDelta);
      for (int i = 0; i < group.inPaths.Count; i++)
      {
        if (group.IsStray(i)) continue;
        Path64 p = group.inPaths[i];
        pathOut = new Path64();
        BuildNormals(p);
        if (p.Count == 1)
        {
          double pointDelta = group.endType == EndType.Polygon ? _delta : Math.Abs(_delta);
          if (DeltaCallback != null)
          {
            pointDelta = DeltaCallback(p, _normals, 0, 0);
            _groupDelta = group.pathsReversed ? -pointDelta : pointDelta;
            absDelta = Math.Abs(_groupDelta);
          }
          if (pointDelta < 1) continue;
          Point64 pt = p[0];
          if (group.joinType == JoinType.Round)
          {
            if (DeltaCallback != null) SetArcSteps(_groupDelta);
            int steps = _stepsPerRad > 0 ? (int)Math.Ceiling(_stepsPerRad * 2 * Math.PI) : 0;
            pathOut = Clipper.Ellipse(pt, absDelta, absDelta, steps);
          }
          else
          {
            int d = (int)Math.Ceiling(absDelta);
            pathOut = new Rect64(pt.X - d, pt.Y - d, pt.X + d, pt.Y + d).AsPath();
          }
#if USINGZ
          pathOut = InternalClipper.SetZ(pathOut, pt.Z);
#endif
          _solution.Add(pathOut);
          continue;
        }
        EndType endType = group.endType;
        if (p.Count == 2 && endType == EndType.Joined)
          endType = group.joinType == JoinType.Round ? EndType.Round : EndType.Square;
        if (endType == EndType.Polygon)
        {
          bool flat = IsFlat(p);
          if (flat && DeltaCallback == null && _delta < 0) continue;
          OffsetPolygon(group, p);
          if (flat && DeltaCallback != null) KeepGrownRegions();
        }
        else if (endType == EndType.Joined) OffsetOpenJoined(group, p);
        else OffsetOpenPath(group, p, endType);
      }
    }
  }
}
