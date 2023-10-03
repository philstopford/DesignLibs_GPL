﻿/*******************************************************************************
* Author    :  Angus Johnson                                                   *
* Date      :  24 September 2023                                               *
* Website   :  http://www.angusj.com                                           *
* Copyright :  Angus Johnson 2010-2023                                         *
* Purpose   :  Path Offset (Inflate/Shrink)                                    *
* License   :  http://www.boost.org/LICENSE_1_0.txt                            *
*******************************************************************************/

using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;

namespace Clipper2Lib
{
  public enum JoinType
  {
    Miter,
    Square,
    Bevel,
    Round
  };

  public enum EndType
  {
    Polygon,
    Joined,
    Butt,
    Square,
    Round
  };

  public class ClipperOffset
  {

    private class Group
    {
      internal Paths64 inPaths;
      internal Path64 outPath;
      internal Paths64 outPaths;
      internal JoinType joinType;
      internal EndType endType;
      internal bool pathsReversed;

      public Group(Paths64 paths, JoinType joinType, EndType endType = EndType.Polygon)
      {
        inPaths = new Paths64(paths);
        this.joinType = joinType;
        this.endType = endType;
        outPath = new Path64();
        outPaths = new Paths64();
        pathsReversed = false;
      }
    }

    private static double Tolerance = 1.0E-12;
    private readonly List<Group> _groupList = new List<Group>();
    private readonly PathD _normals = new PathD();
    private readonly Paths64 _solution = new Paths64();
    private double _groupDelta; //*0.5 for open paths; *-1.0 for negative areas
    private double _delta;
    private double _mitLimSqr;
    private double _stepsPerRad;
    private double _stepSin;
    private double _stepCos;
    private JoinType _joinType;
    private EndType _endType;
    public double ArcTolerance { get; set; }
    public bool MergeGroups { get; set; }
    public double MiterLimit { get; set; }
    public bool PreserveCollinear { get; set; }
    public bool ReverseSolution { get; set; }

    public delegate double DeltaCallback64(Path64 path,
      PathD path_norms, int currPt, int prevPt);
    public ClipperOffset.DeltaCallback64? DeltaCallback { get; set; }

#if USINGZ
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

    private void ExecuteInternal(double delta)
    {
      _solution.Clear();
      if (_groupList.Count == 0) return;

      if (Math.Abs(delta) < 0.5)
      {
        foreach (Group group in _groupList)
          foreach (Path64 path in group.inPaths)
            _solution.Add(path);
      }
      else
      { 
        _delta = delta;
        _mitLimSqr = (MiterLimit <= 1 ?
          2.0 : 2.0 / Clipper.Sqr(MiterLimit));

        foreach (Group group in _groupList)
          DoGroupOffset(group);
      }
    }

    public void Execute(double delta, Paths64 solution)
    {
      solution.Clear();
      ExecuteInternal(delta);
      if (_groupList.Count == 0) return;

      // clean up self-intersections ...
      Clipper64 c = new Clipper64()
      {
        PreserveCollinear = PreserveCollinear,
        // the solution should retain the orientation of the input
        ReverseSolution = ReverseSolution != _groupList[0].pathsReversed
      };
#if USINGZ
      c.ZCallback = ZCallback;
#endif
      c.AddSubject(_solution);
      if (_groupList[0].pathsReversed)
        c.Execute(ClipType.Union, FillRule.Negative, solution);
      else
        c.Execute(ClipType.Union, FillRule.Positive, solution);
    }

    public void Execute(double delta, PolyTree64 polytree)
    {
      polytree.Clear();
      ExecuteInternal(delta);
      if (_groupList.Count == 0) return;

      // clean up self-intersections ...
      Clipper64 c = new Clipper64()
      {
        PreserveCollinear = PreserveCollinear,
        // the solution should retain the orientation of the input
        ReverseSolution = ReverseSolution != _groupList[0].pathsReversed
      };
#if USINGZ
      c.ZCallback = ZCallback;
#endif
      c.AddSubject(_solution);
      if (_groupList[0].pathsReversed)
        c.Execute(ClipType.Union, FillRule.Negative, polytree);
      else
        c.Execute(ClipType.Union, FillRule.Positive, polytree);
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
    private static void GetBoundsAndLowestPolyIdx(Paths64 paths,
      out int index, out Rect64 rec)
    {
      rec = new Rect64(false); // ie invalid rect
      long lpX = long.MinValue;
      index = -1;
      for (int i = 0; i < paths.Count; i++)
        foreach (Point64 pt in paths[i])
        {
          if (pt.Y >= rec.bottom)
          {
            if (pt.Y > rec.bottom || pt.X < lpX)
            {
              index = i;
              lpX = pt.X;
              rec.bottom = pt.Y;
            }
          }
          else if (pt.Y < rec.top) rec.top = pt.Y;
          if (pt.X > rec.right) rec.right = pt.X;
          else if (pt.X < rec.left) rec.left = pt.X;
        }
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
    private static PointD IntersectPoint(PointD pt1a, PointD pt1b, PointD pt2a, PointD pt2b)
    {
      if (InternalClipper.IsAlmostZero(pt1a.x - pt1b.x)) //vertical
      {
        if (InternalClipper.IsAlmostZero(pt2a.x - pt2b.x)) return new PointD(0, 0);
        double m2 = (pt2b.y - pt2a.y) / (pt2b.x - pt2a.x);
        double b2 = pt2a.y - m2 * pt2a.x;
        return new PointD(pt1a.x, m2* pt1a.x + b2);
      }

      if (InternalClipper.IsAlmostZero(pt2a.x - pt2b.x)) //vertical
      {
        double m1 = (pt1b.y - pt1a.y) / (pt1b.x - pt1a.x);
        double b1 = pt1a.y - m1 * pt1a.x;
        return new PointD(pt2a.x, m1* pt2a.x + b1);
      }
      else
      {
        double m1 = (pt1b.y - pt1a.y) / (pt1b.x - pt1a.x);
        double b1 = pt1a.y - m1 * pt1a.x;
        double m2 = (pt2b.y - pt2a.y) / (pt2b.x - pt2a.x);
        double b2 = pt2a.y - m2 * pt2a.x;
        if (InternalClipper.IsAlmostZero(m1 - m2)) return new PointD(0, 0);
        double x = (b2 - b1) / (m1 - m2);
        return new PointD(x, m1 * x + b1);
      }
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
    private void DoBevel(Group group, Path64 path, int j, int k)
    {
      Point64 pt1, pt2;
      if (j == k)
      {
        double absDelta = Math.Abs(_groupDelta);
        pt1 = new Point64(path[j].X - absDelta * _normals[j].x, path[j].Y - absDelta * _normals[j].y);
        pt2 = new Point64(path[j].X + absDelta * _normals[j].x, path[j].Y + absDelta * _normals[j].y);
      }
      else
      {
        pt1 = new Point64(path[j].X + _groupDelta * _normals[k].x, path[j].Y + _groupDelta * _normals[k].y);
        pt2 = new Point64(path[j].X + _groupDelta * _normals[j].x, path[j].Y + _groupDelta * _normals[j].y);
      }
      group.outPath.Add(pt1);
      group.outPath.Add(pt2);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void DoSquare(Group group, Path64 path, int j, int k)
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
        PointD pt = IntersectPoint(pt1, pt2, pt3, pt4);
#if USINGZ
        pt.z = ptQ.z;
#endif    
        //get the second intersect point through reflecion
        group.outPath.Add(new Point64(ReflectPoint(pt, ptQ)));
        group.outPath.Add(new Point64(pt));
      }
      else
      {
        PointD pt4 = GetPerpendicD(path[j], _normals[k]);
        PointD pt = IntersectPoint(pt1, pt2, pt3, pt4);
#if USINGZ
        pt.z = ptQ.z;
#endif
        group.outPath.Add(new Point64(pt));
        //get the second intersect point through reflecion
        group.outPath.Add(new Point64(ReflectPoint(pt, ptQ)));
      }
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void DoMiter(Group group, Path64 path, int j, int k, double cosA)
    {
      double q = _groupDelta / (cosA + 1);
#if USINGZ
      group.outPath.Add(new Point64(
          path[j].X + (_normals[k].x + _normals[j].x) * q,
          path[j].Y + (_normals[k].y + _normals[j].y) * q,
          path[j].Z));
#else
      group.outPath.Add(new Point64(
          path[j].X + (_normals[k].x + _normals[j].x) * q,
          path[j].Y + (_normals[k].y + _normals[j].y) * q));
#endif
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void DoRound(Group group, Path64 path, int j, int k, double angle)
    {
      if (DeltaCallback != null)
      {
        // when DeltaCallback is assigned, _groupDelta won't be constant,
        // so we'll need to do the following calculations for *every* vertex.
        double absDelta = Math.Abs(_groupDelta);
        double arcTol = ArcTolerance > 0.01 ?
          ArcTolerance :
          Math.Log10(2 + absDelta) * InternalClipper.defaultArcTolerance;
        double stepsPer360 = Math.PI / Math.Acos(1 - arcTol / absDelta);
        _stepSin = Math.Sin((2 * Math.PI) / stepsPer360);
        _stepCos = Math.Cos((2 * Math.PI) / stepsPer360);
        if (_groupDelta < 0.0) _stepSin = -_stepSin;
        _stepsPerRad = stepsPer360 / (2 * Math.PI);
      }

      Point64 pt = path[j];
      PointD offsetVec = new PointD(_normals[k].x * _groupDelta, _normals[k].y * _groupDelta);
      if (j == k) offsetVec.Negate();
#if USINGZ
      group.outPath.Add(new Point64(pt.X + offsetVec.x, pt.Y + offsetVec.y, pt.Z));
#else
      group.outPath.Add(new Point64(pt.X + offsetVec.x, pt.Y + offsetVec.y));
#endif
      int steps = (int) Math.Ceiling(_stepsPerRad * Math.Abs(angle));
      for (int i = 1; i < steps; i++) // ie 1 less than steps
      {
        offsetVec = new PointD(offsetVec.x * _stepCos - _stepSin * offsetVec.y,
            offsetVec.x * _stepSin + offsetVec.y * _stepCos);
#if USINGZ
        group.outPath.Add(new Point64(pt.X + offsetVec.x, pt.Y + offsetVec.y, pt.Z));
#else
        group.outPath.Add(new Point64(pt.X + offsetVec.x, pt.Y + offsetVec.y));
#endif
      }
      group.outPath.Add(GetPerpendic(pt, _normals[j]));
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void BuildNormals(Path64 path)
    {
      int cnt = path.Count;
      _normals.Clear();
      _normals.Capacity = cnt;

      for (int i = 0; i < cnt - 1; i++)
        _normals.Add(GetUnitNormal(path[i], path[i + 1]));
      _normals.Add(GetUnitNormal(path[cnt - 1], path[0]));
    }

    private void OffsetPoint(Group group, Path64 path, int j, ref int k)
    {
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
        group.outPath.Add(path[j]);
        return;
      }

      if (cosA > -0.99 && (sinA * _groupDelta < 0)) // test for concavity first (#593)
      {
        // is concave
        group.outPath.Add(GetPerpendic(path[j], _normals[k]));
        // this extra point is the only (simple) way to ensure that
        // path reversals are fully cleaned with the trailing clipper
        group.outPath.Add(path[j]); // (#405)
        group.outPath.Add(GetPerpendic(path[j], _normals[j]));
      }
      else if (cosA > 0.999)
        DoMiter(group, path, j, k, cosA);
      else if (_joinType == JoinType.Miter)
      {
        // miter unless the angle is so acute the miter would exceeds ML
        if (cosA > _mitLimSqr - 1) DoMiter(group, path, j, k, cosA);
        else DoSquare(group, path, j, k);
      }
      else if (cosA > 0.99 || _joinType == JoinType.Bevel)
        //angle less than 8 degrees or a squared join
        DoBevel(group, path, j, k);
      else if (_joinType == JoinType.Round)
        DoRound(group, path, j, k, Math.Atan2(sinA, cosA));
      else
        DoSquare(group, path, j, k);

      k = j;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void OffsetPolygon(Group group, Path64 path)
    {
      // when the path is contracting, make sure 
      // there is sufficient space to do so.              //#593
      //nb: this will have a small impact on performance
      double a = Clipper.Area(path);
      if ((a < 0) != (_groupDelta < 0))
      {
        Rect64 rec = Clipper.GetBounds(path);
        double offsetMinDim = Math.Abs(_groupDelta) * 2;
        if (offsetMinDim > rec.Width || offsetMinDim > rec.Height) return;
      }

      group.outPath = new Path64();
      int cnt = path.Count, prev = cnt - 1;
      for (int i = 0; i < cnt; i++)
        OffsetPoint(group, path, i, ref prev);
      group.outPaths.Add(group.outPath);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private void OffsetOpenJoined(Group group, Path64 path)
    {
      OffsetPolygon(group, path);
      path = Clipper.ReversePath(path);
      BuildNormals(path);
      OffsetPolygon(group, path);
    }

    private void OffsetOpenPath(Group group, Path64 path)
    {
      group.outPath = new Path64();
      int highI = path.Count - 1;

      if (DeltaCallback != null) 
        _groupDelta = DeltaCallback(path, _normals, 0, 0);

      // do the line start cap
      if (Math.Abs(_groupDelta) < Tolerance)
        group.outPath.Add(path[0]);
      else
        switch (_endType)
        {
          case EndType.Butt:
            DoBevel(group, path, 0, 0);
            break;
          case EndType.Round:
            DoRound(group, path, 0, 0, Math.PI);
            break;
          default:
            DoSquare(group, path, 0, 0);
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
        group.outPath.Add(path[highI]);
      else
        switch (_endType)
        {
          case EndType.Butt:
            DoBevel(group, path, highI, highI);
            break;
          case EndType.Round:
            DoRound(group, path, highI, highI, Math.PI);
            break;
          default:
            DoSquare(group, path, highI, highI);
            break;
        }

      // offset the left side going back
      for (int i = highI, k = 0; i > 0; i--)
        OffsetPoint(group, path, i, ref k);

      group.outPaths.Add(group.outPath);
    }

    private void DoGroupOffset(Group group)
    {
      if (group.endType == EndType.Polygon)
      {
        // the lowermost polygon must be an outer polygon. So we can use that as the
        // designated orientation for outer polygons (needed for tidy-up clipping)
        GetBoundsAndLowestPolyIdx(group.inPaths,
          out int lowestIdx, out Rect64 grpBounds);
        if (lowestIdx < 0) return;
        double area = Clipper.Area(group.inPaths[lowestIdx]);
        //if (area == 0) return; // this is probably unhelpful (#430)
        group.pathsReversed = (area < 0);
        if (group.pathsReversed) _groupDelta = -_delta;
        else _groupDelta = _delta;
      }
      else
      {
        group.pathsReversed = false;
        _groupDelta = Math.Abs(_delta) * 0.5;
      }
      double absDelta = Math.Abs(_groupDelta);
      _joinType = group.joinType;
      _endType = group.endType;

      if (DeltaCallback == null && 
        (group.joinType == JoinType.Round || group.endType == EndType.Round))
      {
        // calculate a sensible number of steps (for 360 deg for the given offset
        // arcTol - when fArcTolerance is undefined (0), the amount of
        // curve imprecision that's allowed is based on the size of the
        // offset (delta). Obviously very large offsets will almost always
        // require much less precision. See also offset_triginometry2.svg
        double arcTol = ArcTolerance > 0.01 ?
          ArcTolerance :              
          Math.Log10(2 + absDelta) * InternalClipper.defaultArcTolerance; 
        double stepsPer360 = Math.PI / Math.Acos(1 - arcTol / absDelta);
        _stepSin = Math.Sin((2 * Math.PI) / stepsPer360);
        _stepCos = Math.Cos((2 * Math.PI) / stepsPer360);
        if (_groupDelta < 0.0) _stepSin = -_stepSin;
        _stepsPerRad = stepsPer360 / (2 * Math.PI);
      }

      bool isJoined =
        (group.endType == EndType.Joined) ||
        (group.endType == EndType.Polygon);

      foreach (Path64 p in group.inPaths)
      {
        Path64 path = Clipper.StripDuplicates(p, isJoined);
        int cnt = path.Count;
        if ((cnt == 0) || ((cnt < 3) && (_endType == EndType.Polygon))) 
          continue;

        if (cnt == 1)
        {
          group.outPath = new Path64();
          // single vertex so build a circle or square ...
          if (group.endType == EndType.Round)
          {
            double r = absDelta;
            int steps = (int)Math.Ceiling(_stepsPerRad * 2 * Math.PI);
            group.outPath = Clipper.Ellipse(path[0], r, r, steps);
#if USINGZ
            group.outPath = InternalClipper.SetZ(group.outPath, path[0].Z);
#endif      
          }
          else
          {
            int d = (int) Math.Ceiling(_groupDelta);
            Rect64 r = new Rect64(path[0].X - d, path[0].Y - d,
              path[0].X - d, path[0].Y - d);
            group.outPath = r.AsPath();
#if USINGZ
            group.outPath = InternalClipper.SetZ(group.outPath, path[0].Z);
#endif
          }
          group.outPaths.Add(group.outPath);
        }
        else
        {
          if (cnt == 2 && group.endType == EndType.Joined)
          { 
            if (group.joinType == JoinType.Round)
              _endType = EndType.Round;
            else
              _endType = EndType.Square;
          }
          BuildNormals(path);
          if (_endType == EndType.Polygon) OffsetPolygon(group, path);
          else if (_endType == EndType.Joined) OffsetOpenJoined(group, path);
          else OffsetOpenPath(group, path);
        }
      }
      _solution.AddRange(group.outPaths);
      group.outPaths.Clear();
    }
  }

} // namespace