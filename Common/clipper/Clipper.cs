﻿/*******************************************************************************
* Author    :  Angus Johnson                                                   *
* Version   :  10.0 (beta) - also known as Clipper2                            *
* Date      :  4 May 2022                                                      *
* Website   :  http://www.angusj.com                                           *
* Copyright :  Angus Johnson 2010-2022                                         *
* Purpose   :  This module contains simple functions that will likely cover    *
*              most polygon boolean and offsetting needs, while also avoiding  *
*              the inherent complexities of the other modules.                 *
* Thanks    :  Special thanks to Thong Nguyen, Guus Kuiper, Phil Stopford,     *
*           :  and Daniel Gosnell for their invaluable assistance with C#.     *
* License   :  http://www.boost.org/LICENSE_1_0.txt                            *
*******************************************************************************/


#nullable enable
using System.Collections.Generic;
using System.Runtime.CompilerServices;

namespace Clipper2Lib
{

  //PRE-COMPILER CONDITIONAL ...
  //USINGZ: For user defined Z-coordinates. See Clipper.SetZ

  using Path64 = List<Point64>;
  using Paths64 = List<List<Point64>>;
  using PathD = List<PointD>;
  using PathsD = List<List<PointD>>;

  public static class ClipperFunc
  {
    public static Paths64 Intersect(Paths64 subject, Paths64 clip, FillRule fillRule)
    {
      return BooleanOp(ClipType.Intersection, fillRule, subject, clip);
    }

    public static PathsD Intersect(PathsD subject, PathsD clip, FillRule fillRule)
    {
      return BooleanOp(ClipType.Intersection, fillRule, subject, clip);
    }

    public static Paths64 Union(Paths64 subject, FillRule fillRule)
    {
      return BooleanOp(ClipType.Union, fillRule, subject, null);
    }

    public static Paths64 Union(Paths64 subject, Paths64 clip, FillRule fillRule)
    {
      return BooleanOp(ClipType.Union, fillRule, subject, clip);
    }

    public static PathsD Union(PathsD subject, FillRule fillRule)
    {
      return BooleanOp(ClipType.Union, fillRule, subject, null);
    }

    public static PathsD Union(PathsD subject, PathsD clip, FillRule fillRule)
    {
      return BooleanOp(ClipType.Union, fillRule, subject, clip);
    }

    public static Paths64 Difference(Paths64 subject, Paths64 clip, FillRule fillRule)
    {
      return BooleanOp(ClipType.Difference, fillRule, subject, clip);
    }

    public static PathsD Difference(PathsD subject, PathsD clip, FillRule fillRule)
    {
      return BooleanOp(ClipType.Difference, fillRule, subject, clip);
    }

    public static Paths64 Xor(Paths64 subject, Paths64 clip, FillRule fillRule)
    {
      return BooleanOp(ClipType.Xor, fillRule, subject, clip);
    }

    public static PathsD Xor(PathsD subject, PathsD clip, FillRule fillRule)
    {
      return BooleanOp(ClipType.Xor, fillRule, subject, clip);
    }

    public static Paths64 BooleanOp(ClipType clipType, FillRule fillRule, 
      Paths64? subject, Paths64? clip)
    {
      Paths64 solution = new Paths64();
      if (subject == null) return solution;
      Clipper c = new Clipper();
      c.AddPaths(subject, PathType.Subject);
      if (clip != null)
        c.AddPaths(clip, PathType.Clip);
      c.Execute(clipType, fillRule, solution);
      return solution;
    }

    public static PathsD BooleanOp(ClipType clipType, FillRule fillRule,
        PathsD subject, PathsD? clip, int roundingDecimalPrecision = 0)
    {
      PathsD solution = new PathsD();
      ClipperD c = new ClipperD(roundingDecimalPrecision);
      c.AddSubject(subject);
      if (clip != null)
        c.AddClip(clip);
      c.Execute(clipType, fillRule, solution);
      return solution;
    }

    public static Paths64 InflatePaths(Paths64 paths, double delta, JoinType joinType, EndType endType)
    {
      ClipperOffset co = new ClipperOffset();
      co.AddPaths(paths, joinType, endType);
      PathsD tmp = co.Execute(delta);
      return Paths64(tmp);
    }

    public static PathsD InflatePaths(PathsD paths, double delta, JoinType joinType, EndType endType)
    {
      ClipperOffset co = new ClipperOffset();
      co.AddPaths(paths, joinType, endType);
      return co.Execute(delta);
    }

    public static double Area(Path64 path)
    {
      double a = 0.0;
      int cnt = path.Count;
      if (cnt < 3) return 0.0;
      Point64 prevPt = path[cnt - 1];
      for (int i = 0; i < path.Count; i++)
      {
        a += (double) (prevPt.Y - path[i].Y) * (prevPt.X + path[i].X);
        prevPt = path[i];
      }
#if REVERSE_ORIENTATION
      return a * -0.5;
#else
      return a * 0.5;
#endif
    }

    public static double Area(Paths64 paths)
    {
      double a = 0.0;
      int cnt = paths.Count;
      for (int i = 0; i < cnt; i++)
        a += Area(paths[i]);
      return a;
    }

    public static double Area(PathD path)
    {
      double a = 0.0;
      int cnt = path.Count;
      if (cnt < 3) return 0.0;
      PointD prevPt = path[cnt - 1];
      for (int i = 0; i < cnt; i++)
      {
        a += (prevPt.y - path[i].y) * (prevPt.x + path[i].x);
        prevPt = path[i];
      }
#if REVERSE_ORIENTATION
      return a * -0.5;
#else
      return a * 0.5;
#endif
    }

    public static double Area(PathsD paths)
    {
      double a = 0.0;
      for (int i = 0; i < paths.Count; i++)
        a += Area(paths[i]);
      return a;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static bool IsClockwise(Path64 poly)
    {
      return Area(poly) >= 0;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static bool IsClockwise(PathD poly)
    {
      return Area(poly) >= 0;
    }

    public static Path64 OffsetPath(Path64 path, long dx, long dy)
    {
      Path64 result = new Path64(path.Count);
      for (int i = 0; i < path.Count; i++)
        result.Add(new Point64(path[i].X + dx, path[i].Y + dy));
      return result;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static PointD ScalePoint(Point64 pt, double scale) 
    {
      PointD result = new PointD()
      {
        x = pt.X * scale,
        y = pt.Y * scale,
      };
      return result;
    }

    public static Path64 ScalePath(Path64 path, double scale)
    {
      if (scale == 1) return path;
      Path64 result = new Path64(path.Count);
      for (int i = 0; i < path.Count; i++)
        result.Add(new Point64(path[i], scale));
      return result;
    }

    public static Paths64 ScalePaths(Paths64 paths, double scale)
    {
      if (scale == 1) return paths;
      Paths64 result = new Paths64(paths.Count);
      for (int i = 0; i < paths.Count; i++)
        result.Add(ScalePath(paths[i], scale));
      return result;
    }

    public static PathD ScalePath(PathD path, double scale)
    {
      if (scale == 1) return path;
      PathD result = new PathD(path.Count);
      for (int i = 0; i < path.Count; i++)
        result.Add(new PointD(path[i], scale));
      return result;
    }

    public static PathsD ScalePaths(PathsD paths, double scale)
    {
      if (scale == 1) return paths;
      PathsD result = new PathsD(paths.Count);
      for (int i = 0; i < paths.Count; i++)
        result.Add(ScalePath(paths[i], scale));
      return result;
    }

    //Unlike ScalePath, both ScalePath64 & ScalePathD also involve type conversion
    public static Path64 ScalePath64(PathD path, double scale)
    {
      int cnt = path.Count;
      Path64 res = new Path64(cnt);
      for (int i = 0; i < cnt; i++)
        res.Add(new Point64(path[i].x * scale, path[i].y * scale));
      return res;
    }

    public static Paths64 ScalePaths64(PathsD paths, double scale)
    {
      int cnt = paths.Count;
      Paths64 res = new Paths64(cnt);
      for (int i = 0; i < cnt; i++)
        res.Add(ScalePath64(paths[i], scale));
      return res;
    }

    public static PathD ScalePathD(Path64 path, double scale)
    {
      int cnt = path.Count;
      PathD res = new PathD(cnt);
      for (int i = 0; i < cnt; i++)
        res.Add(new PointD(path[i].X * scale, path[i].Y * scale));
      return res;
    }

    public static PathsD ScalePathsD(Paths64 paths, double scale)
    {
      int cnt = paths.Count;
      PathsD res = new PathsD(cnt);
      for (int i = 0; i < cnt; i++)
        res.Add(ScalePathD(paths[i], scale));
      return res;
    }

    //The static functions Path64 and PathD convert path types without scaling
    public static Path64 Path64(PathD path)
    {
      Path64 result = new Path64(path.Count);
      for (int i = 0; i < path.Count; i++)
        result.Add(new Point64(path[i]));
      return result;
    }

    public static Paths64 Paths64(PathsD paths)
    {
      Paths64 result = new Paths64(paths.Count);
      for (int i = 0; i < paths.Count; i++)
        result.Add(Path64(paths[i]));
      return result;
    }

    public static PathsD PathsD(Paths64 paths)
    {
      PathsD result = new PathsD(paths.Count);
      for (int i = 0; i < paths.Count; i++)
        result.Add(PathD(paths[i]));
      return result;
    }

    public static PathD PathD(Path64 path)
    {
      PathD result = new PathD(path.Count);
      for (int i = 0; i < path.Count; i++)
        result.Add(new PointD(path[i]));
      return result;
    }

    public static Paths64 OffsetPaths(Paths64 paths, long dx, long dy)
    {
      Paths64 result = new Paths64(paths.Count);
      for (int i = 0; i < paths.Count; i++)
        result.Add(OffsetPath(paths[i], dx, dy));
      return result;
    }

    public static PathD OffsetPath(PathD path, long dx, long dy)
    {
      PathD result = new PathD(path.Count);
      for (int i = 0; i < path.Count; i++)
        result.Add(new PointD(path[i].y + dx, path[i].y + dy));
      return result;
    }

    public static PathsD OffsetPaths(PathsD paths, long dx, long dy)
    {
      PathsD result = new PathsD(paths.Count);
      for (int i = 0; i < paths.Count; i++)
        result.Add(OffsetPath(paths[i], dx, dy));
      return result;
    }

    public static Path64 ReversePath(Path64 path)
    {
      Path64 result = new Path64(path.Count);
      for (int i = path.Count - 1; i >= 0; i--)
        result.Add(path[i]);
      return result;
    }

    public static PathD ReversePath(PathD path)
    {
      PathD result = new PathD(path.Count);
      for (int i = path.Count - 1; i >= 0; i--)
        result.Add(path[i]);
      return result;
    }

    public static Paths64 ReversePaths(Paths64 paths)
    {
      Paths64 result = new Paths64(paths.Count);
      for (int i = 0; i < paths.Count; i++)
        result.Add(ReversePath(paths[i]));
      return result;
    }

    public static PathsD ReversePaths(PathsD paths)
    {
      PathsD result = new PathsD(paths.Count);
      for (int i = 0; i < paths.Count; i++)
        result.Add(ReversePath(paths[i]));
      return result;
    }

    public static Rect64 GetBounds(Paths64 paths)
    {
      Rect64 result = new Rect64(
        long.MaxValue, long.MaxValue, long.MinValue, long.MinValue);
      for (int i = 0; i < paths.Count; i++)
        for (int j = 0; j < paths[i].Count; j++)
        {
          Point64 pt = paths[i][j];
          if (pt.X < result.left) result.left = pt.X;
          if (pt.X > result.right) result.right = pt.X;
          if (pt.Y < result.top) result.top = pt.Y;
          if (pt.Y > result.bottom) result.bottom = pt.Y;
        }
      return result.IsEmpty() ? new Rect64() : result;
    }

    public static RectD GetBounds(PathsD paths)
    {
      RectD result = new RectD(double.MaxValue, -double.MaxValue,
        -double.MaxValue, -double.MaxValue);
      for (int i = 0; i < paths.Count; i++)
        for (int j = 0; j < paths[i].Count; j++)
        {
          PointD pt = paths[i][j];
          if (pt.x < result.left) result.left = pt.x;
          if (pt.x > result.right) result.right = pt.x;
          if (pt.y < result.top) result.top = pt.y;
          if (pt.y > result.bottom) result.bottom = pt.y;
        }
      return result.IsEmpty() ? new RectD() : result;
    }

    public static Path64 MakePath(int[] arr)
    {
      int len = arr.Length / 2;
      Path64 p = new Path64(len);
      for (int i = 0; i < len; i++)
        p.Add(new Point64(arr[i * 2], arr[i * 2 + 1]));
      return p;
    }

    public static Path64 MakePath(long[] arr)
    {
      int len = arr.Length / 2;
      Path64 p = new Path64(len);
      for (int i = 0; i < len; i++)
        p.Add(new Point64(arr[i * 2], arr[i * 2 + 1]));
      return p;
    }

    public static PathD MakePath(double[] arr)
    {
      int len = arr.Length / 2;
      PathD p = new PathD(len);
      for (int i = 0; i < len; i++)
        p.Add(new PointD(arr[i * 2], arr[i * 2 + 1]));
      return p;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double Sqr(double value)
    {
      return value * value;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static bool PointsNearEqual(PointD pt1, PointD pt2, double distanceSqrd)
    {
      return Sqr(pt1.x - pt2.x) + Sqr(pt1.y - pt2.y) < distanceSqrd;
    }

    public static PathD StripNearDuplicates(PathD path,
        double minEdgeLenSqrd, bool isClosedPath)
    {
      int cnt = path.Count;
      PathD result = new PathD(cnt);
      if (cnt == 0) return result;
      PointD lastPt = path[0];
      result.Add(lastPt);
      for (int i = 1; i < cnt; i++)
        if (!PointsNearEqual(lastPt, path[i], minEdgeLenSqrd))
        {
          lastPt = path[i];
          result.Add(lastPt);
        }

      if (isClosedPath && PointsNearEqual(lastPt, result[0], minEdgeLenSqrd))
      {
        result.RemoveAt(result.Count - 1);
      }

      return result;
    }

    private static void AddPolyNodeToPaths(PolyPath polyPath, Paths64 paths)
    {
      if (polyPath.Polygon!.Count > 0)
        paths.Add(polyPath.Polygon);
      for (int i = 0; i < polyPath.ChildCount; i++)
        AddPolyNodeToPaths((PolyPath) polyPath._childs[i], paths);
    }

    public static Paths64 PolyTreeToPaths(PolyTree polyTree)
    {
      Paths64 result = new Paths64();
      for (int i = 0; i < polyTree.ChildCount; i++)
        AddPolyNodeToPaths((PolyPath) polyTree._childs[i], result);
      return result;
    }

    public static void AddPolyNodeToPathsD(PolyPathD polyPath, PathsD paths)
    {
      if (polyPath.Polygon!.Count > 0)
        paths.Add(polyPath.Polygon);
      for (int i = 0; i < polyPath.ChildCount; i++)
        AddPolyNodeToPathsD((PolyPathD) polyPath._childs[i], paths);
    }

    public static PathsD PolyTreeToPathsD(PolyTreeD polyTree)
    {
      PathsD result = new PathsD();
      for (int i = 0; i < polyTree.ChildCount; i++)
        AddPolyNodeToPathsD((PolyPathD) polyTree._childs[i], result);
      return result;
    }
  }
} //namespace