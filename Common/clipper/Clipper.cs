/*******************************************************************************
* Author    :  Angus Johnson                                                   *
* Version   :  10.0 (release candidate 1) - also known as Clipper2             *
* Date      :  19 February 2022                                                *
* Website   :  http://www.angusj.com                                           *
* Copyright :  Angus Johnson 2010-2022                                         *
* Purpose   :  This module contains simple functions that will likely cover    *
*              most polygon boolean and offsetting needs, while also avoiding  *
*              the complexities of the other modules in this clipping library. *
*Thanks     :  Special thanks to Thong Nguyen (https://nguyen.mn/) and to      *
*              Guus Kuiper (https://www.guuskuiper.nl/) for their invaluable   *
*              assistance in this C# Clipper port.                             *
*License   :  http://www.boost.org/LICENSE_1_0.txt                             *
*******************************************************************************/


using System.Collections.Generic;

namespace ClipperLib2
{
	using Path = List<Point64>;
	using Paths = List<List<Point64>>;
	using PathD = List<PointD>;
	using PathsD = List<List<PointD>>;

  public class ClipperFunc
  {
    public static Paths Intersect(Paths subject, Paths clip, FillRule fillRule)
    {
      return BooleanOp(ClipType.Intersection, fillRule, subject, clip);
    }

    public static PathsD Intersect(PathsD subject, PathsD clip, FillRule fillRule)
    {
      return BooleanOp(ClipType.Intersection, fillRule, subject, clip);
    }

    public static Paths Union(Paths subject, Paths clip, FillRule fillRule)
    {
      return BooleanOp(ClipType.Union, fillRule, subject, clip);
    }

    public static PathsD Union(PathsD subject, PathsD clip, FillRule fillRule)
    {
      return BooleanOp(ClipType.Union, fillRule, subject, clip);
    }

    public static Paths Difference(Paths subject, Paths clip, FillRule fillRule)
    {
      return BooleanOp(ClipType.Difference, fillRule, subject, clip);
    }

    public static PathsD Difference(PathsD subject, PathsD clip, FillRule fillRule)
    {
      return BooleanOp(ClipType.Difference, fillRule, subject, clip);
    }

    public static Paths Xor(Paths subject, Paths clip, FillRule fillRule)
    {
      return BooleanOp(ClipType.Xor, fillRule, subject, clip);
    }

    public static PathsD Xor(PathsD subject, PathsD clip, FillRule fillRule)
    {
      return BooleanOp(ClipType.Xor, fillRule, subject, clip);
    }

    public static Paths BooleanOp(ClipType clipType, FillRule fillRule, Paths subject, Paths clip)
    {
      if (subject == null) return null;
        Clipper c = new Clipper();
      c.AddPaths(subject, PathType.Subject);
      if (clip != null)
        c.AddPaths(clip, PathType.Clip);
      c.Execute(clipType, fillRule, out Paths solution);
      return solution;
    }

    public static PathsD BooleanOp(ClipType clipType, FillRule fillRule, PathsD subject, PathsD clip)
    {
      ClipperD c = new ClipperD();
      c.AddSubject(subject);
      c.AddClip(clip);
      c.Execute(clipType, fillRule, out PathsD solution);
      return solution;
    }

    public static PathsD InflatePaths(PathsD paths, double delta, bool isOpen)
    {
      ClipperOffset co = new ClipperOffset();
      if (isOpen)
        co.AddPaths(paths, JoinType.Round, EndType.Round);
      else
        co.AddPaths(paths, JoinType.Round, EndType.Polygon);
      return co.Execute(delta);
    }

    public static double Area(Path path)
    {
      double a = 0.0;
      int cnt = path.Count, j = cnt - 1;
      if (cnt < 3) return 0.0;
      for (int i = 0; i < cnt; i++)
      {
        double d = (path[j].X + path[i].X);
        a += d * (path[j].Y - path[i].Y);
        j = i;
      }
      return -a * 0.5;
    }

    public static double Area(Paths paths)
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
      int cnt = path.Count, j = cnt - 1;
      if (cnt < 3) return 0.0;
      for (int i = 0; i < cnt; i++)
      {
        double d = (path[j].x + path[i].x);
        a += d * (path[j].y - path[i].y);
        j = i;
      }
      return -a * 0.5;
    }

    public static double Area(PathsD paths)
    {
      double a = 0.0;
      int cnt = paths.Count;
      for (int i = 0; i < cnt; i++)
        a += Area(paths[i]);
      return a;
    }

    public static Path ScalePath(PathD path, double scale)
    {
      if (scale == 0) scale = 1;
      int cnt = path.Count;
      Path res = new Path(cnt);
      for (int i = 0; i < cnt; i++)
        res.Add(new Point64(path[i].x * scale, path[i].y * scale));
      return res;
    }

    public static Paths ScalePaths(PathsD paths, double scale)
    {
      if (scale == 0) scale = 1;
      int cnt = paths.Count;
      Paths res = new Paths(cnt);
      for (int i = 0; i < cnt; i++)
        res.Add(ScalePath(paths[i], scale));
      return res;
    }

    public static PathD ScalePath(Path path, double scale)
    {
      if (scale == 0) scale = 1;
      int cnt = path.Count;
      PathD res = new PathD(cnt);
      for (int i = 0; i < cnt; i++)
        res.Add(new PointD(path[i].X * scale, path[i].Y * scale));
      return res;
    }

    public static PathsD ScalePaths(Paths paths, double scale)
    {
      if (scale == 0) scale = 1;
      int cnt = paths.Count;
      PathsD res = new PathsD(cnt);
      for (int i = 0; i < cnt; i++)
        res.Add(ScalePath(paths[i], scale));
      return res;
    }

    public static PathD PathD(Path path)
    {
      int cnt = path.Count;
      PathD res = new PathD(cnt);
      for (int i = 0; i < cnt; i++)
        res.Add(new PointD(path[i].X, path[i].Y));
      return res;
    }

    public static PathsD PathsD(Paths paths)
    {
      int cnt = paths.Count;
      PathsD res = new PathsD(cnt);
      for (int i = 0; i < cnt; i++)
        res.Add(PathD(paths[i]));
      return res;
    }

    public static double Sqr(double value) { return value * value; }

    public static bool PointsNearEqual(PointD pt1, PointD pt2, double distanceSqrd)
    {
      return Sqr(pt1.x - pt2.x) + Sqr(pt1.y - pt2.y) < distanceSqrd;
    }

    public static PathD StripNearDuplicates(PathD path, double minEdgeLength, bool isClosedPath)
    {
      int cnt = path.Count;
      PathD result = new PathD(cnt);
      if (cnt == 0) return result;
      PointD lastPt = path[0];
      result.Add(lastPt);
      double minLengthSqrd = minEdgeLength * minEdgeLength;
      for (int i = 1; i < cnt; i++)
        if (!PointsNearEqual(lastPt, path[i], minLengthSqrd))
        {
          lastPt = path[i];
          result.Add(lastPt);
        }
      if (isClosedPath && PointsNearEqual(lastPt, result[0], minLengthSqrd)) result.RemoveAt(result.Count - 1);
      return result;
    }

    private static void AddPolyNodeToPaths(PolyPath polyPath, Paths paths)
    {
      if (polyPath.Polygon.Count > 0) paths.Add(polyPath.Polygon);
      for (int i = 0; i < polyPath.ChildCount; i++)
        AddPolyNodeToPaths((PolyPath)polyPath._childs[i], paths);
    }

    public static Paths PolyTreeToPaths(PolyTree polyTree)
    {
      Paths result = new Paths();
      for (int i = 0; i < polyTree.ChildCount; i++)
        AddPolyNodeToPaths((PolyPath)polyTree._childs[i], result);
      return result;
    }

    public static void AddPolyNodeToPathsD(PolyPathD polyPath, PathsD paths)
    {
      if (polyPath.Polygon.Count > 0) paths.Add(polyPath.Polygon);
      for (int i = 0; i < polyPath.ChildCount; i++)
        AddPolyNodeToPathsD((PolyPathD)polyPath._childs[i], paths);
    }

    public static PathsD PolyTreeToPathsD(PolyTreeD polyTree)
    {
      PathsD result = new PathsD();
      for (int i = 0; i < polyTree.ChildCount; i++)
        AddPolyNodeToPathsD((PolyPathD)polyTree._childs[i], result);
      return result;
    }

  }

} //namespace