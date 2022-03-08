using ClipperLib2;
using geoLib;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using utility;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static GeoLibPoint[] Rotate(GeoLibPoint pivot, GeoLibPoint[] pointList, double angleDegree)
    {
        return pRotate(pivot, pointList, angleDegree);
    }

    private static GeoLibPoint[] pRotate(GeoLibPoint pivot, GeoLibPoint[] pointList, double angleDegree)
    {
        switch (Math.Abs(angleDegree))
        {
            // essentially a zero rotation clamp at 0.01 degrees
            case < 1E-2:
                return pointList;
        }

        int pointListLength = pointList.Length;
        GeoLibPoint[] returnList = new GeoLibPoint[pointListLength];

#if !GWSINGLETHREADED
        Parallel.For(0, pointListLength, i =>
#else
            for (Int32 i = 0; i < pointList.Count(); i++)
#endif
            {
                returnList[i] = Rotate(pivot, pointList[i], angleDegree);
            }
#if !GWSINGLETHREADED
        );
#endif
        return returnList;
    }

    public static GeoLibPointF[] Rotate(GeoLibPointF pivot, GeoLibPointF[] pointList, double angleDegree)
    {
        return pRotate(pivot, pointList, angleDegree);
    }

    private static GeoLibPointF[] pRotate(GeoLibPointF pivot, GeoLibPointF[] pointList, double angleDegree)
    {
        switch (Math.Abs(angleDegree))
        {
            // essentially a zero rotation clamp at 0.01 degrees
            case < 1E-2:
                return pointList;
        }

        int pointListLength = pointList.Length;
        GeoLibPointF[] returnList = new GeoLibPointF[pointListLength];

#if !GWSINGLETHREADED
        Parallel.For(0, pointListLength, i =>
#else
            for (Int32 i = 0; i < pointList.Count(); i++)
#endif
            {
                returnList[i] = Rotate(pivot, pointList[i], angleDegree);
            }
#if !GWSINGLETHREADED
        );
#endif
        return returnList;
    }

    public static List<GeoLibPointF> Rotate(GeoLibPointF pivot, List<GeoLibPointF> pointList, double angleDegree)
    {
        return pRotate(pivot, pointList, angleDegree);
    }

    private static List<GeoLibPointF> pRotate(GeoLibPointF pivot, List<GeoLibPointF> pointList, double angleDegree)
    {
        switch (Math.Abs(angleDegree))
        {
            // essentially a zero rotation clamp at 0.01 degrees
            case < 1E-2:
                return pointList.ToList();
        }

        return pointList.Select(t => Rotate(pivot, t, angleDegree)).ToList();
    }

    public static GeoLibPointF Rotate(GeoLibPointF pivot, GeoLibPointF point, double angleDegree)
    {
        return pRotate(pivot, point, angleDegree);
    }

    private static GeoLibPointF pRotate(GeoLibPointF pivot, GeoLibPointF point, double angleDegree)
    {
        switch (Math.Abs(angleDegree))
        {
            // essentially a zero rotation clamp at 0.01 degrees
            case < 1E-2:
                return point;
        }

        double angle = Utils.toRadians(angleDegree);
        double x = pivot.X + ((point.X - pivot.X) * Math.Cos(angle) -
                              (point.Y - pivot.Y) * Math.Sin(angle));
        double y = pivot.Y + ((point.X - pivot.X) * Math.Sin(angle) +
                              (point.Y - pivot.Y) * Math.Cos(angle));

        GeoLibPointF rotated = new(x, y);
        return rotated;
    }

    public static GeoLibPoint Rotate(GeoLibPoint pivot, GeoLibPoint point, double angleDegree)
    {
        return pRotate(pivot, point, angleDegree);
    }

    private static GeoLibPoint pRotate(GeoLibPoint pivot, GeoLibPoint point, double angleDegree)
    {
        switch (Math.Abs(angleDegree))
        {
            // essentially a zero rotation clamp at 0.01 degrees
            case < 1E-2:
                return point;
        }

        double angle = Utils.toRadians(angleDegree);
        double x = pivot.X + ((point.X - pivot.X) * Math.Cos(angle) -
                              (point.Y - pivot.Y) * Math.Sin(angle));
        double y = pivot.Y + ((point.X - pivot.X) * Math.Sin(angle) +
                              (point.Y - pivot.Y) * Math.Cos(angle));

        GeoLibPoint rotated = new(x, y);
        return rotated;
    }

    public static Point64 Rotate(Point64 pivot, Point64 point, double angleDegree)
    {
        return pRotate(pivot, point, angleDegree);
    }

    private static Point64 pRotate(Point64 pivot, Point64 point, double angleDegree)
    {
        switch (Math.Abs(angleDegree))
        {
            // essentially a zero rotation clamp at 0.01 degrees
            case < 1E-2:
                return point;
        }

        double angle = Utils.toRadians(angleDegree);
        double x = pivot.X + ((point.X - pivot.X) * Math.Cos(angle) -
                              (point.Y - pivot.Y) * Math.Sin(angle));
        double y = pivot.Y + ((point.X - pivot.X) * Math.Sin(angle) +
                              (point.Y - pivot.Y) * Math.Cos(angle));

        Point64 rotated = new(x, y, point.Z);
        return rotated;
    }

    public static GeoLibPointF[] flip(bool H, bool V, bool alignX, bool alignY, GeoLibPointF pivot, GeoLibPointF[] source)
    {
        return pFlip(H, V, alignX, alignY, pivot, source);
    }

    public static List<GeoLibPointF> flip(bool H, bool V, bool alignX, bool alignY, GeoLibPointF pivot, List<GeoLibPointF> source)
    {
        return pFlip(H, V, alignX, alignY, pivot, source.ToArray()).ToList();
    }

    private static GeoLibPointF[] pFlip(bool H, bool V, bool alignX, bool alignY, GeoLibPointF pivot, GeoLibPointF[] source)
    {
        int sLength = source.Length;
        GeoLibPointF[] ret = new GeoLibPointF[sLength];
#if !GWSINGLETHREADED
        Parallel.For(0, sLength, pt =>
#else
            for (Int32 pt = 0; pt < source.Length; pt++)
#endif
            {
                double newX = source[pt].X;
                switch (H)
                {
                    case true:
                    {
                        newX = -newX;
                        switch (alignX)
                        {
                            case true:
                                newX += 2 * pivot.X;
                                break;
                        }

                        break;
                    }
                }
                double newY = source[pt].Y;
                switch (V)
                {
                    case true:
                    {
                        newY = -newY;
                        switch (alignY)
                        {
                            case true:
                                newY += 2 * pivot.Y;
                                break;
                        }

                        break;
                    }
                }
                ret[pt] = new GeoLibPointF(newX, newY);
            }
#if !GWSINGLETHREADED
        );
#endif
        switch (H ^ V)
        {
            case true:
                ret = ret.Reverse().ToArray(); // preserve ordering.
                break;
        }

        return ret;
    }

}