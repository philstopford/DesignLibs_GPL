using Clipper2Lib;
using geoLib;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using utility;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static Path64 Rotate(Point64 pivot, Path64 pointList, double angleDegree)
    {
        return pRotate(pivot, pointList, angleDegree);
    }

    private static Path64 pRotate(Point64 pivot, Path64 pointList, double angleDegree)
    {
        switch (Math.Abs(angleDegree))
        {
            // essentially a zero rotation clamp at 0.01 degrees
            case < 1E-2:
                return pointList;
        }

        int pointListLength = pointList.Capacity;
        Path64 returnList = new (pointListLength);

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

    public static PathD Rotate(PointD pivot, PathD pointList, double angleDegree)
    {
        return pRotate(pivot, pointList, angleDegree);
    }

    private static PathD pRotate(PointD pivot, PathD pointList, double angleDegree)
    {
        switch (Math.Abs(angleDegree))
        {
            // essentially a zero rotation clamp at 0.01 degrees
            case < 1E-2:
                return pointList;
        }

        int pointListLength = pointList.Capacity;
        PathD returnList = new (pointListLength);

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
    public static PointD Rotate(PointD pivot, PointD point, double angleDegree)
    {
        return pRotate(pivot, point, angleDegree);
    }

    private static PointD pRotate(PointD pivot, PointD point, double angleDegree)
    {
        switch (Math.Abs(angleDegree))
        {
            // essentially a zero rotation clamp at 0.01 degrees
            case < 1E-2:
                return point;
        }

        double angle = Utils.toRadians(angleDegree);
        double x = pivot.x + ((point.x - pivot.x) * Math.Cos(angle) -
                              (point.y - pivot.y) * Math.Sin(angle));
        double y = pivot.y + ((point.x - pivot.x) * Math.Sin(angle) +
                              (point.y - pivot.y) * Math.Cos(angle));

        PointD rotated = new(x, y);
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

        Point64 rotated = new(x, y);
        return rotated;
    }

    private static PathD pFlip(bool H, bool V, bool alignX, bool alignY, PointD pivot, PathD source)
    {
        int sLength = source.Capacity;
        PathD ret = new (sLength);
#if !GWSINGLETHREADED
        Parallel.For(0, sLength, pt =>
#else
            for (Int32 pt = 0; pt < source.Length; pt++)
#endif
            {
                double newX = source[pt].x;
                switch (H)
                {
                    case true:
                    {
                        newX = -newX;
                        switch (alignX)
                        {
                            case true:
                                newX += 2 * pivot.x;
                                break;
                        }

                        break;
                    }
                }
                double newY = source[pt].y;
                switch (V)
                {
                    case true:
                    {
                        newY = -newY;
                        switch (alignY)
                        {
                            case true:
                                newY += 2 * pivot.y;
                                break;
                        }

                        break;
                    }
                }
                ret[pt] = new (newX, newY);
            }
#if !GWSINGLETHREADED
        );
#endif
        if (H ^ V)
        {
            ret.Reverse(); // preserve ordering.
        }

        return ret;
    }

}