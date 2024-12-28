using Clipper2Lib;
using System;
using System.Threading.Tasks;
using utility;

namespace geoWrangler;

public static partial class GeoWrangler
{
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
                return new PathD(pointList);
        }

        int pointListLength = pointList.Count;
        PathD returnList = Helper.initedPathD(pointListLength);

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
                return new PointD(point);
        }

        double angle = Utils.toRadians(angleDegree);
        double x = pivot.x + ((point.x - pivot.x) * Math.Cos(angle) -
                              (point.y - pivot.y) * Math.Sin(angle));
        double y = pivot.y + ((point.x - pivot.x) * Math.Sin(angle) +
                              (point.y - pivot.y) * Math.Cos(angle));

        return new PointD(x, y);
    }
    
    public static PathD flip(bool H, bool V, bool alignX, bool alignY, PointD pivot, PathD source)
    {
        return pFlip(H, V, alignX, alignY, pivot, source);
    }
    
    private static PathD pFlip(bool H, bool V, bool alignX, bool alignY, PointD pivot_, PathD source)
    {
        PointD pivot = new(pivot_);
        if (double.IsNaN(pivot.x) || double.IsNaN(pivot.y))
        {
            pivot = midPoint(source);
        }
        int sLength = source.Count;
        PathD ret = Helper.initedPathD(sLength);
#if !GWSINGLETHREADED
        Parallel.For(0, sLength, pt =>
#else
            for (Int32 pt = 0; pt < sLength; pt++)
#endif
            {
                double newX = source[pt].x;
                if (H)
                {
                    newX = -newX;
                    if (alignX)
                    {
                        newX += 2 * pivot.x;
                    }
                }

                double newY = source[pt].y;
                if (V)
                {
                    newY = -newY;
                    if (alignY)
                    {
                        newY += 2 * pivot.y;
                    }
                }

                ret[pt] = new PointD(newX, newY);
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