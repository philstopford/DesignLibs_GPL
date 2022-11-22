using System;
using System.Threading.Tasks;
using Clipper2Lib;
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

        int pointListLength = pointList.Count;
        Path64 returnList = Helper.initedPath64(pointListLength);

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
}