using Clipper2Lib;
using System.Linq;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static int MinX(PathD iPoints)
    {
        return pMinX(iPoints);
    }

    private static int pMinX(PathD iPoints)
    {
        double min = iPoints[0].x;
        int minIndex = 0;

        for (int i = 1; i < iPoints.Count; ++i)
        {
            if (!(iPoints[i].x < min))
            {
                continue;
            }

            min = iPoints[i].x;
            minIndex = i;
        }

        return minIndex;
    }

    public static int MinY(PathD iPoints)
    {
        return pMinY(iPoints);
    }

    private static int pMinY(PathD iPoints)
    {
        double min = iPoints[0].y;
        int minIndex = 0;

        for (int i = 1; i < iPoints.Count; ++i)
        {
            if (!(iPoints[i].y < min))
            {
                continue;
            }

            min = iPoints[i].y;
            minIndex = i;
        }

        return minIndex;
    }

    public static int MaxX(PathD iPoints)
    {
        return pMaxX(iPoints);
    }

    private static int pMaxX(PathD iPoints)
    {
        double max = iPoints[0].x;
        int maxIndex = 0;

        for (int i = 1; i < iPoints.Count; ++i)
        {
            if (!(iPoints[i].x > max))
            {
                continue;
            }

            max = iPoints[i].x;
            maxIndex = i;
        }

        return maxIndex;
    }

    public static int MaxY(PathD iPoints)
    {
        return pMaxY(iPoints);
    }

    private static int pMaxY(PathD iPoints)
    {
        double max = iPoints[0].y;
        int maxIndex = 0;

        for (int i = 1; i < iPoints.Count; ++i)
        {
            if (!(iPoints[i].y > max))
            {
                continue;
            }

            max = iPoints[i].y;
            maxIndex = i;
        }

        return maxIndex;
    }
    
    public static PointD getMinimumPoint(PathD iPoints)
    {
        return pGetMinimumPoint(iPoints);
    }

    private static PointD pGetMinimumPoint(PathD iPoints)
    {
        double x = iPoints.Min(p => p.x);
        double y = iPoints.Min(p => p.y);

        return new PointD(x, y);
    }

    public static PointD getMaximumPoint(PathD iPoints)
    {
        return pGetMaximumPoint(iPoints);
    }

    private static PointD pGetMaximumPoint(PathD iPoints)
    {
        double x = iPoints.Max(p => p.x);
        double y = iPoints.Max(p => p.y);

        return new PointD(x, y);
    }
}