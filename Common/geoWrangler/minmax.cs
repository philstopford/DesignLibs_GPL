using Clipper2Lib;
using geoLib;
using System.Collections.Generic;
using System.Linq;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static int MinX(Path64 sourcePath)
    {
        return pMinX(sourcePath);
    }

    private static int pMinX(Path64 sourcePath)
    {
        long min = sourcePath[0].X;
        int minIndex = 0;

        for (int i = 1; i < sourcePath.Count; ++i)
        {
            if (sourcePath[i].X >= min)
            {
                continue;
            }

            min = sourcePath[i].X;
            minIndex = i;
        }

        return minIndex;
    }

    public static int MaxX(Path64 sourcePath)
    {
        return pMaxX(sourcePath);
    }

    private static int pMaxX(Path64 sourcePath)
    {
        long max = sourcePath[0].X;
        int maxIndex = 0;

        for (int i = 1; i < sourcePath.Count; ++i)
        {
            if (sourcePath[i].X <= max)
            {
                continue;
            }

            max = sourcePath[i].X;
            maxIndex = i;
        }

        return maxIndex;
    }

    public static int MinY(Path64 sourcePath)
    {
        return pMinY(sourcePath);
    }

    private static int pMinY(Path64 sourcePath)
    {
        long min = sourcePath[0].Y;
        int minIndex = 0;

        for (int i = 1; i < sourcePath.Count; ++i)
        {
            if (sourcePath[i].Y >= min)
            {
                continue;
            }

            min = sourcePath[i].Y;
            minIndex = i;
        }

        return minIndex;
    }

    public static int MaxY(Path64 sourcePath)
    {
        return pMaxY(sourcePath);
    }

    private static int pMaxY(Path64 sourcePath)
    {
        long max = sourcePath[0].Y;
        int maxIndex = 0;

        for (int i = 1; i < sourcePath.Count; ++i)
        {
            if (sourcePath[i].Y <= max)
            {
                continue;
            }

            max = sourcePath[i].Y;
            maxIndex = i;
        }

        return maxIndex;
    }

    public static int MinX(PathD iPoints)
    {
        return pMinX(iPoints);
    }

    private static int pMinX(PathD iPoints)
    {
        double min = iPoints[0].x;
        int minIndex = 0;

        for (int i = 1; i < iPoints.Capacity; ++i)
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

        for (int i = 1; i < iPoints.Capacity; ++i)
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

        for (int i = 1; i < iPoints.Capacity; ++i)
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

        for (int i = 1; i < iPoints.Capacity; ++i)
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

    public static Point64 getMinimumPoint(Path64 iPoints)
    {
        return pGetMinimumPoint(iPoints);
    }

    private static Point64 pGetMinimumPoint(Path64 iPoints)
    {
        int x = (int)iPoints.Min(p => p.X);
        int y = (int)iPoints.Min(p => p.Y);

        return new (x, y);
    }

    public static Point64 getMaximumPoint(Path64 iPoints)
    {
        return pGetMaximumPoint(iPoints);
    }

    private static Point64 pGetMaximumPoint(Path64 iPoints)
    {
        int x = (int)iPoints.Max(p => p.X);
        int y = (int)iPoints.Max(p => p.Y);

        return new (x, y);
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