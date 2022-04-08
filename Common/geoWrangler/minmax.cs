using Clipper2Lib;
using geoLib;
using System.Collections.Generic;
using System.Linq;

namespace geoWrangler;

using Path = List<Point64>;

public static partial class GeoWrangler
{
    public static int MinX(Path sourcePath)
    {
        return pMinX(sourcePath);
    }

    private static int pMinX(Path sourcePath)
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

    public static int MaxX(Path sourcePath)
    {
        return pMaxX(sourcePath);
    }

    private static int pMaxX(Path sourcePath)
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

    public static int MinY(Path sourcePath)
    {
        return pMinY(sourcePath);
    }

    private static int pMinY(Path sourcePath)
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

    public static int MaxY(Path sourcePath)
    {
        return pMaxY(sourcePath);
    }

    private static int pMaxY(Path sourcePath)
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

    public static int MinX(List<GeoLibPoint> iPoints)
    {
        return pMinX(iPoints.ToArray());
    }
    
    public static int MinX(GeoLibPoint[] iPoints)
    {
        return pMinX(iPoints);
    }

    private static int pMinX(GeoLibPoint[] iPoints)
    {
        long min = iPoints[0].X;
        int minIndex = 0;

        for (int i = 1; i < iPoints.Length; ++i)
        {
            if (iPoints[i].X >= min)
            {
                continue;
            }

            min = iPoints[i].X;
            minIndex = i;
        }

        return minIndex;
    }

    public static int MinY(GeoLibPoint[] iPoints)
    {
        return pMinY(iPoints);
    }

    private static int pMinY(GeoLibPoint[] iPoints)
    {
        long min = iPoints[0].Y;
        int minIndex = 0;

        for (int i = 1; i < iPoints.Length; ++i)
        {
            if (iPoints[i].Y >= min)
            {
                continue;
            }

            min = iPoints[i].Y;
            minIndex = i;
        }

        return minIndex;
    }

    public static int MaxX(GeoLibPoint[] iPoints)
    {
        return pMaxX(iPoints);
    }

    private static int pMaxX(GeoLibPoint[] iPoints)
    {
        long max = iPoints[0].X;
        int maxIndex = 0;

        for (int i = 1; i < iPoints.Length; ++i)
        {
            if (iPoints[i].X <= max)
            {
                continue;
            }

            max = iPoints[i].X;
            maxIndex = i;
        }

        return maxIndex;
    }

    public static int MaxY(GeoLibPoint[] iPoints)
    {
        return pMaxY(iPoints);
    }

    private static int pMaxY(GeoLibPoint[] iPoints)
    {
        long max = iPoints[0].Y;
        int maxIndex = 0;

        for (int i = 1; i < iPoints.Length; ++i)
        {
            if (iPoints[i].Y <= max)
            {
                continue;
            }

            max = iPoints[i].Y;
            maxIndex = i;
        }

        return maxIndex;
    }

    public static int MinX(GeoLibPointF[] iPoints)
    {
        return pMinX(iPoints);
    }

    private static int pMinX(GeoLibPointF[] iPoints)
    {
        double min = iPoints[0].X;
        int minIndex = 0;

        for (int i = 1; i < iPoints.Length; ++i)
        {
            if (!(iPoints[i].X < min))
            {
                continue;
            }

            min = iPoints[i].X;
            minIndex = i;
        }

        return minIndex;
    }

    public static int MinY(GeoLibPointF[] iPoints)
    {
        return pMinY(iPoints);
    }

    private static int pMinY(GeoLibPointF[] iPoints)
    {
        double min = iPoints[0].Y;
        int minIndex = 0;

        for (int i = 1; i < iPoints.Length; ++i)
        {
            if (!(iPoints[i].Y < min))
            {
                continue;
            }

            min = iPoints[i].Y;
            minIndex = i;
        }

        return minIndex;
    }

    public static int MaxX(GeoLibPointF[] iPoints)
    {
        return pMaxX(iPoints);
    }

    private static int pMaxX(GeoLibPointF[] iPoints)
    {
        double max = iPoints[0].X;
        int maxIndex = 0;

        for (int i = 1; i < iPoints.Length; ++i)
        {
            if (!(iPoints[i].X > max))
            {
                continue;
            }

            max = iPoints[i].X;
            maxIndex = i;
        }

        return maxIndex;
    }

    public static int MaxY(GeoLibPointF[] iPoints)
    {
        return pMaxY(iPoints);
    }

    private static int pMaxY(GeoLibPointF[] iPoints)
    {
        double max = iPoints[0].Y;
        int maxIndex = 0;

        for (int i = 1; i < iPoints.Length; ++i)
        {
            if (!(iPoints[i].Y > max))
            {
                continue;
            }

            max = iPoints[i].Y;
            maxIndex = i;
        }

        return maxIndex;
    }

    public static GeoLibPoint getMinimumPoint(GeoLibPoint[] iPoints)
    {
        return pGetMinimumPoint(iPoints);
    }

    private static GeoLibPoint pGetMinimumPoint(GeoLibPoint[] iPoints)
    {
        int x = iPoints.Min(p => p.X);
        int y = iPoints.Min(p => p.Y);

        return new GeoLibPoint(x, y);
    }

    public static GeoLibPoint getMaximumPoint(GeoLibPoint[] iPoints)
    {
        return pGetMaximumPoint(iPoints);
    }

    private static GeoLibPoint pGetMaximumPoint(GeoLibPoint[] iPoints)
    {
        int x = iPoints.Max(p => p.X);
        int y = iPoints.Max(p => p.Y);

        return new GeoLibPoint(x, y);
    }

    public static GeoLibPointF getMinimumPoint(GeoLibPointF[] iPoints)
    {
        return pGetMinimumPoint(iPoints);
    }

    private static GeoLibPointF pGetMinimumPoint(GeoLibPointF[] iPoints)
    {
        double x = iPoints.Min(p => p.X);
        double y = iPoints.Min(p => p.Y);

        return new GeoLibPointF(x, y);
    }

    public static GeoLibPointF getMaximumPoint(GeoLibPointF[] iPoints)
    {
        return pGetMaximumPoint(iPoints);
    }

    private static GeoLibPointF pGetMaximumPoint(GeoLibPointF[] iPoints)
    {
        double x = iPoints.Max(p => p.X);
        double y = iPoints.Max(p => p.Y);

        return new GeoLibPointF(x, y);
    }
}