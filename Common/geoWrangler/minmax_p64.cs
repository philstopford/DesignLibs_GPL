using Clipper2Lib;
using System.Linq;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static int MinX(Path64 iPoints)
    {
        return pMinX(iPoints);
    }

    private static int pMinX(Path64 iPoints)
    {
        double min = iPoints[0].X;
        int minIndex = 0;

        for (int i = 1; i < iPoints.Count; ++i)
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

    public static int MinY(Path64 iPoints)
    {
        return pMinY(iPoints);
    }

    private static int pMinY(Path64 iPoints)
    {
        double min = iPoints[0].Y;
        int minIndex = 0;

        for (int i = 1; i < iPoints.Count; ++i)
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

    public static int MaxX(Path64 iPoints)
    {
        return pMaxX(iPoints);
    }

    private static int pMaxX(Path64 iPoints)
    {
        double max = iPoints[0].X;
        int maxIndex = 0;

        for (int i = 1; i < iPoints.Count; ++i)
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

    public static int MaxY(Path64 iPoints)
    {
        return pMaxY(iPoints);
    }

    private static int pMaxY(Path64 iPoints)
    {
        double max = iPoints[0].Y;
        int maxIndex = 0;

        for (int i = 1; i < iPoints.Count; ++i)
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
    
    public static Point64 getMinimumPoint(Path64 iPoints)
    {
        return pGetMinimumPoint(iPoints);
    }

    private static Point64 pGetMinimumPoint(Path64 iPoints)
    {
        double x = iPoints.Min(p => p.X);
        double y = iPoints.Min(p => p.Y);

        return new Point64(x, y);
    }

    public static Point64 getMaximumPoint(Path64 iPoints)
    {
        return pGetMaximumPoint(iPoints);
    }

    private static Point64 pGetMaximumPoint(Path64 iPoints)
    {
        double x = iPoints.Max(p => p.X);
        double y = iPoints.Max(p => p.Y);

        return new Point64(x, y);
    }
}