using Clipper2Lib;
using System;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static Point64 midPoint(Paths64 source)
    {
        return pMidPoint(source);
    }

    private static Point64 pMidPoint(Paths64 source)
    {
        Path64 bounds = getBounds(source[0]);

        double minX = bounds[0].X;
        double minY = bounds[0].Y;
        double maxX = bounds[1].X;
        double maxY = bounds[1].Y;

        for (int i = 1; i < source.Count; i++)
        {
            Path64 tmp = getBounds(source[i]);
            minX = Math.Min(minX, tmp[0].X);
            minY = Math.Min(minY, tmp[0].Y);
            maxX = Math.Max(maxX, tmp[1].X);
            maxY = Math.Max(maxY, tmp[1].Y);
        }

        double avX = minX + (maxX - minX) / 2.0f;
        double avY = minY + (maxY - minY) / 2.0f;

        return new (avX, avY);
    }

    public static PointD midPoint(Path64 source)
    {
        return pMidPoint(source);
    }

    private static PointD pMidPoint(Path64 source)
    {
        Path64 bounds = getBounds(source);

        double minX = bounds[0].X;
        double minY = bounds[0].Y;
        double maxX = bounds[1].X;
        double maxY = bounds[1].Y;

        double avX = minX + (maxX - minX) / 2.0f;
        double avY = minY + (maxY - minY) / 2.0f;

        return new (avX, avY);
    }
}