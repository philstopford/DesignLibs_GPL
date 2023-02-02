using Clipper2Lib;
using System;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static PointD midPoint(PathsD source)
    {
        return pMidPoint(source);
    }

    private static PointD pMidPoint(PathsD source)
    {
        PathD bounds = getBounds(source[0]);

        double minX = bounds[0].x;
        double minY = bounds[0].y;
        double maxX = bounds[1].x;
        double maxY = bounds[1].y;

        for (int i = 1; i < source.Count; i++)
        {
            PathD tmp = getBounds(source[i]);
            minX = Math.Min(minX, tmp[0].x);
            minY = Math.Min(minY, tmp[0].y);
            maxX = Math.Max(maxX, tmp[1].x);
            maxY = Math.Max(maxY, tmp[1].y);
        }

        double avX = minX + (maxX - minX) / 2.0f;
        double avY = minY + (maxY - minY) / 2.0f;

        return new (avX, avY);
    }

    public static PointD midPoint(PathD source)
    {
        return pMidPoint(source);
    }

    private static PointD pMidPoint(PathD source)
    {
        PathD bounds = getBounds(source);

        double minX = bounds[0].x;
        double minY = bounds[0].y;
        double maxX = bounds[1].x;
        double maxY = bounds[1].y;

        double avX = minX + (maxX - minX) / 2.0f;
        double avY = minY + (maxY - minY) / 2.0f;

        return new (avX, avY);
    }
}