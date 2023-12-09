using Clipper2Lib;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static PointD getExtents(PathD source)
    {
        PathD bounds = pGetBounds(source);
        PointD extents = pDistanceBetweenPoints_point(bounds[1], bounds[0]);

        return extents;
    }
    
    public static PathD getBounds(PathD source)
    {
        return pGetBounds(source);
    }

    private static PathD pGetBounds(PathD source)
    {
        double minX = 0;
        double minY = 0;
        double maxX = 0;
        double maxY = 0;

        if (source.Count > 0)
        {
            try
            {
                minX = source[MinX(source)].x;
                maxX = source[MaxX(source)].x;
                minY = source[MinY(source)].y;
                maxY = source[MaxY(source)].y;
            }
            catch
            {
                // ignored
            }
        }

        return new () { new (minX, minY), new (maxX, maxY) };
    }
}