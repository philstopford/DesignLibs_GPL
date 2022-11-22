using Clipper2Lib;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static Point64 getExtents(Path64 source)
    {
        Path64 bounds = pGetBounds(source);
        Point64 extents = Point64_distanceBetweenPoints(bounds[0], bounds[2]);

        return extents;
    }
    
    public static Path64 getBounds(Path64 source)
    {
        return pGetBounds(source);
    }

    private static Path64 pGetBounds(Path64 source)
    {
        double minX = 0;
        double minY = 0;
        double maxX = 0;
        double maxY = 0;

        if (source.Count > 0)
        {
            try
            {
                minX = source[MinX(source)].X;
                maxX = source[MaxX(source)].X;
                minY = source[MinY(source)].Y;
                maxY = source[MaxY(source)].Y;
            }
            catch
            {
                // ignored
            }
        }

        return new () { new (minX, minY), new (maxX, maxY) };
    }
}