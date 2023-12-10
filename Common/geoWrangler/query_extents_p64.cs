using Clipper2Lib;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static Point64 getExtents(Path64 source)
    {
        Path64 bounds = pGetBounds(source);
        Point64 extents = Point64_distanceBetweenPoints(bounds[1], bounds[0]);

        return extents;
    }
}