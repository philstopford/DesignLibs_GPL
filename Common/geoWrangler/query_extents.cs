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
}