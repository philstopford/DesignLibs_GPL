using Clipper2Lib;
using System;
using utility;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static PointD PointD_distanceBetweenPoints(PointD pt1, PointD pt2)
    {
        return pPointD_distanceBetweenPoints(pt1, pt2);
    }

    private static PointD pPointD_distanceBetweenPoints(PointD pt1, PointD pt2)
    {
        double x_Distance = pt1.x - pt2.x;
        double y_Distance = pt1.y - pt2.y;
        return new (x_Distance, y_Distance);
    }
    
    public static PointD distanceBetweenPoints_point(PointD pt1, PointD pt2)
    {
        return pDistanceBetweenPoints_point(pt1, pt2);
    }

    private static PointD pDistanceBetweenPoints_point(PointD pt1, PointD pt2)
    {
        PointD ret = new(pt1.x - pt2.x, pt1.y - pt2.y);
        return ret;
    }
}