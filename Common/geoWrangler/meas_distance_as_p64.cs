using Clipper2Lib;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static Point64 Point64_distanceBetweenPoints(Point64 pt1, Point64 pt2)
    {
        return pPoint64_distanceBetweenPoints(pt1, pt2);
    }

    private static Point64 pPoint64_distanceBetweenPoints(Point64 pt1, Point64 pt2)
    {
        long x_Distance = pt1.X - pt2.X;
        long y_Distance = pt1.Y - pt2.Y;
        return new (x_Distance, y_Distance);
    }
    
    public static PointD distanceBetweenPoints_point(Point64 pt1, Point64 pt2)
    {
        return pDistanceBetweenPoints_point(pt1, pt2);
    }

    private static PointD pDistanceBetweenPoints_point(Point64 pt1, Point64 pt2)
    {
        PointD ret = new(pt1.X - pt2.X, pt1.Y - pt2.Y);
        return ret;
    }
}