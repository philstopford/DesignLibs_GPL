using System;
using Clipper2Lib;
using utility;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static double distanceBetweenPoints(Point64 pt1, Point64 pt2)
    {
        return pDistanceBetweenPoints(pt1, pt2);
    }

    private static double pDistanceBetweenPoints(Point64 pt1, Point64 pt2)
    {
        return Math.Sqrt(Utils.myPow(pt1.X - pt2.X, 2) + Utils.myPow(pt1.Y - pt2.Y, 2));
    }
    public static double distanceBetweenPoints(PointD pt1, PointD pt2)
    {
        return pDistanceBetweenPoints(pt1, pt2);
    }

    private static double pDistanceBetweenPoints(PointD pt1, PointD pt2)
    {
        return Math.Sqrt(Utils.myPow(pt1.x - pt2.x, 2) + Utils.myPow(pt1.y - pt2.y, 2));
    }

    public static double distanceBetweenPoints(PointD pt1, Point64 pt2)
    {
        return pDistanceBetweenPoints(pt1, pt2);
    }

    private static double pDistanceBetweenPoints(PointD pt1, Point64 pt2)
    {
        return Math.Sqrt(Utils.myPow(pt1.x - pt2.X, 2) + Utils.myPow(pt1.y - pt2.Y, 2));
    }
}