using System;
using Clipper2Lib;

namespace geoWrangler;

public static partial class GeoWrangler{
    private static double CrossProduct(Point64 pt1, Point64 pt2, Point64 pt3)
    {
        // typecast to double to avoid potential int overflow
        return ((double) (pt2.X - pt1.X) * (pt3.Y - pt2.Y) -
                (double) (pt2.Y - pt1.Y) * (pt3.X - pt2.X));
    }

    private static double DotProduct(Point64 pt1, Point64 pt2, Point64 pt3)
    {
        // typecast to double to avoid potential int overflow
        return ((double) (pt2.X - pt1.X) * (pt3.X - pt2.X) +
                (double) (pt2.Y - pt1.Y) * (pt3.Y - pt2.Y));
    }

    private static double CrossProduct(PointD vec1, PointD vec2)
    {
        return (vec1.y * vec2.x - vec2.y * vec1.x);
    }

    private static double DotProduct(PointD vec1, PointD vec2)
    {
        return (vec1.x * vec2.x + vec1.y * vec2.y);
    }

    private static double CrossProduct(PointD pt1, PointD pt2, PointD pt3)
    {
        // typecast to double to avoid potential int overflow
        return ((pt2.x - pt1.x) * (pt3.y - pt2.y) -
                (pt2.y - pt1.y) * (pt3.x - pt2.x));
    }

    private static double DotProduct(PointD pt1, PointD pt2, PointD pt3)
    {
        // typecast to double to avoid potential int overflow
        return ((pt2.x - pt1.x) * (pt3.x - pt2.x) +
                (pt2.y - pt1.y) * (pt3.y - pt2.y));
    }

    private static double varOffset(Path64 path, PathD norms, int pt, int prevPt)
    {
        // sin(A) < 0: right turning
        // cos(A) < 0: change in angle is more than 90 degree
        int nextPt = pt + 1;
        if (nextPt == path.Count)
        {
            nextPt = 0;
        }
        double sinA = CrossProduct(norms[nextPt], norms[pt], norms[prevPt]);
        double cosA = DotProduct(norms[nextPt], norms[pt], norms[prevPt]);

        return Math.Pow(norms[pt].y + norms[prevPt].y, 2) * 10;
    }
}