using System;
using System.Runtime.CompilerServices;
using Clipper2Lib;
using utility;

namespace geoWrangler;

public static partial class GeoWrangler
{
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double distanceBetweenPoints(Point64 pt1, Point64 pt2)
    {
        return pDistanceBetweenPoints(pt1, pt2);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double pDistanceBetweenPoints(Point64 pt1, Point64 pt2)
    {
        // Optimized: avoid Utils.myPow for simple case of squaring
        double dx = pt1.X - pt2.X;
        double dy = pt1.Y - pt2.Y;
        return Math.Sqrt(dx * dx + dy * dy);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double distanceBetweenPoints(PointD pt1, PointD pt2)
    {
        return pDistanceBetweenPoints(pt1, pt2);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double pDistanceBetweenPoints(PointD pt1, PointD pt2)
    {
        // Optimized: avoid Utils.myPow for simple case of squaring
        double dx = pt1.x - pt2.x;
        double dy = pt1.y - pt2.y;
        return Math.Sqrt(dx * dx + dy * dy);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double distanceBetweenPoints(PointD pt1, Point64 pt2)
    {
        return pDistanceBetweenPoints(pt1, pt2);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double pDistanceBetweenPoints(PointD pt1, Point64 pt2)
    {
        // Optimized: avoid Utils.myPow for simple case of squaring
        double dx = pt1.x - pt2.X;
        double dy = pt1.y - pt2.Y;
        return Math.Sqrt(dx * dx + dy * dy);
    }

    /// <summary>
    /// Optimized distance squared calculation - avoids expensive sqrt operation
    /// Use when you only need to compare distances
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double distanceSquaredBetweenPoints(Point64 pt1, Point64 pt2)
    {
        double dx = pt1.X - pt2.X;
        double dy = pt1.Y - pt2.Y;
        return dx * dx + dy * dy;
    }

    /// <summary>
    /// Optimized distance squared calculation - avoids expensive sqrt operation
    /// Use when you only need to compare distances
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double distanceSquaredBetweenPoints(PointD pt1, PointD pt2)
    {
        double dx = pt1.x - pt2.x;
        double dy = pt1.y - pt2.y;
        return dx * dx + dy * dy;
    }
}