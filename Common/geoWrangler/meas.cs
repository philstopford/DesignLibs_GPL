using ClipperLib2;
using geoLib;
using System;
using utility;

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
        return new Point64(x_Distance, y_Distance);
    }

    public static double distanceBetweenPoints(Point64 pt1, Point64 pt2)
    {
        return pDistanceBetweenPoints(pt1, pt2);
    }

    private static double pDistanceBetweenPoints(Point64 pt1, Point64 pt2)
    {
        return Math.Sqrt(Utils.myPow(pt1.X - pt2.X, 2) + Utils.myPow(pt1.Y - pt2.Y, 2));
    }

    public static double distanceBetweenPoints(GeoLibPointF pt1, GeoLibPointF pt2)
    {
        return pDistanceBetweenPoints(pt1, pt2);
    }

    private static double pDistanceBetweenPoints(GeoLibPointF pt1, GeoLibPointF pt2)
    {
        return Math.Sqrt(Utils.myPow(pt1.X - pt2.X, 2) + Utils.myPow(pt1.Y - pt2.Y, 2));
    }

    public static double distanceBetweenPoints(GeoLibPoint pt1, GeoLibPoint pt2)
    {
        return pDistanceBetweenPoints(pt1, pt2);
    }

    private static double pDistanceBetweenPoints(GeoLibPoint pt1, GeoLibPoint pt2)
    {
        return Math.Sqrt(Utils.myPow(pt1.X - pt2.X, 2) + Utils.myPow(pt1.Y - pt2.Y, 2));
    }

    public static GeoLibPointF distanceBetweenPoints_point(GeoLibPointF pt1, GeoLibPointF pt2)
    {
        return pDistanceBetweenPoints_point(pt1, pt2);
    }

    private static GeoLibPointF pDistanceBetweenPoints_point(GeoLibPointF pt1, GeoLibPointF pt2)
    {
        GeoLibPointF ret = new(pt1.X - pt2.X, pt1.Y - pt2.Y);
        return ret;
    }

    public static GeoLibPointF distanceBetweenPoints_point(GeoLibPoint pt1, GeoLibPoint pt2)
    {
        return pDistanceBetweenPoints_point(pt1, pt2);
    }

    private static GeoLibPointF pDistanceBetweenPoints_point(GeoLibPoint pt1, GeoLibPoint pt2)
    {
        GeoLibPointF ret = new(pt1.X - pt2.X, pt1.Y - pt2.Y);
        return ret;
    }

    public static double angleBetweenPoints(Point64 interSection_A, Point64 interSection_B, Point64 interSection_C, bool allowNegative = false)
    {
        return pAngleBetweenPoints(interSection_A, interSection_B, interSection_C, allowNegative);
    }

    private static double pAngleBetweenPoints(Point64 interSection_A, Point64 interSection_B, Point64 interSection_C, bool allowNegative)
    {
        Point64 cBVector = new(interSection_B.X - interSection_C.X, interSection_B.Y - interSection_C.Y);
        Point64 cAVector = new(interSection_A.X - interSection_C.X, interSection_A.Y - interSection_C.Y);

        long xComponents = cBVector.X * cAVector.X;
        long yComponents = cBVector.Y * cAVector.Y;

        long scalarProduct = xComponents + yComponents;

        return angleCalc(cAVector.X, cAVector.Y, cBVector.X, cBVector.Y, scalarProduct, allowNegative);
    }

    private static double angleCalc(long aVX, long aVY, long bVX, long bVY, long scalarProduct, bool allowNegative)
    {
        double cBMagnitude = Math.Sqrt(Utils.myPow(bVX, 2) + Utils.myPow(bVY, 2));
        double cAMagnitude = Math.Sqrt(Utils.myPow(aVX, 2) + Utils.myPow(aVY, 2));

        double theta = Utils.toDegrees(Math.Acos(scalarProduct / (cBMagnitude * cAMagnitude)));

        switch (allowNegative)
        {
            case false:
                // Avoid falling into a trap with negative angles.
                theta = Math.Abs(theta);
                break;
        }
        return theta;
    }

    public static double angleBetweenPoints(GeoLibPoint interSection_A, GeoLibPoint interSection_B, GeoLibPoint interSection_C, bool allowNegative = false)
    {
        return pAngleBetweenPoints(interSection_A, interSection_B, interSection_C, allowNegative);
    }

    private static double pAngleBetweenPoints(GeoLibPoint interSection_A, GeoLibPoint interSection_B, GeoLibPoint interSection_C, bool allowNegative)
    {
        GeoLibPoint cBVector = new(interSection_B.X - interSection_C.X, interSection_B.Y - interSection_C.Y);
        GeoLibPoint cAVector = new(interSection_A.X - interSection_C.X, interSection_A.Y - interSection_C.Y);

        long xComponents = cBVector.X * cAVector.X;
        long yComponents = cBVector.Y * cAVector.Y;

        long scalarProduct = xComponents + yComponents;

        return angleCalc(cAVector.X, cAVector.Y, cBVector.X, cBVector.Y, scalarProduct, allowNegative);
    }

    public static double angleBetweenPoints(GeoLibPointF interSection_A, GeoLibPointF interSection_B, GeoLibPointF interSection_C, bool allowNegative = false)
    {
        return pAngleBetweenPoints(interSection_A, interSection_B, interSection_C, allowNegative);
    }

    private static double pAngleBetweenPoints(GeoLibPointF interSection_A, GeoLibPointF interSection_B, GeoLibPointF interSection_C, bool allowNegative)
    {
        GeoLibPointF cBVector = new(interSection_B.X - interSection_C.X, interSection_B.Y - interSection_C.Y);
        GeoLibPointF cAVector = new(interSection_A.X - interSection_C.X, interSection_A.Y - interSection_C.Y);

        double xComponents = cBVector.X * cAVector.X;
        double yComponents = cBVector.Y * cAVector.Y;

        double scalarProduct = xComponents + yComponents;

        return angleCalc(cAVector.X, cAVector.Y, cBVector.X, cBVector.Y, scalarProduct, allowNegative);
    }

    private static double angleCalc(double aVX, double aVY, double bVX, double bVY, double scalarProduct, bool allowNegative)
    {
        double cBMagnitude = Math.Sqrt(Utils.myPow(bVX, 2) + Utils.myPow(bVY, 2));
        double cAMagnitude = Math.Sqrt(Utils.myPow(aVX, 2) + Utils.myPow(aVY, 2));

        double theta = Math.Abs(Utils.toDegrees(Math.Acos(scalarProduct / (cBMagnitude * cAMagnitude)))); // Avoid falling into a trap with negative angles.

        switch (allowNegative)
        {
            case false:
                // Avoid falling into a trap with negative angles.
                theta = Math.Abs(theta);
                break;
        }

        return theta;
    }

}