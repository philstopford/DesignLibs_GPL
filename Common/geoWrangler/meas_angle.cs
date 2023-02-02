using System;
using Clipper2Lib;
using utility;

namespace geoWrangler;

public static partial class GeoWrangler
{
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

    public static double angleBetweenPoints(PointD interSection_A, PointD interSection_B, PointD interSection_C, bool allowNegative = false)
    {
        return pAngleBetweenPoints(interSection_A, interSection_B, interSection_C, allowNegative);
    }

    private static double pAngleBetweenPoints(PointD interSection_A, PointD interSection_B, PointD interSection_C, bool allowNegative)
    {
        PointD cBVector = new(interSection_B.x - interSection_C.x, interSection_B.y - interSection_C.y);
        PointD cAVector = new(interSection_A.x - interSection_C.x, interSection_A.y - interSection_C.y);

        double xComponents = cBVector.x * cAVector.x;
        double yComponents = cBVector.y * cAVector.y;

        double scalarProduct = xComponents + yComponents;

        return angleCalc(cAVector.x, cAVector.y, cBVector.x, cBVector.y, scalarProduct, allowNegative);
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