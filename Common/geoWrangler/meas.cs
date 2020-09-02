using ClipperLib;
using geoLib;
using System;
using utility;

namespace geoWrangler
{
    public static partial class GeoWrangler
    {
        public static IntPoint intPoint_distanceBetweenPoints(IntPoint pt1, IntPoint pt2)
        {
            return pIntPoint_distanceBetweenPoints(pt1, pt2);
        }

        static IntPoint pIntPoint_distanceBetweenPoints(IntPoint pt1, IntPoint pt2)
        {
            Int64 x_Distance = pt1.X - pt2.X;
            Int64 y_Distance = pt1.Y - pt2.Y;
            return new IntPoint(x_Distance, y_Distance);
        }

        public static double distanceBetweenPoints(IntPoint pt1, IntPoint pt2)
        {
            return pDistanceBetweenPoints(pt1, pt2);
        }

        static double pDistanceBetweenPoints(IntPoint pt1, IntPoint pt2)
        {
            return Math.Sqrt(Utils.myPow(pt1.X - pt2.X, 2) + Utils.myPow(pt1.Y - pt2.Y, 2));
        }

        public static double distanceBetweenPoints(GeoLibPointF pt1, GeoLibPointF pt2)
        {
            return pDistanceBetweenPoints(pt1, pt2);
        }

        static double pDistanceBetweenPoints(GeoLibPointF pt1, GeoLibPointF pt2)
        {
            return Math.Sqrt(Utils.myPow(pt1.X - pt2.X, 2) + Utils.myPow(pt1.Y - pt2.Y, 2));
        }

        public static double distanceBetweenPoints(GeoLibPoint pt1, GeoLibPoint pt2)
        {
            return pDistanceBetweenPoints(pt1, pt2);
        }

        static double pDistanceBetweenPoints(GeoLibPoint pt1, GeoLibPoint pt2)
        {
            return Math.Sqrt(Utils.myPow(pt1.X - pt2.X, 2) + Utils.myPow(pt1.Y - pt2.Y, 2));
        }


        public static GeoLibPointF distanceBetweenPoints_point(GeoLibPointF pt1, GeoLibPointF pt2)
        {
            return pDistanceBetweenPoints_point(pt1, pt2);
        }

        static GeoLibPointF pDistanceBetweenPoints_point(GeoLibPointF pt1, GeoLibPointF pt2)
        {
            GeoLibPointF ret = new GeoLibPointF(pt1.X - pt2.X, pt1.Y - pt2.Y);
            return ret;
        }

        public static GeoLibPointF distanceBetweenPoints_point(GeoLibPoint pt1, GeoLibPoint pt2)
        {
            return pDistanceBetweenPoints_point(pt1, pt2);
        }

        static GeoLibPointF pDistanceBetweenPoints_point(GeoLibPoint pt1, GeoLibPoint pt2)
        {
            GeoLibPointF ret = new GeoLibPointF(pt1.X - pt2.X, pt1.Y - pt2.Y);
            return ret;
        }

        public static double angleBetweenPoints(IntPoint interSection_A, IntPoint interSection_B, IntPoint interSection_C, bool allowNegative = false)
        {
            return pAngleBetweenPoints(interSection_A, interSection_B, interSection_C, allowNegative);
        }

        static double pAngleBetweenPoints(IntPoint interSection_A, IntPoint interSection_B, IntPoint interSection_C, bool allowNegative)
        {
            IntPoint cBVector = new IntPoint(interSection_B.X - interSection_C.X, interSection_B.Y - interSection_C.Y);
            IntPoint cAVector = new IntPoint(interSection_A.X - interSection_C.X, interSection_A.Y - interSection_C.Y);

            Int64 xComponents = cBVector.X * cAVector.X;
            Int64 yComponents = cBVector.Y * cAVector.Y;

            Int64 scalarProduct = xComponents + yComponents;

            return angleCalc(cAVector.X, cAVector.Y, cBVector.X, cBVector.Y, scalarProduct, allowNegative);
        }

        static double angleCalc(Int64 aVX, Int64 aVY, Int64 bVX, Int64 bVY, Int64 scalarProduct, bool allowNegative)
        {
            double cBMagnitude = (Math.Sqrt(Utils.myPow(bVX, 2) + Utils.myPow(bVY, 2)));
            double cAMagnitude = (Math.Sqrt(Utils.myPow(aVX, 2) + Utils.myPow(aVY, 2)));

            double theta = Utils.toDegrees(Math.Acos((scalarProduct) / (cBMagnitude * cAMagnitude)));

            if (!allowNegative)
            {
                // Avoid falling into a trap with negative angles.
                theta = Math.Abs(theta);
            }
            return theta;
        }

        public static double angleBetweenPoints(GeoLibPoint interSection_A, GeoLibPoint interSection_B, GeoLibPoint interSection_C, bool allowNegative = false)
        {
            return pAngleBetweenPoints(interSection_A, interSection_B, interSection_C, allowNegative);
        }

        static double pAngleBetweenPoints(GeoLibPoint interSection_A, GeoLibPoint interSection_B, GeoLibPoint interSection_C, bool allowNegative)
        {
            GeoLibPoint cBVector = new GeoLibPoint(interSection_B.X - interSection_C.X, interSection_B.Y - interSection_C.Y);
            GeoLibPoint cAVector = new GeoLibPoint(interSection_A.X - interSection_C.X, interSection_A.Y - interSection_C.Y);

            Int64 xComponents = cBVector.X * cAVector.X;
            Int64 yComponents = cBVector.Y * cAVector.Y;

            Int64 scalarProduct = xComponents + yComponents;

            return angleCalc(cAVector.X, cAVector.Y, cBVector.X, cBVector.Y, scalarProduct, allowNegative);
        }

        public static double angleBetweenPoints(GeoLibPointF interSection_A, GeoLibPointF interSection_B, GeoLibPointF interSection_C, bool allowNegative = false)
        {
            return pAngleBetweenPoints(interSection_A, interSection_B, interSection_C, allowNegative);
        }

        static double pAngleBetweenPoints(GeoLibPointF interSection_A, GeoLibPointF interSection_B, GeoLibPointF interSection_C, bool allowNegative)
        {
            GeoLibPointF cBVector = new GeoLibPointF(interSection_B.X - interSection_C.X, interSection_B.Y - interSection_C.Y);
            GeoLibPointF cAVector = new GeoLibPointF(interSection_A.X - interSection_C.X, interSection_A.Y - interSection_C.Y);

            double xComponents = cBVector.X * cAVector.X;
            double yComponents = cBVector.Y * cAVector.Y;

            double scalarProduct = xComponents + yComponents;

            return angleCalc(cAVector.X, cAVector.Y, cBVector.X, cBVector.Y, scalarProduct, allowNegative);
        }

        static double angleCalc(double aVX, double aVY, double bVX, double bVY, double scalarProduct, bool allowNegative)
        {
            double cBMagnitude = (Math.Sqrt(Utils.myPow(bVX, 2) + Utils.myPow(bVY, 2)));
            double cAMagnitude = (Math.Sqrt(Utils.myPow(aVX, 2) + Utils.myPow(aVY, 2)));

            double theta = Math.Abs(Utils.toDegrees(Math.Acos((scalarProduct) / (cBMagnitude * cAMagnitude)))); // Avoid falling into a trap with negative angles.

            if (!allowNegative)
            {
                // Avoid falling into a trap with negative angles.
                theta = Math.Abs(theta);
            }

            return theta;
        }

    }
}
