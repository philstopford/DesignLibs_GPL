using ClipperLib;
using geoLib;
using System;
using System.Collections.Generic;

namespace geoWrangler
{
    using Path = List<IntPoint>;
    using Paths = List<List<IntPoint>>;

    public static partial class GeoWrangler
    {
        public static GeoLibPointF midPoint(List<GeoLibPoint[]> source)
        {
            return pMidPoint(source);
        }
        static GeoLibPointF pMidPoint(List<GeoLibPoint[]> source)
        {
            GeoLibPoint[] bounds = getBounds(source[0]);

            double minX = bounds[0].X;
            double minY = bounds[0].Y;
            double maxX = bounds[1].X;
            double maxY = bounds[1].Y;

            for (int i = 1; i < source.Count; i++)
            {
                GeoLibPoint[] tmp = getBounds(source[i]);
                minX = Math.Min(minX, tmp[0].X);
                minY = Math.Min(minY, tmp[0].Y);
                maxX = Math.Max(maxX, tmp[1].X);
                maxY = Math.Max(maxY, tmp[1].Y);
            }

            double avX = minX + (maxX - minX) / 2.0f;
            double avY = minY + (maxY - minY) / 2.0f;

            return new GeoLibPointF(avX, avY);
        }

        public static GeoLibPointF midPoint(GeoLibPoint[] source)
        {
            return pMidPoint(source);
        }

        static GeoLibPointF pMidPoint(GeoLibPoint[] source)
        {
            GeoLibPoint[] bounds = getBounds(source);

            double minX = bounds[0].X;
            double minY = bounds[0].Y;
            double maxX = bounds[1].X;
            double maxY = bounds[1].Y;

            double avX = minX + (maxX - minX) / 2.0f;
            double avY = minY + (maxY - minY) / 2.0f;

            return new GeoLibPointF(avX, avY);
        }

        public static GeoLibPointF midPoint(List<GeoLibPointF[]> source)
        {
            return pMidPoint(source);
        }
        static GeoLibPointF pMidPoint(List<GeoLibPointF[]> source)
        {
            GeoLibPointF[] bounds = getBounds(source[0]);

            double minX = bounds[0].X;
            double minY = bounds[0].Y;
            double maxX = bounds[1].X;
            double maxY = bounds[1].Y;

            for (int i = 1; i < source.Count; i++)
            {
                GeoLibPointF[] tmp = getBounds(source[i]);
                minX = Math.Min(minX, tmp[0].X);
                minY = Math.Min(minY, tmp[0].Y);
                maxX = Math.Max(maxX, tmp[1].X);
                maxY = Math.Max(maxY, tmp[1].Y);
            }

            double avX = minX + (maxX - minX) / 2.0f;
            double avY = minY + (maxY - minY) / 2.0f;

            return new GeoLibPointF(avX, avY);
        }

        public static GeoLibPointF midPoint(GeoLibPointF[] source)
        {
            return pMidPoint(source);
        }

        static GeoLibPointF pMidPoint(GeoLibPointF[] source)
        {
            GeoLibPointF[] bounds = getBounds(source);

            double minX = bounds[0].X;
            double minY = bounds[0].Y;
            double maxX = bounds[1].X;
            double maxY = bounds[1].Y;

            double avX = minX + (maxX - minX) / 2.0f;
            double avY = minY + (maxY - minY) / 2.0f;

            return new GeoLibPointF(avX, avY);
        }

        public static GeoLibPointF getExtents(GeoLibPointF[] source)
        {
            GeoLibPointF[] bounds = pGetBounds(source);
            GeoLibPointF extents = pDistanceBetweenPoints_point(source[0], source[2]);

            return extents;
        }

        public static GeoLibPoint[] getBounds(GeoLibPoint[] source)
        {
            return pGetBounds(source);
        }

        static GeoLibPoint[] pGetBounds(GeoLibPoint[] source)
        {
            double minX = 0;
            double minY = 0;
            double maxX = 0;
            double maxY = 0;

            try
            {
                minX = source[MinX(source)].X;
                maxX = source[MaxX(source)].X;
                minY = source[MinY(source)].Y;
                maxY = source[MaxY(source)].Y;
            }
            catch (Exception)
            {

            }

            return new GeoLibPoint[] { new GeoLibPoint(minX, minY), new GeoLibPoint(maxX, maxY) };
        }

        public static GeoLibPointF[] getBounds(GeoLibPointF[] source)
        {
            return pGetBounds(source);
        }

        static GeoLibPointF[] pGetBounds(GeoLibPointF[] source)
        {
            double minX = 0;
            double minY = 0;
            double maxX = 0;
            double maxY = 0;

            try
            {
                minX = source[MinX(source)].X;
                maxX = source[MaxX(source)].X;
                minY = source[MinY(source)].Y;
                maxY = source[MaxY(source)].Y;
            }
            catch (Exception)
            {

            }

            return new GeoLibPointF[] { new GeoLibPointF(minX, minY), new GeoLibPointF(maxX, maxY) };
        }

        public static bool enclosed(Path a, Paths b, bool strict = false)
        {
            return pEnclosed(a, b, strict);
        }

        static bool pEnclosed(Path a, Paths b, bool strict = false)
        {
            Paths aPath = new Paths();
            aPath.Add(a);
            return pEnclosed(aPath, b, strict);
        }

        // To handle this situation, a gap removal is performed. This strips the keyhole (if present) and yields enclosed cutters that we can detect and flag for rigorous processing.
        public static bool enclosed(Paths source, double customSizing, double extension, bool strict = false)
        {
            return pEnclosed(pGetOutersAndCutters(pRemoveFragments(source, customSizing, extension)), strict);
        }

        public static bool enclosed(Paths a, Paths b, bool strict = false)
        {
            return pEnclosed(a, b, strict);
        }

        public static bool enclosed(Paths[] source, bool strict = false)
        {
            return pEnclosed(source, strict);
        }
        static bool pEnclosed(Paths[] source, bool strict = false)
        {
            return pEnclosed(source[(int)outerCutterIndex.cutter], source[(int)outerCutterIndex.outer], strict);
        }

        static bool pEnclosed(Paths a, Paths b, bool strict)
        {
            bool result = false;

            if ((a.Count == 0) || (b.Count == 0))
            {
                return result;
            }

            Clipper c = new Clipper();

            Paths rationalizedFirstLayer = new Paths();
            // Force to clockwise as a safety measure.
            for (int i = 0; i < a.Count; i++)
            {
                rationalizedFirstLayer.Add(GeoWrangler.clockwise(/*Clipper.CleanPolygon*/(a[i])));
            }

            Paths rationalizedSecondLayer = new Paths();
            for (int i = 0; i < b.Count; i++)
            {
                rationalizedSecondLayer.Add(GeoWrangler.clockwise(/*Clipper.CleanPolygon*/(b[i])));
            }

            // Intersection should not matter based on order.
            Paths intersectionPaths = new Paths();
            c.AddPaths(rationalizedSecondLayer, PolyType.ptClip, true);
            c.AddPaths(rationalizedFirstLayer, PolyType.ptSubject, true);
            c.Execute(ClipType.ctUnion, intersectionPaths);

            // Force clockwise.
            for (int i = 0; i < intersectionPaths.Count; i++)
            {
                // Fix point order to ensure we can compare easily.
                Path intersectionPath = GeoWrangler.clockwise(intersectionPaths[i]);

                // Compare hashes.
                string intersectionPathHash = utility.Utils.GetMD5Hash(intersectionPath.ToString());

                bool polyMatchFound = false;

                for (int poly = 0; poly < rationalizedSecondLayer.Count; poly++)
                {
                    string myHash = utility.Utils.GetMD5Hash(rationalizedSecondLayer[poly].ToString());

                    if (myHash == intersectionPathHash)
                    {
                        polyMatchFound = true;
                    }

                    if (polyMatchFound)
                    {
                        break;
                    }
                }

                // Strict requires a to be enclosed by b and not to check the case that b is enclosed by a.
                if (!strict && !polyMatchFound)
                {
                    for (int poly = 0; poly < rationalizedFirstLayer.Count; poly++)
                    {
                        string myHash = utility.Utils.GetMD5Hash(rationalizedFirstLayer[poly].ToString());

                        if (myHash == intersectionPathHash)
                        {
                            polyMatchFound = true;
                        }

                        if (polyMatchFound)
                        {
                            break;
                        }
                    }
                }

                result = polyMatchFound;

                // Check based on area. This seems to be needed - the above doesn't do the right thing every time.
                double overlapArea = Math.Abs(Clipper.Area(intersectionPath));
                double clipArea = Math.Abs(Clipper.Area(rationalizedSecondLayer[0]));
                double subjArea = Math.Abs(Clipper.Area(rationalizedFirstLayer[0]));

                if ((overlapArea == clipArea) || (overlapArea == subjArea))
                {
                    result = true;
                }
                else
                {
                    result = false;
                }
                if (result)
                {
                    break;
                }
            }

            return result;
        }

        public static bool orthogonal(GeoLibPoint[] sourcePoly, double angularTolerance)
        {
            return pOrthogonal(sourcePoly, angularTolerance);
        }

        static bool pOrthogonal(GeoLibPoint[] sourcePoly, double angularTolerance)
        {
            bool isOrthogonal = true;

            double[] _angles = angles(sourcePoly, allowNegative: true);

            for (int i = 0; i < _angles.Length; i++)
            {
                if (Math.Abs(Math.Abs(_angles[i]) - 90.0) > angularTolerance)
                {
                    isOrthogonal = false;
                    break;
                }
            }

            return isOrthogonal;
        }

        public static double[] angles(GeoLibPoint[] sourcePoly, bool allowNegative)
        {
            return pAngles(sourcePoly, allowNegative);
        }

        static double[] pAngles(GeoLibPoint[] sourcePoly, bool allowNegative)
        {
            GeoLibPoint interSection_A, interSection_B, interSection_C;
            GeoLibPoint[] stripped = pStripTerminators(sourcePoly, true);
            int finalIndex = stripped.Length - 2;

            double[] angles = new double[finalIndex + 1];

            for (int pt = 0; pt <= finalIndex; pt++)
            {
                // Assess angle.
                if (pt == 0)
                {
                    interSection_B = stripped[finalIndex]; // map to last point
                    interSection_C = stripped[pt];
                    interSection_A = stripped[pt + 1];
                }
                else if (pt == finalIndex) // last point in the list
                {
                    interSection_B = stripped[pt - 1];
                    interSection_C = stripped[pt];
                    interSection_A = stripped[0]; // map to the first point
                }
                else
                {
                    interSection_B = stripped[pt - 1];
                    interSection_C = stripped[pt];
                    interSection_A = stripped[pt + 1];
                }

                double theta = pAngleBetweenPoints(interSection_A, interSection_B, interSection_C, allowNegative);

                angles[pt] = theta;
            }

            return angles;
        }

        public static bool orthogonal(GeoLibPointF[] sourcePoly, double angularTolerance)
        {
            return pOrthogonal(sourcePoly, angularTolerance);
        }

        static bool pOrthogonal(GeoLibPointF[] sourcePoly, double angularTolerance)
        {
            bool isOrthogonal = true;

            double[] _angles = angles(sourcePoly, allowNegative: false);

            for (int i = 0; i < _angles.Length; i++)
            {
                if (Math.Abs(Math.Abs(_angles[i]) - 90.0) > angularTolerance)
                {
                    isOrthogonal = false;
                    break;
                }
            }

            return isOrthogonal;
        }

        public static double[] angles(GeoLibPointF[] sourcePoly, bool allowNegative)
        {
            return pAngles(sourcePoly, allowNegative);
        }

        static double[] pAngles(GeoLibPointF[] sourcePoly, bool allowNegative)
        {
            GeoLibPointF interSection_A, interSection_B, interSection_C;
            GeoLibPointF[] stripped = pStripTerminators(sourcePoly, false);
            int finalIndex = stripped.Length - 1;

            double[] angles = new double[stripped.Length];

            for (int pt = 0; pt <= finalIndex; pt++)
            {
                // Assess angle.
                if (pt == 0)
                {
                    interSection_B = stripped[finalIndex]; // map to last point
                    interSection_C = stripped[pt];
                    interSection_A = stripped[pt + 1];
                }
                else if (pt == finalIndex) // last point in the list
                {
                    interSection_B = stripped[pt - 1];
                    interSection_C = stripped[pt];
                    interSection_A = stripped[0]; // map to the first point
                }
                else
                {
                    interSection_B = stripped[pt - 1];
                    interSection_C = stripped[pt];
                    interSection_A = stripped[pt + 1];
                }

                double theta = pAngleBetweenPoints(interSection_A, interSection_B, interSection_C, allowNegative);

                angles[pt] = theta;
            }

            return angles;
        }

        public static bool orthogonal(Path sourcePoly, double angularTolerance)
        {
            return pOrthogonal(sourcePoly, angularTolerance);
        }

        static bool pOrthogonal(Path sourcePoly, double angularTolerance)
        {
            bool isOrthogonal = true;

            double[] _angles = angles(sourcePoly, allowNegative: false);

            for (int i = 0; i < _angles.Length; i++)
            {
                if (Math.Abs(Math.Abs(_angles[i]) - 90.0) > angularTolerance)
                {
                    isOrthogonal = false;
                    break;
                }
            }

            return isOrthogonal;
        }

        public static double[] angles(Path sourcePoly, bool allowNegative)
        {
            return pAngles(sourcePoly, allowNegative);
        }

        static double[] pAngles(Path sourcePoly, bool allowNegative)
        {
            IntPoint interSection_A, interSection_B, interSection_C;

            Path stripped = pStripTerminators(sourcePoly, false);
            int finalIndex = stripped.Count - 1;

            double[] angles = new double[stripped.Count];

            for (int pt = 0; pt <= finalIndex; pt++)
            {
                // Assess angle.
                if (pt == 0)
                {
                    interSection_B = stripped[finalIndex]; // map to last point
                    interSection_C = stripped[pt];
                    interSection_A = stripped[pt + 1];
                }
                else if (pt == finalIndex) // last point in the list
                {
                    interSection_B = stripped[pt - 1];
                    interSection_C = stripped[pt];
                    interSection_A = stripped[0]; // map to the first point
                }
                else
                {
                    interSection_B = stripped[pt - 1];
                    interSection_C = stripped[pt];
                    interSection_A = stripped[pt + 1];
                }

                double theta = pAngleBetweenPoints(interSection_A, interSection_B, interSection_C, allowNegative);

                angles[pt] = theta;
            }

            return angles;
        }

        public static bool isClockwise(GeoLibPoint[] points)
        {
            return pIsClockwise(points);
        }

        static bool pIsClockwise(GeoLibPoint[] points)
        {
            // Based on stackoverflow.com/questions/1165647/how-to-determine-if-a-list-of-polygon-points-are-in-clockwise-order
            // Shoelace formula.

            Int64 delta = 0;

            for (Int32 pt = 0; pt < points.Length; pt++)
            {
                Int64 deltaX = 0;
                Int64 deltaY = 0;
                if (pt == points.Length - 1)
                {
                    deltaX = (points[0].X - points[pt].X);
                    deltaY = (points[0].Y + points[pt].Y);
                }
                else
                {
                    deltaX = (points[pt + 1].X - points[pt].X);
                    deltaY = (points[pt + 1].Y + points[pt].Y);
                }

                delta += deltaX * deltaY;
            }

            if (delta > 0)
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        public static bool isClockwise(GeoLibPointF[] points)
        {
            return pIsClockwise(points);
        }

        static bool pIsClockwise(GeoLibPointF[] points)
        {
            // Based on stackoverflow.com/questions/1165647/how-to-determine-if-a-list-of-polygon-points-are-in-clockwise-order
            // Shoelace formula.

            double delta = 0;

            for (Int32 pt = 0; pt < points.Length; pt++)
            {
                double deltaX = 0;
                double deltaY = 0;
                if (pt == points.Length - 1)
                {
                    deltaX = (points[0].X - points[pt].X);
                    deltaY = (points[0].Y + points[pt].Y);
                }
                else
                {
                    deltaX = (points[pt + 1].X - points[pt].X);
                    deltaY = (points[pt + 1].Y + points[pt].Y);
                }

                delta += deltaX * deltaY;
            }

            if (delta > 0)
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        public static bool isClockwise(Path points)
        {
            return pIsClockwise(points);
        }

        static bool pIsClockwise(Path points)
        {
            // Based on stackoverflow.com/questions/1165647/how-to-determine-if-a-list-of-polygon-points-are-in-clockwise-order
            // Shoelace formula.

            Int64 delta = 0;

            for (Int32 pt = 0; pt < points.Count; pt++)
            {
                Int64 deltaX = 0;
                Int64 deltaY = 0;
                if (pt == points.Count - 1)
                {
                    deltaX = (points[0].X - points[pt].X);
                    deltaY = (points[0].Y + points[pt].Y);
                }
                else
                {
                    deltaX = (points[pt + 1].X - points[pt].X);
                    deltaY = (points[pt + 1].Y + points[pt].Y);
                }

                delta += deltaX * deltaY;
            }

            if (delta > 0)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
    }
}
