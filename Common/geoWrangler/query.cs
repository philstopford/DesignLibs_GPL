using ClipperLib2;
using geoLib;
using System;
using System.Collections.Generic;
using System.Linq;

namespace geoWrangler;

using Path = List<Point64>;
using Paths = List<List<Point64>>;

public static partial class GeoWrangler
{
    public static GeoLibPointF midPoint(List<GeoLibPoint[]> source)
    {
        return pMidPoint(source);
    }

    private static GeoLibPointF pMidPoint(List<GeoLibPoint[]> source)
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

    private static GeoLibPointF pMidPoint(GeoLibPoint[] source)
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

    private static GeoLibPointF pMidPoint(List<GeoLibPointF[]> source)
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

    private static GeoLibPointF pMidPoint(GeoLibPointF[] source)
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
        GeoLibPointF extents = pDistanceBetweenPoints_point(bounds[0], bounds[2]);

        return extents;
    }

    public static GeoLibPoint[] getBounds(GeoLibPoint[] source)
    {
        return pGetBounds(source);
    }

    private static GeoLibPoint[] pGetBounds(GeoLibPoint[] source)
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

        return new [] { new GeoLibPoint(minX, minY), new GeoLibPoint(maxX, maxY) };
    }

    public static GeoLibPointF[] getBounds(GeoLibPointF[] source)
    {
        return pGetBounds(source);
    }

    private static GeoLibPointF[] pGetBounds(GeoLibPointF[] source)
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

        return new [] { new GeoLibPointF(minX, minY), new GeoLibPointF(maxX, maxY) };
    }

    public static bool enclosed(Path a, Paths b, bool strict = false)
    {
        return pEnclosed(a, b, strict);
    }

    private static bool pEnclosed(Path a, Paths b, bool strict = false)
    {
        Paths aPath = new() {a};
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

    private static bool pEnclosed(Paths[] source, bool strict = false)
    {
        return pEnclosed(source[(int)outerCutterIndex.cutter], source[(int)outerCutterIndex.outer], strict);
    }

    private static bool pEnclosed(Paths a, Paths b, bool strict)
    {

        if (a.Count == 0 || b.Count == 0)
        {
            return false;
        }

        bool result = false;
            
        Clipper c = new();

        Paths rationalizedFirstLayer = a.Select(t => clockwise( /*Clipper.CleanPolygon*/t)).ToList();
        // Force to clockwise as a safety measure.

        Paths rationalizedSecondLayer = b.Select(t => clockwise( /*Clipper.CleanPolygon*/t)).ToList();

        // Intersection should not matter based on order.
        PolyTree pt = new();
        c.AddClip(rationalizedSecondLayer);
        c.AddSubject(rationalizedFirstLayer);
        c.Execute(ClipType.Union, FillRule.EvenOdd, pt);
        Paths intersectionPaths = ClipperFunc.PolyTreeToPaths(pt);

        // Force clockwise.
        foreach (Path t in intersectionPaths)
        {
            // Fix point order to ensure we can compare easily.
            Path intersectionPath = clockwise(t);

            // Compare hashes.
            string intersectionPathHash = utility.Utils.GetMD5Hash(intersectionPath.ToString());

            bool polyMatchFound = false;

            foreach (string myHash in rationalizedSecondLayer.Select(t1 => utility.Utils.GetMD5Hash(t1.ToString())))
            {
                if (myHash == intersectionPathHash)
                {
                    polyMatchFound = true;
                }

                if (polyMatchFound)
                {
                    break;
                }
            }

            switch (strict)
            {
                // Strict requires a to be enclosed by b and not to check the case that b is enclosed by a.
                case false when !polyMatchFound:
                {
                    foreach (string myHash in rationalizedFirstLayer.Select(t1 => utility.Utils.GetMD5Hash(t1.ToString())))
                    {
                        if (myHash == intersectionPathHash)
                        {
                            polyMatchFound = true;
                        }

                        if (polyMatchFound)
                        {
                            break;
                        }
                    }

                    break;
                }
            }
                
            // Check based on area. This seems to be needed - the above doesn't do the right thing every time.
            double overlapArea = Math.Abs(ClipperFunc.Area(intersectionPath));
            double clipArea = Math.Abs(ClipperFunc.Area(rationalizedSecondLayer[0]));
            double subjArea = Math.Abs(ClipperFunc.Area(rationalizedFirstLayer[0]));

            if (Math.Abs(overlapArea - clipArea) <= double.Epsilon || Math.Abs(overlapArea - subjArea) <= double.Epsilon)
            {
                result = true;
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

    private static bool pOrthogonal(GeoLibPoint[] sourcePoly, double angularTolerance)
    {
        double[] _angles = angles(sourcePoly, allowNegative: true);

        return _angles.All(t => !(Math.Abs(Math.Abs(t) - 90.0) > angularTolerance));
    }

    public static double[] angles(GeoLibPoint[] sourcePoly, bool allowNegative)
    {
        return pAngles(sourcePoly, allowNegative);
    }

    private static double[] pAngles(GeoLibPoint[] sourcePoly, bool allowNegative)
    {
        GeoLibPoint[] stripped = pStripTerminators(sourcePoly, true);
        int finalIndex = stripped.Length - 2;

        double[] angles = new double[finalIndex + 1];

        for (int pt = 0; pt <= finalIndex; pt++)
        {
            // Assess angle.
            GeoLibPoint interSection_A;
            GeoLibPoint interSection_B;
            GeoLibPoint interSection_C;
            switch (pt)
            {
                case 0:
                    interSection_B = stripped[finalIndex]; // map to last point
                    interSection_C = stripped[pt];
                    interSection_A = stripped[pt + 1];
                    break;
                default:
                {
                    if (pt == finalIndex) // last point in the list
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

                    break;
                }
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

    private static bool pOrthogonal(GeoLibPointF[] sourcePoly, double angularTolerance)
    {
        double[] _angles = angles(sourcePoly, allowNegative: false);

        return _angles.All(t => !(Math.Abs(Math.Abs(t) - 90.0) > angularTolerance));
    }

    public static double[] angles(GeoLibPointF[] sourcePoly, bool allowNegative)
    {
        return pAngles(sourcePoly, allowNegative);
    }

    private static double[] pAngles(GeoLibPointF[] sourcePoly, bool allowNegative)
    {
        GeoLibPointF[] stripped = pStripTerminators(sourcePoly, false);
        int finalIndex = stripped.Length - 1;

        double[] angles = new double[stripped.Length];

        for (int pt = 0; pt <= finalIndex; pt++)
        {
            // Assess angle.
            GeoLibPointF interSection_A;
            GeoLibPointF interSection_B;
            GeoLibPointF interSection_C;
            switch (pt)
            {
                case 0:
                    interSection_B = stripped[finalIndex]; // map to last point
                    interSection_C = stripped[pt];
                    interSection_A = stripped[pt + 1];
                    break;
                default:
                {
                    if (pt == finalIndex) // last point in the list
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

                    break;
                }
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

    private static bool pOrthogonal(Path sourcePoly, double angularTolerance)
    {
        bool isOrthogonal = true;

        double[] _angles = angles(sourcePoly, allowNegative: false);

        foreach (double t in _angles)
        {
            if (Math.Abs(Math.Abs(t) - 90.0) > angularTolerance)
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

    private static double[] pAngles(Path sourcePoly, bool allowNegative)
    {
        Path stripped = pStripTerminators(sourcePoly, false);
        int finalIndex = stripped.Count - 1;

        double[] angles = new double[stripped.Count];

        for (int pt = 0; pt <= finalIndex; pt++)
        {
            // Assess angle.
            Point64 interSection_A;
            Point64 interSection_B;
            Point64 interSection_C;
            switch (pt)
            {
                case 0:
                    interSection_B = stripped[finalIndex]; // map to last point
                    interSection_C = stripped[pt];
                    interSection_A = stripped[pt + 1];
                    break;
                default:
                {
                    if (pt == finalIndex) // last point in the list
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

                    break;
                }
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

    private static bool pIsClockwise(GeoLibPoint[] points)
    {
        // Based on stackoverflow.com/questions/1165647/how-to-determine-if-a-list-of-polygon-points-are-in-clockwise-order
        // Shoelace formula.

        long delta = 0;

        for (int pt = 0; pt < points.Length; pt++)
        {
            long deltaX;
            long deltaY;
            if (pt == points.Length - 1)
            {
                deltaX = points[0].X - points[pt].X;
                deltaY = points[0].Y + points[pt].Y;
            }
            else
            {
                deltaX = points[pt + 1].X - points[pt].X;
                deltaY = points[pt + 1].Y + points[pt].Y;
            }

            delta += deltaX * deltaY;
        }

        return delta switch
        {
            > 0 => true,
            _ => false
        };
    }

    public static bool isClockwise(GeoLibPointF[] points)
    {
        return pIsClockwise(points);
    }

    private static bool pIsClockwise(GeoLibPointF[] points)
    {
        // Based on stackoverflow.com/questions/1165647/how-to-determine-if-a-list-of-polygon-points-are-in-clockwise-order
        // Shoelace formula.

        double delta = 0;

        for (int pt = 0; pt < points.Length; pt++)
        {
            double deltaX;
            double deltaY;
            if (pt == points.Length - 1)
            {
                deltaX = points[0].X - points[pt].X;
                deltaY = points[0].Y + points[pt].Y;
            }
            else
            {
                deltaX = points[pt + 1].X - points[pt].X;
                deltaY = points[pt + 1].Y + points[pt].Y;
            }

            delta += deltaX * deltaY;
        }

        return delta switch
        {
            > 0 => true,
            _ => false
        };
    }

    public static bool isClockwise(Path points)
    {
        return pIsClockwise(points);
    }

    private static bool pIsClockwise(Path points)
    {
        // Based on stackoverflow.com/questions/1165647/how-to-determine-if-a-list-of-polygon-points-are-in-clockwise-order
        // Shoelace formula.

        long delta = 0;

        for (int pt = 0; pt < points.Count; pt++)
        {
            long deltaX;
            long deltaY;
            if (pt == points.Count - 1)
            {
                deltaX = points[0].X - points[pt].X;
                deltaY = points[0].Y + points[pt].Y;
            }
            else
            {
                deltaX = points[pt + 1].X - points[pt].X;
                deltaY = points[pt + 1].Y + points[pt].Y;
            }

            delta += deltaX * deltaY;
        }

        return delta switch
        {
            > 0 => true,
            _ => false
        };
    }
}