using Clipper2Lib;
using geoLib;
using System;
using System.Collections.Generic;
using System.Linq;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static PointD midPoint(Paths64 source)
    {
        return pMidPoint(source);
    }

    private static PointD pMidPoint(Paths64 source)
    {
        Path64 bounds = getBounds(source[0]);

        double minX = bounds[0].X;
        double minY = bounds[0].Y;
        double maxX = bounds[1].X;
        double maxY = bounds[1].Y;

        for (int i = 1; i < source.Count; i++)
        {
            Path64 tmp = getBounds(source[i]);
            minX = Math.Min(minX, tmp[0].X);
            minY = Math.Min(minY, tmp[0].Y);
            maxX = Math.Max(maxX, tmp[1].X);
            maxY = Math.Max(maxY, tmp[1].Y);
        }

        double avX = minX + (maxX - minX) / 2.0f;
        double avY = minY + (maxY - minY) / 2.0f;

        return new (avX, avY);
    }

    public static PointD midPoint(Path64 source)
    {
        return pMidPoint(source);
    }

    private static PointD pMidPoint(Path64 source)
    {
        Path64 bounds = getBounds(source);

        double minX = bounds[0].X;
        double minY = bounds[0].Y;
        double maxX = bounds[1].X;
        double maxY = bounds[1].Y;

        double avX = minX + (maxX - minX) / 2.0f;
        double avY = minY + (maxY - minY) / 2.0f;

        return new (avX, avY);
    }

    public static PointD midPoint(PathsD source)
    {
        return pMidPoint(source);
    }

    private static PointD pMidPoint(PathsD source)
    {
        PathD bounds = getBounds(source[0]);

        double minX = bounds[0].x;
        double minY = bounds[0].y;
        double maxX = bounds[1].x;
        double maxY = bounds[1].y;

        for (int i = 1; i < source.Count; i++)
        {
            PathD tmp = getBounds(source[i]);
            minX = Math.Min(minX, tmp[0].x);
            minY = Math.Min(minY, tmp[0].y);
            maxX = Math.Max(maxX, tmp[1].x);
            maxY = Math.Max(maxY, tmp[1].y);
        }

        double avX = minX + (maxX - minX) / 2.0f;
        double avY = minY + (maxY - minY) / 2.0f;

        return new (avX, avY);
    }

    public static PointD midPoint(PathD source)
    {
        return pMidPoint(source);
    }

    private static PointD pMidPoint(PathD source)
    {
        PathD bounds = getBounds(source);

        double minX = bounds[0].x;
        double minY = bounds[0].y;
        double maxX = bounds[1].x;
        double maxY = bounds[1].y;

        double avX = minX + (maxX - minX) / 2.0f;
        double avY = minY + (maxY - minY) / 2.0f;

        return new (avX, avY);
    }

    public static PointD getExtents(PathD source)
    {
        PathD bounds = pGetBounds(source);
        PointD extents = pDistanceBetweenPoints_point(bounds[0], bounds[2]);

        return extents;
    }

    public static Path64 getBounds(Path64 source)
    {
        return pGetBounds(source);
    }

    private static Path64 pGetBounds(Path64 source)
    {
        double minX = 0;
        double minY = 0;
        double maxX = 0;
        double maxY = 0;

        if (source.Count > 0)
        {
            try
            {
                minX = source[MinX(source)].X;
                maxX = source[MaxX(source)].X;
                minY = source[MinY(source)].Y;
                maxY = source[MaxY(source)].Y;
            }
            catch
            {
                // ignored
            }
        }

        return new () { new (minX, minY), new (maxX, maxY) };
    }

    public static PathD getBounds(PathD source)
    {
        return pGetBounds(source);
    }

    private static PathD pGetBounds(PathD source)
    {
        double minX = 0;
        double minY = 0;
        double maxX = 0;
        double maxY = 0;

        if (source.Count > 0)
        {
            try
            {
                minX = source[MinX(source)].x;
                maxX = source[MaxX(source)].x;
                minY = source[MinY(source)].y;
                maxY = source[MaxY(source)].y;
            }
            catch
            {
                // ignored
            }
        }

        return new () { new (minX, minY), new (maxX, maxY) };
    }

    public static bool enclosed(Path64 a, Paths64 b, bool strict = false)
    {
        return pEnclosed(a, b, strict);
    }

    private static bool pEnclosed(Path64 a, Paths64 b, bool strict = false)
    {
        Paths64 aPath = new() {a};
        return pEnclosed(aPath, b, strict);
    }

    // To handle this situation, a gap removal is performed. This strips the keyhole (if present) and yields enclosed cutters that we can detect and flag for rigorous processing.
    public static bool enclosed(Paths64 source, double customSizing, double extension, bool strict = false)
    {
        return pEnclosed(pGetOutersAndCutters(pRemoveFragments(source, customSizing, extension)), strict);
    }

    public static bool enclosed(Paths64 a, Paths64 b, bool strict = false)
    {
        return pEnclosed(a, b, strict);
    }

    public static bool enclosed(Paths64[] source, bool strict = false)
    {
        return pEnclosed(source, strict);
    }

    private static bool pEnclosed(Paths64[] source, bool strict = false)
    {
        return pEnclosed(source[(int)outerCutterIndex.cutter], source[(int)outerCutterIndex.outer], strict);
    }

    private static bool pEnclosed(Paths64 a, Paths64 b, bool strict)
    {

        if (a.Count == 0 || b.Count == 0)
        {
            return false;
        }

        bool result = false;
            
        Clipper64 c = new();

        Paths64 rationalizedFirstLayer = clockwise(a);
        // Force to clockwise as a safety measure.

        Paths64 rationalizedSecondLayer = clockwise(b);

        // Intersection should not matter based on order.
        Paths64 intersectionPaths = new();
        c.AddClip(rationalizedSecondLayer);
        c.AddSubject(rationalizedFirstLayer);
        c.Execute(ClipType.Union, FillRule.EvenOdd, intersectionPaths);

        intersectionPaths = pReorderXY(intersectionPaths);

        // Force clockwise.
        foreach (Path64 t in intersectionPaths)
        {
            // Fix point order to ensure we can compare easily.
            Path64 intersectionPath = clockwise(t);

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
            double overlapArea = Math.Abs(Clipper.Area(intersectionPath));
            double clipArea = Math.Abs(Clipper.Area(rationalizedSecondLayer[0]));
            double subjArea = Math.Abs(Clipper.Area(rationalizedFirstLayer[0]));

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

    public static bool orthogonal(Path64 sourcePoly, double angularTolerance)
    {
        return pOrthogonal(sourcePoly, angularTolerance);
    }

    private static bool pOrthogonal(Path64 sourcePoly, double angularTolerance)
    {
        double[] _angles = angles(sourcePoly, allowNegative: true);

        return _angles.All(t => !(Math.Abs(Math.Abs(t) - 90.0) > angularTolerance));
    }

    public static double[] angles(Path64 sourcePoly, bool allowNegative)
    {
        return pAngles(sourcePoly, allowNegative);
    }

    private static double[] pAngles(Path64 sourcePoly, bool allowNegative)
    {
        Path64 stripped = pStripTerminators(sourcePoly, true);
        int finalIndex = stripped.Count - 2;

        double[] angles = new double[finalIndex + 1];

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

    public static bool orthogonal(PathD sourcePoly, double angularTolerance)
    {
        return pOrthogonal(sourcePoly, angularTolerance);
    }

    private static bool pOrthogonal(PathD sourcePoly, double angularTolerance)
    {
        double[] _angles = angles(sourcePoly, allowNegative: false);

        return _angles.All(t => !(Math.Abs(Math.Abs(t) - 90.0) > angularTolerance));
    }

    public static double[] angles(PathD sourcePoly, bool allowNegative)
    {
        return pAngles(sourcePoly, allowNegative);
    }

    private static double[] pAngles(PathD sourcePoly, bool allowNegative)
    {
        PathD stripped = pStripTerminators(sourcePoly, false);
        int finalIndex = stripped.Count - 1;

        double[] angles = new double[stripped.Count];

        for (int pt = 0; pt <= finalIndex; pt++)
        {
            // Assess angle.
            PointD interSection_A;
            PointD interSection_B;
            PointD interSection_C;
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

    public static bool orthogonal(Paths64 sourcePoly, double angularTolerance)
    {
        return pOrthogonal(sourcePoly, angularTolerance);
    }

    public static bool pOrthogonal(Paths64 sourcePoly, double angularTolerance)
    {
        bool ret = true;
        foreach (Path64 p in sourcePoly)
        {
            ret = ret && pOrthogonal(p, angularTolerance);
        }

        return ret;
    }
    
    public static bool isClockwise(Path64 points)
    {
        return pIsClockwise(points);
    }

    private static bool pIsClockwise(Path64 points)
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

    public static bool isClockwise(PathD points)
    {
        return pIsClockwise(points);
    }

    private static bool pIsClockwise(PathD points)
    {
        // Based on stackoverflow.com/questions/1165647/how-to-determine-if-a-list-of-polygon-points-are-in-clockwise-order
        // Shoelace formula.

        double delta = 0;

        for (int pt = 0; pt < points.Count; pt++)
        {
            double deltaX;
            double deltaY;
            if (pt == points.Count - 1)
            {
                deltaX = points[0].x - points[pt].x;
                deltaY = points[0].y + points[pt].y;
            }
            else
            {
                deltaX = points[pt + 1].x - points[pt].x;
                deltaY = points[pt + 1].y + points[pt].y;
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