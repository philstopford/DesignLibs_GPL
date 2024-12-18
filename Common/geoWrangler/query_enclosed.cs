using System;
using System.Linq;
using Clipper2Lib;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static bool anyPartialOverlap(PathsD a, PathsD b)
    {
        double a_area = Clipper.Area(a);
        double b_area = Clipper.Area(b);
        ClipperD c = new(Constants.roundingDecimalPrecision);

        PathsD rationalizedFirstLayer = clockwise(a);
        // Force to clockwise as a safety measure.

        PathsD rationalizedSecondLayer = clockwise(b);

        // Intersection should not matter based on order.
        PathsD intersectionPaths = [];
        c.AddClip(rationalizedFirstLayer);
        c.AddSubject(rationalizedSecondLayer);
        c.Execute(ClipType.Intersection, FillRule.EvenOdd, intersectionPaths);

        double intersectionArea = Math.Abs(Clipper.Area(intersectionPaths));
        // We have an intersection, let's check for any full enclosure scenario.
        if (intersectionArea > 0)
        {
            return (Math.Abs(Math.Abs(intersectionArea) - Math.Abs(a_area)) > Constants.tolerance) && 
                   (Math.Abs(Math.Abs(intersectionArea) - Math.Abs(b_area)) > Constants.tolerance);
            // Totally enclosed.
            // Partial overlap
        }

        // Default to no partial overlap.
        return false;
    }
    public static bool enclosed(PathD a, PathsD b, bool strict = false)
    {
        return pEnclosed(a, b, strict);
    }

    private static bool pEnclosed(PathD a, PathsD b, bool strict = false)
    {
        PathsD aPath = [a];
        return pEnclosed(aPath, b, strict);
    }

    // To handle this situation, a gap removal is performed. This strips the keyhole (if present) and yields enclosed cutters that we can detect and flag for rigorous processing.
    public static bool enclosed(PathsD source, double customSizing, double extension, bool strict = false)
    {
        return pEnclosed(pGetOutersAndCutters(pRemoveFragments(source, customSizing, extension)), strict);
    }

    public static bool enclosed(PathsD a, PathsD b, bool strict = false)
    {
        return pEnclosed(a, b, strict);
    }

    public static bool enclosed(PathsD[] source, bool strict = false)
    {
        return pEnclosed(source, strict);
    }

    private static bool pEnclosed(PathsD[] source, bool strict = false)
    {
        return pEnclosed(source[(int)outerCutterIndex.cutter], source[(int)outerCutterIndex.outer], strict);
    }

    private static bool pEnclosed(PathsD a, PathsD b, bool strict)
    {
        if (a.Count == 0 || b.Count == 0)
        {
            return false;
        }

        // Force to clockwise as a safety measure.
        PathsD rationalizedFirstLayer = clockwise(a);
        PathsD rationalizedSecondLayer = clockwise(b);

        if (anyPartialOverlap(rationalizedFirstLayer, rationalizedSecondLayer))
        {
            return false;
        }

        bool result = false;
        
        ClipperD c = new(Constants.roundingDecimalPrecision);

        // Intersection should not matter based on order.
        PathsD intersectionPaths = [];
        c.AddClip(rationalizedSecondLayer);
        c.AddSubject(rationalizedFirstLayer);
        c.Execute(ClipType.Intersection, FillRule.EvenOdd, intersectionPaths);
        if (!intersectionPaths.Any())
        {
            return false;
        }

        intersectionPaths = pReorderXY(intersectionPaths);

        // Force clockwise.
        foreach (PathD t in intersectionPaths)
        {
            // Fix point order to ensure we can compare easily.
            PathD intersectionPath = clockwise(t);

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

            if (Math.Abs(overlapArea - clipArea) <= Constants.tolerance || Math.Abs(overlapArea - subjArea) <= Constants.tolerance)
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
}