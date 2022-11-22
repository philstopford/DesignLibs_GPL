using Clipper2Lib;
using geoLib;
using System.Collections.Generic;

namespace geoWrangler;

using Path = Path64;
using Paths = Paths64;

public static partial class GeoWrangler
{
    // Bounds will force the negation to only work against the extents of the shape. Useful for capturing islands in negative tone.
    public static PathsD invertTone(PathD source, bool preserveColinear, bool useTriangulation = false, bool useBounds = false)
    {
        return pPathsDFromPaths64(pInvertTone(pPaths64FromPathsD(new PathsD { source }), preserveColinear: preserveColinear, useTriangulation: useTriangulation, useBounds: useBounds));
    }

    public static PathsD invertTone(PathsD source, bool preserveColinear, bool useTriangulation = false, bool useBounds = false)
    {
        return pPathsDFromPaths64(pInvertTone(pPaths64FromPathsD(source), preserveColinear: preserveColinear, useTriangulation: useTriangulation, useBounds: useBounds));
    }

    public static Paths invertTone(Path sourcePath, bool preserveColinear, bool useTriangulation = false, bool useBounds = false)
    {
        Paths t = new() {sourcePath};
        return invertTone(t, preserveColinear: preserveColinear, useTriangulation: useTriangulation, useBounds:useBounds);
    }

    public static Paths invertTone(Paths sourcePaths, bool preserveColinear, bool useTriangulation = false, bool useBounds = false)
    {
        return pInvertTone(sourcePaths, preserveColinear: preserveColinear, useTriangulation: useTriangulation, useBounds: useBounds);
    }

    private static Paths pInvertTone(Paths sourcePaths, bool preserveColinear, bool useTriangulation, bool useBounds)
    {
        Path firstLayerBP = new();
        switch (useBounds)
        {
            case false:
                firstLayerBP.Add(new Point64(-int.MaxValue, -int.MaxValue));
                firstLayerBP.Add(new Point64(-int.MaxValue, int.MaxValue));
                firstLayerBP.Add(new Point64(int.MaxValue, int.MaxValue));
                firstLayerBP.Add(new Point64(int.MaxValue, -int.MaxValue));
                firstLayerBP.Add(new Point64(-int.MaxValue, -int.MaxValue));
                break;
            default:
            {
                Rect64 bounds = Clipper.GetBounds(sourcePaths);
                firstLayerBP.Add(new Point64(bounds.left, bounds.bottom));
                firstLayerBP.Add(new Point64(bounds.left, bounds.top));
                firstLayerBP.Add(new Point64(bounds.right, bounds.top));
                firstLayerBP.Add(new Point64(bounds.right, bounds.bottom));
                firstLayerBP.Add(new Point64(bounds.left, bounds.bottom));
                break;
            }
        }

        Clipper64 c = new() {PreserveCollinear = preserveColinear};

        c.AddSubject(firstLayerBP);
        // Add hole polygons from our paths
        c.AddClip(sourcePaths);

        Paths cutters = new();
        c.Execute(ClipType.Difference, FillRule.EvenOdd, cutters);

        cutters = pReorderXY(cutters);
        
        switch (useTriangulation)
        {
            case false:
                sourcePaths = cutters;
                return sourcePaths;
            default:
                return pMakeKeyHole(sourcePaths, reverseEval:false, biDirectionalEval:true);
        }
    }
}