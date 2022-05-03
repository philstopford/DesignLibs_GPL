using Clipper2Lib;
using geoLib;
using System.Collections.Generic;

namespace geoWrangler;

public static partial class GeoWrangler
{
    // Bounds will force the negation to only work against the extents of the shape. Useful for capturing islands in negative tone.
    public static List<GeoLibPointF[]> invertTone(GeoLibPointF[] source, long scaleFactor, bool preserveColinear, bool useTriangulation = false, bool useBounds = false)
    {
        return pPointFsFromPaths(pInvertTone(pPathsFromPointFs(new List<GeoLibPointF[]> { source }, scaleFactor), preserveColinear: preserveColinear, useTriangulation: useTriangulation, useBounds: useBounds), scaleFactor);
    }

    public static List<GeoLibPointF[]> invertTone(List<GeoLibPointF[]> source, long scaleFactor, bool preserveColinear, bool useTriangulation = false, bool useBounds = false)
    {
        return pPointFsFromPaths(pInvertTone(pPathsFromPointFs(source, scaleFactor), preserveColinear: preserveColinear, useTriangulation: useTriangulation, useBounds: useBounds), scaleFactor);
    }

    public static Paths64 invertTone(Path64 sourcePath, bool preserveColinear, bool useTriangulation = false, bool useBounds = false)
    {
        Paths64 t = new() {sourcePath};
        return invertTone(t, preserveColinear: preserveColinear, useTriangulation: useTriangulation, useBounds:useBounds);
    }

    public static Paths64 invertTone(Paths64 sourcePaths, bool preserveColinear, bool useTriangulation = false, bool useBounds = false)
    {
        return pInvertTone(sourcePaths, preserveColinear: preserveColinear, useTriangulation: useTriangulation, useBounds: useBounds);
    }

    private static Paths64 pInvertTone(Paths64 sourcePaths, bool preserveColinear, bool useTriangulation, bool useBounds)
    {
        Path64 firstLayerBP = new();
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
                Rect64 bounds = ClipperFunc.GetBounds(sourcePaths);
                firstLayerBP.Add(new Point64(bounds.left, bounds.bottom));
                firstLayerBP.Add(new Point64(bounds.left, bounds.top));
                firstLayerBP.Add(new Point64(bounds.right, bounds.top));
                firstLayerBP.Add(new Point64(bounds.right, bounds.bottom));
                firstLayerBP.Add(new Point64(bounds.left, bounds.bottom));
                break;
            }
        }

        Clipper c = new() {PreserveCollinear = preserveColinear};

        c.AddSubject(firstLayerBP);
        // Add hole polygons from our paths
        c.AddClip(sourcePaths);

        Paths64 cutters = new();
        c.Execute(ClipType.Difference, FillRule.EvenOdd, cutters);

        cutters = pReorder(cutters);
        
        switch (useTriangulation)
        {
            case false:
                sourcePaths = cutters;
                return sourcePaths;
            default:
                return pMakeKeyHole(sourcePaths);
        }
    }
}