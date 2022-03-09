using ClipperLib2;
using geoLib;
using System.Collections.Generic;
using System.Linq;

namespace geoWrangler;

using Path = List<Point64>;
using Paths = List<List<Point64>>;

public static partial class GeoWrangler
{
    // Bounds will force the negation to only work against the extents of the shape. Useful for capturing islands in negative tone.
    public static List<GeoLibPointF[]> invertTone(GeoLibPointF[] source, long scaleFactor, bool useTriangulation = false, bool useBounds = false)
    {
        return pPointFsFromPaths(pInvertTone(pPathsFromPointFs(new List<GeoLibPointF[]> { source }, scaleFactor), useTriangulation, useBounds), scaleFactor);
    }

    public static List<GeoLibPointF[]> invertTone(List<GeoLibPointF[]> source, long scaleFactor, bool useTriangulation = false, bool useBounds = false)
    {
        return pPointFsFromPaths(pInvertTone(pPathsFromPointFs(source, scaleFactor), useTriangulation, useBounds), scaleFactor);
    }

    public static Paths invertTone(Path sourcePath, bool useTriangulation = false, bool useBounds = false)
    {
        Paths t = new() {sourcePath};
        return invertTone(t, useTriangulation, useBounds);
    }

    public static Paths invertTone(Paths sourcePaths, bool useTriangulation = false, bool useBounds = false)
    {
        return pInvertTone(sourcePaths, useTriangulation, useBounds);
    }

    private static Paths pInvertTone(Paths sourcePaths, bool useTriangulation, bool useBounds)
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
                Rect64 bounds = ClipperFunc.GetBounds(sourcePaths);
                firstLayerBP.Add(new Point64(bounds.left, bounds.bottom));
                firstLayerBP.Add(new Point64(bounds.left, bounds.top));
                firstLayerBP.Add(new Point64(bounds.right, bounds.top));
                firstLayerBP.Add(new Point64(bounds.right, bounds.bottom));
                firstLayerBP.Add(new Point64(bounds.left, bounds.bottom));
                break;
            }
        }

        Clipper c = new();

        c.AddSubject(firstLayerBP);
        // Add hole polygons from our paths
        c.AddClip(sourcePaths);

        Paths cutters = new();
        c.Execute(ClipType.Difference, FillRule.EvenOdd, cutters);

        cutters = pStripColinear(cutters);

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