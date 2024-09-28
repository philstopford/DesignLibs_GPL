using Clipper2Lib;

namespace geoWrangler;

public static partial class GeoWrangler
{
    // Bounds will force the negation to only work against the extents of the shape. Useful for capturing islands in negative tone.
    public static PathsD invertTone(PathD source, bool preserveCollinear, bool useTriangulation = false, bool useBounds = false)
    {
        return pInvertTone([source], preserveCollinear: preserveCollinear, useTriangulation: useTriangulation, useBounds: useBounds);
    }

    public static PathsD invertTone(PathsD source, bool preserveCollinear, bool useTriangulation = false, bool useBounds = false)
    {
        return pInvertTone(source, preserveCollinear: preserveCollinear, useTriangulation: useTriangulation, useBounds: useBounds);
    }
    
    private static PathsD pInvertTone(PathsD sourcePaths, bool preserveCollinear, bool useTriangulation, bool useBounds)
    {
        PathD firstLayerBP = [];
        switch (useBounds)
        {
            case false:
                firstLayerBP.Add(new PointD(-int.MaxValue, -int.MaxValue));
                firstLayerBP.Add(new PointD(-int.MaxValue, int.MaxValue));
                firstLayerBP.Add(new PointD(int.MaxValue, int.MaxValue));
                firstLayerBP.Add(new PointD(int.MaxValue, -int.MaxValue));
                firstLayerBP.Add(new PointD(-int.MaxValue, -int.MaxValue));
                break;
            default:
            {
                RectD bounds = Clipper.GetBounds(sourcePaths);
                firstLayerBP.Add(new PointD(bounds.left, bounds.bottom));
                firstLayerBP.Add(new PointD(bounds.left, bounds.top));
                firstLayerBP.Add(new PointD(bounds.right, bounds.top));
                firstLayerBP.Add(new PointD(bounds.right, bounds.bottom));
                firstLayerBP.Add(new PointD(bounds.left, bounds.bottom));
                break;
            }
        }

        ClipperD c = new(Constants.roundingDecimalPrecision) {PreserveCollinear = preserveCollinear};

        c.AddSubject(firstLayerBP);
        // Add hole polygons from our paths
        c.AddClip(sourcePaths);

        PathsD cutters = [];
        c.Execute(ClipType.Difference, FillRule.EvenOdd, cutters);

        cutters = pReorderXY(cutters);

        return useTriangulation switch
        {
            false => cutters,
            _ => pMakeKeyHole(cutters, reverseEval: false, biDirectionalEval: true)
        };
    }
}