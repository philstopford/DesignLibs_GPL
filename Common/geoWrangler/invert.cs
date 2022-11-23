using Clipper2Lib;

namespace geoWrangler;

public static partial class GeoWrangler
{
    // Bounds will force the negation to only work against the extents of the shape. Useful for capturing islands in negative tone.
    public static PathsD invertTone(PathD source, bool preserveColinear, bool useTriangulation = false, bool useBounds = false)
    {
        return pInvertTone(new PathsD { source }, preserveColinear: preserveColinear, useTriangulation: useTriangulation, useBounds: useBounds);
    }

    public static PathsD invertTone(PathsD source, bool preserveColinear, bool useTriangulation = false, bool useBounds = false)
    {
        return pInvertTone(source, preserveColinear: preserveColinear, useTriangulation: useTriangulation, useBounds: useBounds);
    }
    
    private static PathsD pInvertTone(PathsD sourcePaths, bool preserveColinear, bool useTriangulation, bool useBounds)
    {
        PathD firstLayerBP = new();
        switch (useBounds)
        {
            case false:
                firstLayerBP.Add(new (-int.MaxValue, -int.MaxValue));
                firstLayerBP.Add(new (-int.MaxValue, int.MaxValue));
                firstLayerBP.Add(new (int.MaxValue, int.MaxValue));
                firstLayerBP.Add(new (int.MaxValue, -int.MaxValue));
                firstLayerBP.Add(new (-int.MaxValue, -int.MaxValue));
                break;
            default:
            {
                RectD bounds = Clipper.GetBounds(sourcePaths);
                firstLayerBP.Add(new (bounds.left, bounds.bottom));
                firstLayerBP.Add(new (bounds.left, bounds.top));
                firstLayerBP.Add(new (bounds.right, bounds.top));
                firstLayerBP.Add(new (bounds.right, bounds.bottom));
                firstLayerBP.Add(new (bounds.left, bounds.bottom));
                break;
            }
        }

        ClipperD c = new() {PreserveCollinear = preserveColinear};

        c.AddSubject(firstLayerBP);
        // Add hole polygons from our paths
        c.AddClip(sourcePaths);

        PathsD cutters = new();
        c.Execute(ClipType.Difference, FillRule.EvenOdd, cutters);

        cutters = pReorderXY(cutters);
        
        switch (useTriangulation)
        {
            case false:
                return cutters;
            default:
                return pMakeKeyHole(sourcePaths, reverseEval:false, biDirectionalEval:true);
        }
    }
}