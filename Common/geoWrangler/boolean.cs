using System;
using System.Threading.Tasks;
using Clipper2Lib;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public enum LayerFlag
    {
        none,
        NOT
    }

    public enum booleanOperation
    {
        AND,
        OR
    }

    public static PathsD customBoolean(int firstLayerOperator, PathsD firstLayer, int secondLayerOperator, PathsD secondLayer, int booleanFlag, double resolution, double extension)
    {
        PathsD ret = pCustomBoolean(firstLayerOperator, firstLayer, secondLayerOperator, secondLayer, booleanFlag, resolution, extension);

        return ret;
    }

    private static PathsD pCustomBoolean(int firstLayerOperator, PathsD firstLayer_, int secondLayerOperator, PathsD secondLayer_, int booleanFlag, double resolution, double extension)
    {
        PathsD firstLayer = new (firstLayer_);
        PathsD secondLayer = new (secondLayer_);
        Fragmenter f = new(resolution);

        // In principle, 'rigorous' handling is only needed where the cutter is fully enclosed by the subject polygon.
        // The challenge is to know whether this is the case or not.
        // Possibility would be an intersection test and a vertex count and location comparison from before and after, to see whether anything changed.
        bool rigorous = enclosed(firstLayer, secondLayer); // this is not a strict check because the enclosure can exist either way for this situation.
        
        // Need a secondary check because keyholed geometry could be problematic.
        // Both paths will be reviewed; first one to have a keyhole will trigger the rigorous process.
        if (!rigorous)
        {
            try
            {
                rigorous = enclosed(firstLayer, customSizing: 1, extension: extension, strict: true); // force a strict check.

                if (!rigorous)
                {
                    // Need a further check because keyholed geometry in B could be problematic.
                    rigorous = enclosed(secondLayer, customSizing: 1, extension: extension, strict: true); // force a strict check.
                }
            }
            catch (Exception)
            {
                // No big deal - carry on.
            }
        }

        // Squash incoming artifacts to allow boolean to resolve holes vs outers.
        // The 0.5 factor here allows for both sides of a keyhole cut to be moved to touch, thus merging.
        double scaling_for_diagonal_keyholes = 2 * keyhole_sizing;
        firstLayer = sliverGapRemoval(firstLayer, customSizing:scaling_for_diagonal_keyholes);
        for (int i = 0;  i < firstLayer.Count; i++)
        {
            firstLayer[i] = stripCollinear(firstLayer[i], 1);
        }
        firstLayer = f.fragmentPaths(firstLayer);
        secondLayer = sliverGapRemoval(secondLayer, customSizing:scaling_for_diagonal_keyholes);
        for (int i = 0;  i < secondLayer.Count; i++)
        {
            secondLayer[i] = stripCollinear(secondLayer[i], 1);
        }
        secondLayer = f.fragmentPaths(secondLayer);
        // Clipper strips terminating points, so force closed.
        firstLayer = close(firstLayer);
        secondLayer = close(secondLayer);
        PathsD ret = pLayerBoolean(firstLayerOperator, firstLayer, secondLayerOperator, secondLayer, booleanFlag, preserveCollinear: false);
        
        // Secondary clean-up of the result. This seems to be needed, so retained for now.
        ret = new PathsD(gapRemoval(ret, customSizing:scaling_for_diagonal_keyholes, extension: extension));

        bool holes = false;

        foreach (PathD t in ret)
        {
            holes = !Clipper.IsPositive(t); // reports false for outers
            if (holes)
            {
                break;
            }
        }

        // Apply the keyholing and rationalize.
        if (holes)
        {
            ret = f.fragmentPaths(ret);
            PathsD merged = makeKeyHole(ret, reverseEval:false, biDirectionalEval:true, extension:extension);

            int count = merged.Count;
#if !GWSINGLETHREADED
            Parallel.For(0, count, i =>
#else
                for (int i = 0; i < count; i++)
#endif
                {
                    merged[i] = clockwise(merged[i]);
                }
#if !GWSINGLETHREADED
            );
#endif
            // Squash any accidental keyholes - not ideal, but best option found so far.
            ClipperD c1 = new(Constants.roundingDecimalPrecision) {PreserveCollinear = true };
            c1.AddSubject(merged);
            c1.Execute(ClipType.Union, FillRule.EvenOdd, ret);
            ret = reOrderXY(ret);
            ret = stripCollinear(ret);
        }

        ret = sliverRemoval(ret, extension: extension); // experimental to try and remove any slivers.

        if (rigorous && !holes)
        {
            int count = ret.Count;
#if !GWSINGLETHREADED
            Parallel.For(0, count, i =>
#else
                for (int i = 0; i < count; i++)
#endif
                {
                    ret[i] = clockwise(ret[i]);
                    ret[i] = close(ret[i]);
                }
#if !GWSINGLETHREADED
            );
#endif
            // Return here because the attempt to rationalize the geometry below also screws things up, it seems.
            return stripCollinear(ret);
        }

        RectD bounds = Clipper.GetBounds(ret);

        PathD bound =
        [
            new PointD(bounds.left, bounds.bottom),
            new PointD(bounds.left, bounds.top),
            new PointD(bounds.right, bounds.top),
            new PointD(bounds.right, bounds.bottom),
            new PointD(bounds.left, bounds.bottom)
        ];

        ClipperD c = new(Constants.roundingDecimalPrecision) {PreserveCollinear = false};

        c.AddSubject(ret);
        c.AddClip(bound);

        PathsD simple = [];
        c.Execute(ClipType.Intersection, FillRule.EvenOdd, simple);
        // ret = reOrderXY(simple);

        return clockwiseAndReorderXY(simple);
    }

    public static PathsD LayerBoolean(int firstLayerOperator, PathsD firstLayerPaths, int secondLayerOperator,
        PathsD secondLayerPaths, int booleanFlag, bool preserveCollinear)
    {
        return pLayerBoolean(firstLayerOperator, firstLayerPaths, secondLayerOperator,
            secondLayerPaths, booleanFlag, preserveCollinear);
    }
    private static PathsD pLayerBoolean(int firstLayerOperator, PathsD firstLayerPaths_, int secondLayerOperator, PathsD secondLayerPaths_, int booleanFlag, bool preserveCollinear)
    {
        PathsD firstLayerPaths = new(firstLayerPaths_);
        if (firstLayerOperator == (int)LayerFlag.NOT)
        {
            try
            {
                firstLayerPaths = new PathsD(invertTone(firstLayerPaths_, preserveCollinear: preserveCollinear));
            }
            catch (Exception)
            {
                // Something blew up.
            }
            firstLayerPaths[0] = close(firstLayerPaths[0]);
        }

        PathsD secondLayerPaths = new(secondLayerPaths_);
        if (secondLayerOperator == (int)LayerFlag.NOT)
        {
            try
            {
                secondLayerPaths = new PathsD(invertTone(secondLayerPaths, preserveCollinear: preserveCollinear));
            }
            catch (Exception)
            {
                // Something blew up.
            }
            secondLayerPaths[0] = close(secondLayerPaths[0]);
        }

        if (firstLayerPaths[0].Count <= 1)
        {
            return new PathsD(secondLayerPaths);
        }
        return secondLayerPaths[0].Count <= 1 ? new PathsD(firstLayerPaths) : pLayerBoolean(firstLayerPaths, secondLayerPaths, booleanFlag, preserveCollinear: preserveCollinear);
    }

    public static PathsD LayerBoolean(PathsD firstPaths, PathsD secondPaths, int booleanFlag, bool preserveCollinear = true)
    {
        return pLayerBoolean(firstPaths, secondPaths, booleanFlag, preserveCollinear);

    }
    private static PathsD pLayerBoolean(PathsD firstPaths, PathsD secondPaths, int booleanFlag, bool preserveCollinear = true)
    {
        // important - if we don't do this, we lose the fragmentation on straight edges.
        ClipperD c = new(Constants.roundingDecimalPrecision) {PreserveCollinear = preserveCollinear};

        c.AddSubject(firstPaths);
        c.AddClip(secondPaths);

        PathsD outputPoints = [];

        switch (booleanFlag)
        {
            case (int)booleanOperation.AND:
                c.Execute(ClipType.Intersection, FillRule.EvenOdd, outputPoints);
                break;
            case (int)booleanOperation.OR:
                c.Execute(ClipType.Union, FillRule.EvenOdd, outputPoints);
                break;
        }

        outputPoints = reOrderXY(outputPoints);

        return outputPoints; // Return our first list of points as the result of the boolean.
    }

}