using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Clipper2Lib;

namespace geoWrangler;

using Path = List<Point64>;
using Paths = List<List<Point64>>;

public static partial class GeoWrangler
{
    public static Paths customBoolean(int firstLayerOperator, Paths firstLayer, int secondLayerOperator, Paths secondLayer, int booleanFlag, double resolution, double extension, int scaling)
    {
        Paths ret = pCustomBoolean(firstLayerOperator, firstLayer, secondLayerOperator, secondLayer, booleanFlag, resolution, extension, scaling);

        return ret;
    }

    private static Paths pCustomBoolean(int firstLayerOperator, Paths firstLayer, int secondLayerOperator, Paths secondLayer, int booleanFlag, double resolution, double extension, int scaling)
    {
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
        firstLayer = sliverGapRemoval(firstLayer, customSizing:0.5*keyhole_sizing);
        secondLayer = sliverGapRemoval(secondLayer, customSizing:0.5*keyhole_sizing);
        // Clipper strips terminating points, so force closed.
        firstLayer = close(firstLayer);
        secondLayer = close(secondLayer);
        Paths ret = pLayerBoolean(firstLayerOperator, firstLayer, secondLayerOperator, secondLayer, booleanFlag, preserveColinear: false);
        
        // Secondary clean-up of the result. This seems to be needed, so retained for now.
        ret = gapRemoval(ret, customSizing:0.5*keyhole_sizing,extension: extension).ToList();

        bool holes = false;

        foreach (Path t in ret)
        {
            holes = !ClipperFunc.IsClockwise(t);
            bool gwHoles = !isClockwise(t);
            if (holes != gwHoles)
            {
            }
            if (holes)
            {
                break;
            }
        }

        // Apply the keyholing and rationalize.
        if (holes)
        {
            Fragmenter f = new(resolution * scaling);
            ret = f.fragmentPaths(ret);
            Paths merged = makeKeyHole(ret, extension:extension);

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
            Clipper c1 = new() {PreserveCollinear = true };
            c1.AddSubject(merged);
            c1.Execute(ClipType.Union, FillRule.EvenOdd, ret);
            ret = reOrderXY(ret);
            ret = stripColinear(ret, 1.0);
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
            return stripColinear(ret, 1.0);
        }

        Rect64 bounds = ClipperFunc.GetBounds(ret);

        Path bound = new()
        {
            new Point64(bounds.left, bounds.bottom),
            new Point64(bounds.left, bounds.top),
            new Point64(bounds.right, bounds.top),
            new Point64(bounds.right, bounds.bottom),
            new Point64(bounds.left, bounds.bottom)
        };

        Clipper c = new() {PreserveCollinear = false};

        c.AddSubject(ret);
        c.AddClip(bound);

        Paths simple = new();
        c.Execute(ClipType.Intersection, FillRule.EvenOdd, simple);
        ret = reOrderXY(simple);

        return clockwiseAndReorderXY(simple);
    }

    public static Paths LayerBoolean(int firstLayerOperator, Paths firstLayerPaths, int secondLayerOperator,
        Paths secondLayerPaths, int booleanFlag, bool preserveColinear)
    {
        return pLayerBoolean(firstLayerOperator, firstLayerPaths, secondLayerOperator,
            secondLayerPaths, booleanFlag, preserveColinear);
    }
    private static Paths pLayerBoolean(int firstLayerOperator, Paths firstLayerPaths, int secondLayerOperator, Paths secondLayerPaths, int booleanFlag, bool preserveColinear)
    {
        if (firstLayerOperator == 1) // NOT layer handling
        {
            try
            {
                firstLayerPaths = invertTone(firstLayerPaths, preserveColinear: preserveColinear).ToList();
            }
            catch (Exception)
            {
                // Something blew up.
            }
            firstLayerPaths[0] = close(firstLayerPaths[0]);
        }

        if (secondLayerOperator == 1) // NOT layer handling
        {
            try
            {
                secondLayerPaths = invertTone(secondLayerPaths, preserveColinear: preserveColinear).ToList();
            }
            catch (Exception)
            {
                // Something blew up.
            }
            secondLayerPaths[0] = close(secondLayerPaths[0]);
        }

        if (firstLayerPaths[0].Count <= 1)
        {
            return secondLayerPaths.ToList();
        }
        return secondLayerPaths[0].Count <= 1 ? firstLayerPaths.ToList() : pLayerBoolean(firstLayerPaths, secondLayerPaths, booleanFlag, preserveColinear: preserveColinear);
    }

    public static Paths LayerBoolean(Paths firstPaths, Paths secondPaths, int booleanFlag, bool preserveColinear = true)
    {
        return pLayerBoolean(firstPaths, secondPaths, booleanFlag, preserveColinear);

    }
    private static Paths pLayerBoolean(Paths firstPaths, Paths secondPaths, int booleanFlag, bool preserveColinear = true)
    {
        string booleanType = "AND";
        if (booleanFlag == 1)
        {
            booleanType = "OR";
        }

        // important - if we don't do this, we lose the fragmentation on straight edges.
        Clipper c = new() {PreserveCollinear = preserveColinear};

        c.AddSubject(firstPaths);
        c.AddClip(secondPaths);

        Paths outputPoints = new();

        switch (booleanType)
        {
            case "AND":
                c.Execute(ClipType.Intersection, FillRule.EvenOdd, outputPoints);
                break;
            case "OR":
                c.Execute(ClipType.Union, FillRule.EvenOdd, outputPoints);
                break;
        }

        outputPoints = reOrderXY(outputPoints);
        
        return outputPoints; // Return our first list of points as the result of the boolean.
    }

}