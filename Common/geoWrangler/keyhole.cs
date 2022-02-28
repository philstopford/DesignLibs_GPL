using ClipperLib2;
using System;
using System.Collections.Generic;
using System.Linq;
using utility;

namespace geoWrangler;

using Path = List<Point64>;
using Paths = List<List<Point64>>;
using PathsD = List<List<PointD>>;

public static partial class GeoWrangler
{
    // Sizing is used to define the keyhole width (default) and will be used for the sliver/gap removal.
    // Use of a custom value will cause headaches.
    public const double keyhole_sizing = 500;
    private const double keyhole_extension_default = 1.03;
    // Clipper 2 will not merge glancing contacts so we need to shift to provide an overlap. Value here is arbitrary and subject to more testing may change.
    public const double clipper_glancingContact_fudge = 1.1;

    public static Paths makeKeyHole(Paths outers, Paths cutters, double customSizing = 0, double extension = 0, double angularTolerance = 0)
    {
        return pMakeKeyHole(outers, cutters, customSizing: customSizing, extension: extension, angularTolerance: angularTolerance);
    }

    public static Paths makeKeyHole(Path source, double customSizing = 0, double extension = 0, double angularTolerance = 0)
    {
        return pMakeKeyHole(new Paths { source }, customSizing, extension, angularTolerance);
    }

    public static Paths makeKeyHole(Paths source, double customSizing = 0, double extension = 0, double angularTolerance = 0)
    {
        return pMakeKeyHole(source, customSizing: customSizing, extension: extension, angularTolerance: angularTolerance);
    }

    private static Paths pMakeKeyHole(Paths source, double customSizing = 0, double extension = 0, double angularTolerance = 0)
    {
        // Reconcile any overlapping geometry.
        Clipper c = new ();
        
        switch (source.Count)
        {
            case < 1:
                return source;
        }

        customSizing = customSizing switch
        {
            0 => keyhole_sizing,
            _ => customSizing
        };

        // Limit the offset used in the removal otherwise we can cause self-intersections that
        // result in lots of artifacts and trouble.
        Paths input = pRemoveFragments(source, customSizing, extension);

        switch (input.Count)
        {
            case < 2:
                return pClose(input);
        }

        Paths[] decomp = pGetDecomposed(input);
        Paths[] odecomp = new Paths[decomp.Length];
        for (int i = 0; i < decomp.Length; i++)
        {
            decomp[i] = pClose(decomp[i]);
            odecomp[i] = decomp[i].ToList();
        }
        // Outer areas
        List<double> outerAreas = new();

        foreach (Path p in decomp[(int)type.outer])
        {
            outerAreas.Add(ClipperFunc.Area(p));
        }
        
        // So here, things get annoying. We can have nested donuts, which means that we have outers fully covered by cutters (from the larger donut).
        // Unless we massage things, these get killed as the cutters are applied en-masse to outers in the keyholer.
        
        // First, we'll run the keyholer as usual.
        Paths ret = pMakeKeyHole(decomp[(int)type.outer], decomp[(int)type.cutter], customSizing, extension, angularTolerance);

        double origArea = 0;
        List<double> origAreas = new();
        foreach (Path t in sliverGapRemoval(source))
        {
            double tArea = ClipperFunc.Area(t);
            origAreas.Add(tArea);
            origArea += tArea;
        }

        double newArea = 0;
        List<double> newAreas = new();
        foreach (Path t in sliverGapRemoval(ret))
        {
            double tArea = ClipperFunc.Area(t);
            newAreas.Add(tArea);
            newArea += tArea;
        }
        
        // If we lost area, we probably had a cutter fully cover up one or more of our polygons.
        double lostArea = origArea - newArea;
        
        if (lostArea > 0)
        {
            // Track whether we are looking at the outer-most outer, and avoid touching it.
            bool bypassOuter = false;
            // We need to find out which cutters might have completely killed one or more outers and figure out a plan.
            for (int oIndex = 0; oIndex < odecomp[(int) type.outer].Count; oIndex++)
            {
                Path tOuter = odecomp[(int) type.outer][oIndex].ToList();
                double outerArea = ClipperFunc.Area(tOuter);
                if (!bypassOuter && (outerArea > lostArea))
                {
                    if (Math.Abs(outerArea - outerAreas.Max()) <= double.Epsilon)
                    {
                        bypassOuter = true;
                        continue;
                    }
                }
                Paths tCutters = new();
                // Do any cutters cover up our outer?
                for (int cIndex = 0; cIndex < odecomp[(int) type.cutter].Count; cIndex++)
                {
                    PolyTree pt = new();
                    c.Clear();
                    c.AddSubject(tOuter);
                    c.AddClip(odecomp[(int) type.cutter][cIndex].ToList());
                    c.Execute(ClipType.Difference, FillRule.EvenOdd, pt);
                    Paths test = ClipperFunc.PolyTreeToPaths(pt);
                    double area = 0;
                    foreach (Path t in test)
                    {
                        area += ClipperFunc.Area(t);
                    }

                    // If area is zero, the cutter fully covered the outer and we need to skip it.
                    if ((area != 0) && (area <= lostArea))
                    {
                        tCutters.Add(odecomp[(int) type.cutter][cIndex]);
                    }

                }

                Paths tOuters = new() {tOuter};

                Paths tRet = pMakeKeyHole(tOuters.ToList(), tCutters.ToList(), customSizing, extension,
                    angularTolerance);

                ret.AddRange(tRet.ToList());
            }
            
            // Screen ret for any duplicates and remove them.
            ret = pClockwiseAndReorder(ret.ToList());
            ret = removeDuplicatePaths(ret.ToList());
        }
        
        switch (ret.Count)
        {
            case 0:
                return source; // something blew up. Send back the original geometry.
            default:

                // Remove any overlapping duplicate polygons.
                PolyTree pt = new();
                c.AddSubject(ret);
                c.Execute(ClipType.Union, FillRule.Positive, pt);
                Paths cleaned = ClipperFunc.PolyTreeToPaths(pt);

                switch (cleaned.Count)
                {
                    default:
                        ret = pClose(cleaned);
                        break;
                    case 0:
                        ret = pClose(ret);
                        break;
                }

                return ret;
        }
    }

    private static Paths pMakeKeyHole(Paths outers, Paths cutters, double customSizing = 0, double extension = 0, double angularTolerance = 0)
    {
        customSizing = customSizing switch
        {
            0 => keyhole_sizing,
            _ => customSizing
        };

        try
        {
            bool outerOrient = ClipperFunc.IsClockwise(outers[0]);
            // Use raycaster to project from holes to outer, to try and find a keyhole path that is minimal length, and ideally orthogonal.
            foreach (Path t in cutters)
            {
                Paths extraCutters = new();
                Path projCheck = pStripTerminators(t, true);
                projCheck = pStripColinear(projCheck);
                //  Strip the terminator again to meet the requirements below.
                projCheck = pStripTerminators(projCheck, false);
                projCheck = pClockwise(projCheck);

                bool projectCorners = orthogonal(projCheck, angularTolerance);
                RayCast rc = new(t, outers, 1000000, invert: true, projectCorners: projectCorners);
                Paths clipped = rc.getClippedRays();
                // Need to find minimal length, ideally orthogonal.
                double minLength = -1;
                bool minLength_ortho = false;
                int cutPathIndex = -1;
                for (int r = 0; r < clipped.Count; r++)
                {
                    bool ray_isOrtho = clipped[r][0].X == clipped[r][1].X || clipped[r][0].Y == clipped[r][1].Y;

                    switch (minLength_ortho)
                    {
                        // Don't replace an orthogonal ray with a non-orthogonal ray.
                        case true when !ray_isOrtho:
                            continue;
                    }

                    double ray_length = Math.Abs(distanceBetweenPoints(clipped[r][0], clipped[r][1]));

                    switch (ray_length)
                    {
                        case <= double.Epsilon:
                            continue;
                    }

                    minLength_ortho = minLength_ortho || ray_isOrtho;

                    // First ray or a smaller distance causes us to make this the keyhole edge.
                    if (r != 0 && !(ray_length < minLength))
                    {
                        continue;
                    }

                    cutPathIndex = r;
                    minLength = ray_length;
                }

                if (cutPathIndex != -1)
                {
                    // Offset our cutter and assign to the clipping scenario.
                    Paths sPaths = pInflateEdge(clipped[cutPathIndex], customSizing);

                    extraCutters.AddRange(new Paths(sPaths));
                }
                
                // Review orientations.
                for (int p = 0; p < cutters.Count; p++)
                {
                    if (ClipperFunc.IsClockwise(cutters[p]))
                    {
                        cutters[p] = ClipperFunc.ReversePath(cutters[p]);
                    }
                }
                for (int p = 0; p < extraCutters.Count; p++)
                {
                    if (ClipperFunc.IsClockwise(extraCutters[p]))
                    {
                        extraCutters[p] = ClipperFunc.ReversePath(extraCutters[p]);
                    }
                }

                Clipper c = new();
                c.AddSubject(cutters);
                c.AddClip(extraCutters);

                Paths mergedCutters = new();
                PolyTree pt = new ();
                c.Execute(ClipType.Union, FillRule.EvenOdd, pt);
                mergedCutters = ClipperFunc.PolyTreeToPaths(pt);

                c.Clear();
                c.AddSubject(outers);
                c.AddClip(mergedCutters);

                // Reduce our geometry back to the simplest form.
                pt.Clear();
                Paths new_outers = new();
                c.Execute(ClipType.Difference, FillRule.EvenOdd, pt);//, PolyFillType.pftNonZero, PolyFillType.pftNegative);
                new_outers = ClipperFunc.PolyTreeToPaths(pt);

                outers.Clear();
                outers.AddRange(new_outers.Where(t1 => ClipperFunc.IsClockwise(t1) == outerOrient));
            }
            
            return pClockwiseAndReorder(outers);
        }
        catch (Exception)
        {
            throw new Exception("pMakeKeyHole error");
        }
    }

    private static Paths pInflateEdge(Path edge, double customSizing)
    {
        customSizing = customSizing switch
        {
            0 => keyhole_sizing,
            _ => customSizing
        };
        // Force clockwise, which should get us something consistent to work with.

        double dTmp0 = pDistanceBetweenPoints(new Point64(0, 0), edge[0]);
        double dTmp1 = pDistanceBetweenPoints(new Point64(0, 0), edge[1]);

        if (dTmp1 < dTmp0)
        {
            edge.Reverse();
        }

        // Get sorted out for dx, dy and normalization.
        double dx = edge[0].X - edge[1].X;
        double dy = edge[0].Y - edge[1].Y;

        double length = Math.Sqrt(Utils.myPow(dx, 2) + Utils.myPow(dy, 2));

        dx /= length;
        dy /= length;

        // Extend the line slightly.
        edge[0] = new Point64((long)(edge[0].X - Math.Abs(clipper_glancingContact_fudge * dx * keyhole_sizing)), (long)(edge[0].Y + Math.Abs(clipper_glancingContact_fudge * dy * keyhole_sizing)));
        edge[1] = new Point64((long)(edge[1].X + Math.Abs(clipper_glancingContact_fudge * dx * keyhole_sizing)), (long)(edge[1].Y - Math.Abs(clipper_glancingContact_fudge * dy * keyhole_sizing)));

        ClipperOffset co = new();
        co.AddPath(edge, JoinType.Square, EndType.Polygon);
        Paths sPaths = ClipperFunc.PathsFromPathsD(co.Execute(2 * customSizing));

        return sPaths;
    }

    public static Paths sliverGapRemoval(Path source, double customSizing = 0, double extension = 0, bool maySimplify = false)
    {
        return pSliverGapRemoval(new Paths { source }, customSizing, extension, maySimplify: maySimplify);
    }

    public static Paths sliverGapRemoval(Paths source, double customSizing = 0, double extension = 0, bool maySimplify = false)
    {
        return pSliverGapRemoval(source, customSizing, extension, maySimplify: maySimplify);
    }

    private static Paths pSliverGapRemoval(Paths source, double customSizing, double extension, bool maySimplify)
    {
        customSizing = customSizing switch
        {
            0 => keyhole_sizing,
            _ => customSizing
        };
        // Remove gaps, then remove slivers. Same process, different direction for sizing.
        Paths ret = pRemoveFragments(pRemoveFragments(source, customSizing, extension, maySimplify), -customSizing, extension, maySimplify: maySimplify);

        return ret.Count switch
        {
            0 =>
                // Probably not what we intended, so return the original as we assume the incoming geometry got crushed out of existence.
                source,
            _ => ret
        };
    }

    public static Paths gapRemoval(Path source, double customSizing = 0, double extension = 0, bool maySimplify = false)
    {
        return gapRemoval(new Paths { source }, customSizing, extension, maySimplify);
    }

    public static Paths gapRemoval(Paths source, double customSizing = 0, double extension = 0, bool maySimplify = false)
    {
        switch (source.Count)
        {
            case < 1:
                return source;
        }

        bool orig_orient_gw = isClockwise(source[0]);
        bool orig_orient_c = ClipperFunc.IsClockwise(source[0]);

        Paths ret = pRemoveFragments(source, customSizing, extension, maySimplify);

        switch (ret.Count)
        {
            case 0:
                return source; // something blew up. Send back the original geometry.
        }

        // Clean-up the geometry.
        PolyTree pt = new();
        Clipper c = new();
        c.AddSubject(ret);
        c.Execute(ClipType.Union, FillRule.EvenOdd, pt);
        ret = ClipperFunc.PolyTreeToPaths(pt);
        for (int p = 0; p < ret.Count; p++)
        {
            ret[p] = pStripColinear(ret[p]);
        }
        
        // Validate orientations.
        bool gR_orient_gw = isClockwise(ret[0]);
        bool gR_orient_c = ClipperFunc.IsClockwise(ret[0]);

        bool reverseNeeded = gR_orient_gw != orig_orient_gw || gR_orient_c != orig_orient_c;

        switch (reverseNeeded)
        {
            // Re-spin if needed. 1st path is always an outer, so we need to ensure the orientation is consistent.
            case true:
            {
                foreach (Path t in ret)
                {
                    t.Reverse();
                }

                break;
            }
        }

        // Clean-up stripped co-linear vertices; need to re-fragment.
        ret = pClose(ret);

        return ret;
    }

    public static Paths sliverRemoval(Paths source, double customSizing = 0, double extension = 0, bool maySimplify = false)
    {
        customSizing = customSizing switch
        {
            0 => keyhole_sizing,
            _ => customSizing
        };
        double oArea = source.Sum(t => ClipperFunc.Area(t));
        Paths ret = pRemoveFragments(source, -customSizing, extension, maySimplify: maySimplify);
        double nArea = ret.Sum(t => ClipperFunc.Area(t));

        return (Math.Abs(oArea) - Math.Abs(nArea)) switch
        {
            < 10000 => source,
            _ => ret
        };
    }

    // Positive incoming value removes gaps (keyholes); negative incoming value will remove slivers.
    private static Paths pRemoveFragments(Paths source, double customSizing, double extension, bool maySimplify = false, JoinType joinType = JoinType.Miter)
    {
        customSizing = customSizing switch
        {
            0 => keyhole_sizing,
            _ => customSizing
        };

        extension = extension switch
        {
            0 => keyhole_extension_default,
            _ => extension
        };

        // Used to try and avoid residual fragments; empirically derived.
        customSizing *= extension;
        
        ClipperOffset co = new();
        co.AddPaths(source, joinType, EndType.Polygon);
        Paths cGeometry = ClipperFunc.PathsFromPathsD(co.Execute(2 * customSizing));
        if (maySimplify)
        {
            Clipper c = new();
            PolyTree pt = new();
            c.AddSubject(cGeometry);
            c.Execute(ClipType.Union, FillRule.EvenOdd, pt);
            cGeometry = ClipperFunc.PolyTreeToPaths(pt);
        }
        co.Clear();
        co.AddPaths(cGeometry.ToList(), joinType, EndType.Polygon);
        cGeometry = ClipperFunc.PathsFromPathsD(co.Execute(-2 * customSizing)); // Size back to original dimensions
        if (maySimplify)
        {
            Clipper c = new();
            PolyTree pt = new();
            c.AddSubject(cGeometry);
            c.Execute(ClipType.Union, FillRule.EvenOdd, pt);
            cGeometry = ClipperFunc.PolyTreeToPaths(pt);
        }

        double newArea = cGeometry.Sum(t => ClipperFunc.Area(t));

        return Math.Abs(newArea) switch
        {
            <= double.Epsilon =>
                // We crushed our geometry, it seems.
                // This was probably not the plan, so send back the original geometry instead.
                source,
            _ => cGeometry
        };
    }

    private static Paths pRemoveFragments(Path source, double customSizing, bool maySimplify = false, JoinType joinType = JoinType.Miter)
    {
        customSizing = customSizing switch
        {
            0 => keyhole_sizing,
            _ => customSizing
        };
        double sourceArea = ClipperFunc.Area(source);

        ClipperOffset co = new();
        co.AddPath(source, joinType, EndType.Polygon);
        Paths cGeometry = ClipperFunc.PathsFromPathsD(co.Execute(2 * customSizing));
        if (maySimplify)
        {
            Clipper c = new();
            PolyTree pt = new();
            c.AddSubject(cGeometry);
            c.Execute(ClipType.Union, FillRule.EvenOdd, pt);
            cGeometry = ClipperFunc.PolyTreeToPaths(pt);
        }
        co.Clear();
        co.AddPaths(cGeometry.ToList(), joinType, EndType.Polygon);
        cGeometry = ClipperFunc.PathsFromPathsD(co.Execute(-2 * customSizing)); // Size back to original dimensions
        if (maySimplify)
        {
            Clipper c = new();
            PolyTree pt = new();
            c.AddSubject(cGeometry);
            c.Execute(ClipType.Union, FillRule.EvenOdd, pt);
            cGeometry = ClipperFunc.PolyTreeToPaths(pt);
        }

        double newArea = 0;
        foreach (Path t in cGeometry)
        {
            newArea += ClipperFunc.Area(t);
        }

        switch (Math.Abs(newArea))
        {
            case <= double.Epsilon:
                // We crushed our geometry, it seems.
                // This was probably not the plan, so send back the original geometry instead.
                cGeometry.Clear();
                cGeometry.Add(source);
                break;
            default:
            {
                switch (sourceArea)
                {
                    // Do we need to flip the direction to match the original orientation?
                    case < 0 when newArea > 0:
                    case > 0 when newArea < 0:
                    {
                        // Multi-path handling gets interesting. The first path is assumed to be the outer. Let's compare that with the original geometry. If the orientation is different, reverse the full set.
                        bool orientation = ClipperFunc.IsClockwise(source);
                        bool origCG0_o = ClipperFunc.IsClockwise(cGeometry[0]);
                        foreach (Path t in cGeometry.Where(t => origCG0_o != orientation))
                        {
                            t.Reverse();
                        }

                        break;
                    }
                }

                break;
            }
        }

        return cGeometry;
    }
}