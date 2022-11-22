using Clipper2Lib;
using System;
using System.Collections.Generic;
using System.Linq;

namespace geoWrangler;

public static partial class GeoWrangler
{
    // Sizing is used to define the keyhole width (default) and will be used for the sliver/gap removal.
    // Use of a custom value will cause headaches.
    public const double keyhole_sizing = 500;
    private const double keyhole_extension_default = 1.00;

    public static PathsD makeKeyHole(PathsD outers, PathsD cutters, bool reverseEval, bool biDirectionalEval, RayCast.inversionMode invert = RayCast.inversionMode.x, double customSizing = 0, double extension = 0, double angularTolerance = 0)
    {
        return pMakeKeyHole(outers, cutters, reverseEval, biDirectionalEval, invert, customSizing: customSizing, extension: extension, angularTolerance: angularTolerance);
    }

    public static PathsD makeKeyHole(PathD source, bool reverseEval, bool biDirectionalEval, RayCast.inversionMode invert = RayCast.inversionMode.x, double customSizing = 0, double extension = 0, double angularTolerance = 0)
    {
        return pMakeKeyHole(new PathsD { source }, reverseEval, biDirectionalEval, invert, customSizing, extension, angularTolerance);
    }

    public static PathsD makeKeyHole(PathsD source, bool reverseEval, bool biDirectionalEval, RayCast.inversionMode invert = RayCast.inversionMode.x, double customSizing = 0, double extension = 0, double angularTolerance = 0)
    {
        return pMakeKeyHole(source, reverseEval, biDirectionalEval, invert, customSizing: customSizing, extension: extension, angularTolerance: angularTolerance);
    }

    private static PathsD pMakeKeyHole(PathsD source, bool reverseEval, bool biDirectionalEval, RayCast.inversionMode invert = RayCast.inversionMode.x, double customSizing = 0, double extension = 0, double angularTolerance = 0)
    {
        // Reconcile any overlapping geometry.
        ClipperD c = new () {PreserveCollinear = true};
        
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
        PathsD input = pRemoveFragments(source, customSizing, extension);
        input = pStripColinear(input);

        switch (input.Count)
        {
            case < 2:
                return pClose(input);
        }

        PathsD[] decomp = pGetDecomposed(input);
        PathsD[] odecomp = new PathsD[decomp.Length];
        for (int i = 0; i < decomp.Length; i++)
        {
            decomp[i] = pClose(decomp[i]);
            odecomp[i] = new(decomp[i]);
        }
        // Outer areas
        List<double> outerAreas = decomp[(int)type.outer].Select(Clipper.Area).ToList();

        // So here, things get annoying. We can have nested donuts, which means that we have outers fully covered by cutters (from the larger donut).
        // Unless we massage things, these get killed as the cutters are applied en-masse to outers in the keyholer.
        
        // First, we'll run the keyholer as usual.
        PathsD ret = pMakeKeyHole(decomp[(int)type.outer], decomp[(int)type.cutter], reverseEval, biDirectionalEval, invert, customSizing, extension, angularTolerance);
        
        double origArea = 0;
        List<double> origAreas = new();
        foreach (PathD t in sliverGapRemoval(source))
        {
            double tArea = Clipper.Area(t);
            origAreas.Add(tArea);
            origArea += tArea;
        }

        double newArea = 0;
        List<double> newAreas = new();
        foreach (PathD t in sliverGapRemoval(ret))
        {
            double tArea = Clipper.Area(t);
            newAreas.Add(tArea);
            newArea += tArea;
        }

        if (newArea < 0)
        {
            newArea = -newArea;
        }
        
        // If we lost area, we probably had a cutter fully cover up one or more of our polygons.
        double lostArea = Math.Abs(origArea - newArea);
        
        if (lostArea > 0)
        {
            // We need to find out which cutters might have completely killed one or more outers and figure out a plan.
            for (int oIndex = 0; oIndex < odecomp[(int) type.outer].Count; oIndex++)
            {
                PathD tOuter = odecomp[(int) type.outer][oIndex];
                double outerArea = Clipper.Area(tOuter);
                if (outerArea > lostArea)
                {
                    if (Math.Abs(outerArea - outerAreas.Max()) <= double.Epsilon)
                    {
                        continue;
                    }
                }
                PathsD tCutters = new();
                // Do any cutters cover up our outer?
                for (int cIndex = 0; cIndex < odecomp[(int) type.cutter].Count; cIndex++)
                {
                    PathsD test = new();
                    c.Clear();
                    c.AddSubject(tOuter);
                    c.AddClip(odecomp[(int) type.cutter][cIndex]);
                    c.Execute(ClipType.Difference, FillRule.EvenOdd, test);
                    test = pReorderXY(test);
                    double area = 0;
                    foreach (PathD t in test)
                    {
                        area += Clipper.Area(t);
                    }

                    // If area is zero, the cutter fully covered the outer and we need to skip it.
                    if ((area != 0) && (area <= lostArea))
                    {
                        tCutters.Add(odecomp[(int) type.cutter][cIndex]);
                    }

                }

                PathsD tOuters = new() {tOuter};

                PathsD tRet = pMakeKeyHole(tCutters, tOuters, reverseEval, biDirectionalEval, invert, customSizing, extension,
                    angularTolerance);

                ret.AddRange(tRet);
            }
            
            // Screen ret for any duplicates and remove them.
            ret = pClockwiseAndReorderXY(ret);
            ret = removeDuplicatePaths(ret);
        }
        
        switch (ret.Count)
        {
            case 0:
                return source; // something blew up. Send back the original geometry.
            default:

                // Remove any overlapping duplicate polygons.
                ret = pClose(ret);

                return ret;
        }
    }

    private static PathsD pClipRays(PathsD outers, PathD cutter, bool reverseEval, bool biDirectionalEval, RayCast.inversionMode invert = RayCast.inversionMode.x, double angularTolerance = 0)
    {
        // Needed due to ordering sequence from ClipperLib2. Have to clean, re-order and re-close to make the geometry work for the raycaster.
        PathD t = new (cutter);
        t = pStripTerminators(t, false); // chop the terminator off to ensure re-ordering doesn't yield a zero-length segment.
        t = pClockwiseAndReorderYX(t); // speculative : might need to be XY, but picked YX for now.
        bool projectCorners = pOrthogonal(t, angularTolerance);
        t = pClose(t); // re-close to make the raycaster happy.

        RayCast rc;
        PathsD clipped = new();

        /*
         * Some explanation is required for the below. Project corners will only shoot a single ray out from each corner. This makes the evaluation
         * sensitive to the direction of travel around the shape. This can lead to missed candidates for ray insertion. To counter this, both directions
         * of travel must be evaluated.
         * The cost of this is somewhat mitigated by heavy use of multithreading in the raycaster, but the desire for a robust calculation makes the cost
         * worth the effort.
         *
         * However, in some cases, the reverse walk is not desired (e.g. for rectangular decomposition in Quilt, where single emission is preferred).
         */
        if (biDirectionalEval)
        {
            // Reverse walk in case there is a better option walking the geometry in the other direction.
            rc = new(t, outers, 1000000, invert: invert, projectCorners: projectCorners);
            clipped = rc.getClippedRays();
        }

        if (!reverseEval || biDirectionalEval)
        {
            t.Reverse();
        }
        rc = new(t, outers, 1000000, invert: invert, projectCorners: projectCorners);
        clipped.AddRange(rc.getClippedRays());

        return clipped;
    }

    private static PathD pFindInsertionCandidate(PathsD edges, bool orthogonalInput)
    {
        PathD ret = null;
        // Need to find minimal length, ideally orthogonal.
        double minLength = -1;
        bool minLength_ortho = false;
        int cutPathIndex = -1;
        for (int r = 0; r < edges.Count; r++)
        {
            // If the input geometry is purely orthogonal, this is a strong preference.
            // Otherwise, we don't care and will take the minimum internal edge.
            if (orthogonalInput)
            {
                bool ray_isOrtho = edges[r][0].x == edges[r][1].x || edges[r][0].y == edges[r][1].y;
                if (minLength_ortho && !ray_isOrtho)
                {
                    // Don't replace an orthogonal ray with a non-orthogonal ray.
                    continue;
                }
                        
                // Update our tracking with the ray ortho flag.
                minLength_ortho = ray_isOrtho;
            }

            double ray_length = Math.Abs(distanceBetweenPoints(edges[r][0], edges[r][1]));

            switch (ray_length)
            {
                case <= double.Epsilon:
                    continue;
            }

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
            ret = edges[cutPathIndex];
        }

        return ret;
    }

    private static PathsD pCutKeyHole(PathsD outers, PathsD cutters, PathsD extraCutters)
    {
        ClipperD c = new();
        c.AddSubject(cutters);
        c.AddClip(extraCutters);

        PathsD mergedCutters = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, mergedCutters);

        mergedCutters = pReorderXY(mergedCutters);

        c.Clear();
        c.AddSubject(outers);
        c.AddClip(mergedCutters);

        // Reduce our geometry back to the simplest form.
        PathsD new_outers = new();
        c.Execute(ClipType.Difference, FillRule.EvenOdd, new_outers);//, PolyFillType.pftNonZero, PolyFillType.pftNegative);

        new_outers = pReorderXY(new_outers);

        return new_outers;
    }
    
    private static PathsD pMakeKeyHole(PathsD outers, PathsD cutters, bool reverseEval, bool biDirectionalEval, RayCast.inversionMode invert = RayCast.inversionMode.x, double customSizing = 0, double extension = 0, double angularTolerance = 0)
    {
        customSizing = customSizing switch
        {
            0 => keyhole_sizing,
            _ => customSizing
        };

        try
        {
            bool outerOrient = Clipper.IsPositive(outers[0]);
            bool orthogonalInput = pOrthogonal(outers, 0);
            orthogonalInput = orthogonalInput && pOrthogonal(cutters, 0);
            // Review orientations.
            for (int p = 0; p < cutters.Count; p++)
            {
                if (Clipper.IsPositive(cutters[p]))
                {
                    cutters[p] = Clipper.ReversePath(cutters[p]);
                }
            }
            // Use raycaster to project from holes to outer, to try and find a keyhole path that is minimal length, and ideally orthogonal.
            foreach (PathD t1 in cutters)
            {
                PathsD clipped = pClipRays(outers, t1, reverseEval, biDirectionalEval, invert, angularTolerance);

                PathD insertionCandidate = pFindInsertionCandidate(clipped, orthogonalInput);

                PathsD extraCutters = new();
                if (insertionCandidate != null)
                {
                    // Offset our cutter and assign to the clipping scenario.
                    PathsD sPaths = pInflateEdge(insertionCandidate, customSizing);

                    extraCutters.AddRange(new PathsD (sPaths));
                }

                for (int p = 0; p < extraCutters.Count; p++)
                {
                    if (Clipper.IsPositive(extraCutters[p]))
                    {
                        extraCutters[p] = Clipper.ReversePath(extraCutters[p]);
                    }
                }

                PathsD new_outers = pCutKeyHole(outers, cutters, extraCutters);

                outers.Clear();
                outers.AddRange(new_outers.Where(t1_ => Clipper.IsPositive(t1_) == outerOrient));
            }
            
            return pClockwiseAndReorderXY(outers);
        }
        catch (Exception)
        {
            throw new Exception("pMakeKeyHole error");
        }
    }

    private static PathsD pInflateEdge(PathD edge, double customSizing)
    {
        customSizing = customSizing switch
        {
            0 => keyhole_sizing,
            _ => customSizing
        };
        // Force clockwise, which should get us something consistent to work with.
        edge = pClockwise(edge);

        // edge = pExtendEdge(edge, keyhole_sizing);

        // Need to workaround missing PathD support in ClipperOffset...
        double scalar = 10000;
        Path64 rescaledSource = _pPath64FromPathD(edge, scalar);

        ClipperOffset co = new() {PreserveCollinear = true};
        co.AddPath(rescaledSource, JoinType.Miter, EndType.Square);
        
        Paths64 tmp = co.Execute(2 * customSizing);

        PathsD sPaths = _pPathsDFromPaths64(tmp, scalar);

        return pReorderXY(sPaths);
    }

    public static PathsD sliverGapRemoval(PathD source, double customSizing = 0, double extension = 0, bool maySimplify = false)
    {
        return pSliverGapRemoval(new () { source }, customSizing, extension, maySimplify: maySimplify);
    }

    public static PathsD sliverGapRemoval(PathsD source, double customSizing = 0, double extension = 0, bool maySimplify = false)
    {
        return pSliverGapRemoval(source, customSizing, extension, maySimplify: maySimplify);
    }

    private static PathsD pSliverGapRemoval(PathsD source, double customSizing, double extension, bool maySimplify)
    {
        customSizing = customSizing switch
        {
            0 => keyhole_sizing,
            _ => customSizing
        };
        // Remove gaps, then remove slivers. Same process, different direction for sizing.
        PathsD ret = pRemoveFragments(pRemoveFragments(source, customSizing, extension, maySimplify), -customSizing, extension, maySimplify: maySimplify);

        return ret.Count switch
        {
            0 =>
                // Probably not what we intended, so return the original as we assume the incoming geometry got crushed out of existence.
                source,
            _ => ret
        };
    }

    public static PathsD gapRemoval(PathD source, double customSizing = 0, double extension = 0, bool maySimplify = false)
    {
        return gapRemoval(new PathsD { source }, customSizing, extension, maySimplify);
    }

    public static PathsD gapRemoval(PathsD source, double customSizing = 0, double extension = 0, bool maySimplify = false)
    {
        switch (source.Count)
        {
            case < 1:
                return source;
        }

        bool orig_orient_gw = isClockwise(source[0]);
        bool orig_orient_c = Clipper.IsPositive(source[0]);

        PathsD ret = pRemoveFragments(source, customSizing, extension, maySimplify);

        switch (ret.Count)
        {
            case 0:
                return source; // something blew up. Send back the original geometry.
        }

        // Clean-up the geometry.
        ClipperD c = new();
        c.AddSubject(ret);
        c.Execute(ClipType.Union, FillRule.EvenOdd, ret);

        ret = pReorderXY(ret);

        // ret = stripColinear(ret);

        // Validate orientations.
        bool gR_orient_gw = isClockwise(ret[0]);
        bool gR_orient_c = Clipper.IsPositive(ret[0]);

        bool reverseNeeded = gR_orient_gw != orig_orient_gw || gR_orient_c != orig_orient_c;

        switch (reverseNeeded)
        {
            // Re-spin if needed. 1st path is always an outer, so we need to ensure the orientation is consistent.
            case true:
            {
                foreach (PathD t in ret)
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

    public static PathsD sliverRemoval(PathsD source, double customSizing = 0, double extension = 0, bool maySimplify = false)
    {
        customSizing = customSizing switch
        {
            0 => keyhole_sizing,
            _ => customSizing
        };
        double oArea = source.Sum(t => Clipper.Area(t));
        PathsD ret = pRemoveFragments(source, -customSizing, extension, maySimplify: maySimplify);
        double nArea = ret.Sum(t => Clipper.Area(t));

        return (Math.Abs(oArea) - Math.Abs(nArea)) switch
        {
            < 10000 => source,
            _ => ret
        };
    }

    // Positive incoming value removes gaps (keyholes); negative incoming value will remove slivers.
    private static PathsD pRemoveFragments(PathsD source, double customSizing, double extension, bool maySimplify = false, JoinType joinType = JoinType.Miter)
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
        
        // Need to workaround missing PathD support in ClipperOffset...
        double scalar = 10000;
        Paths64 rescaledSource = _pPaths64FromPathsD(source, scalar);

        ClipperOffset co = new() {PreserveCollinear = true};
        co.AddPaths(rescaledSource, joinType, EndType.Polygon);
        Paths64 tmp = co.Execute(customSizing);
        co.Clear();
        co.AddPaths(new(tmp), joinType, EndType.Polygon);
        tmp.Clear();
        tmp = co.Execute(-customSizing); // Size back to original dimensions

        PathsD cGeometry = _pPathsDFromPaths64(tmp, scalar);
        
        cGeometry = pReorderXY(cGeometry);
        
        if (maySimplify)
        {
            cGeometry = stripColinear(cGeometry);
        }

        double newArea = cGeometry.Sum(t => Clipper.Area(t));

        return Math.Abs(newArea) switch
        {
            <= double.Epsilon =>
                // We crushed our geometry, it seems.
                // This was probably not the plan, so send back the original geometry instead.
                source,
            _ => pReorderXY(cGeometry)
        };
    }

    private static PathsD pRemoveFragments(PathD source, double customSizing, bool maySimplify = false, JoinType joinType = JoinType.Miter)
    {
        customSizing = customSizing switch
        {
            0 => keyhole_sizing,
            _ => customSizing
        };
        double sourceArea = Clipper.Area(source);

        // Need to workaround missing PathD support in ClipperOffset...
        double scalar = 10000;
        Path64 rescaledSource = _pPath64FromPathD(source, scalar);

        ClipperOffset co = new() {PreserveCollinear = true};
        co.AddPath(rescaledSource, joinType, EndType.Polygon);
        Paths64 tmp = co.Execute(customSizing);
        co.Clear();
        co.AddPaths(new (tmp), joinType, EndType.Polygon);
        tmp.Clear();
        tmp = co.Execute(-customSizing); // Size back to original dimensions

        PathsD cGeometry = _pPathsDFromPaths64(tmp, scalar);

        cGeometry = pReorderXY(cGeometry);
        
        if (maySimplify)
        {
            cGeometry = stripColinear(cGeometry);
        }

        double newArea = 0;
        foreach (PathD t in cGeometry)
        {
            newArea += Clipper.Area(t);
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
                        bool orientation = Clipper.IsPositive(source);
                        bool origCG0_o = Clipper.IsPositive(cGeometry[0]);
                        foreach (PathD t in cGeometry.Where(_ => origCG0_o != orientation))
                        {
                            t.Reverse();
                        }

                        break;
                    }
                }

                break;
            }
        }

        return pReorderXY(cGeometry);
    }
}