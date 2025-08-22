using Clipper2Lib;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;

namespace geoWrangler;

public static partial class GeoWrangler
{
    // Sizing is used to define the keyhole width (default) and will be used for the sliver/gap removal.
    // Use of a custom value will cause headaches.
    public const double keyhole_sizing = 2.5;
    private const double keyhole_extension_default = 2.06;

    public static PathsD makeKeyHole(PathsD outers, PathsD cutters, bool reverseEval, bool biDirectionalEval, RayCast.inversionMode invert = RayCast.inversionMode.x, double customSizing = 0, double extension = 0, double angularTolerance = 0)
    {
        return pMakeKeyHole(outers, cutters, reverseEval, biDirectionalEval, invert, customSizing, extension, angularTolerance);
    }

    public static PathsD makeKeyHole(PathD source, bool reverseEval, bool biDirectionalEval, RayCast.inversionMode invert = RayCast.inversionMode.x, double customSizing = 0, double extension = 0, double angularTolerance = 0)
    {
        return pMakeKeyHole([source], reverseEval, biDirectionalEval, invert, customSizing, extension, angularTolerance);
    }

    public static PathsD makeKeyHole(PathsD source, bool reverseEval, bool biDirectionalEval, RayCast.inversionMode invert = RayCast.inversionMode.x, double customSizing = 0, double extension = 0, double angularTolerance = 0)
    {
        return pMakeKeyHole(source, reverseEval, biDirectionalEval, invert, customSizing, extension, angularTolerance);
    }

    private static PathsD pMakeKeyHole(PathsD source, bool reverseEval, bool biDirectionalEval, RayCast.inversionMode invert = RayCast.inversionMode.x, double customSizing = 0, double extension = 0, double angularTolerance = 0)
    {
        // Reconcile any overlapping geometry.
        ClipperD c = new (Constants.roundingDecimalPrecision) {PreserveCollinear = true};
        
        switch (source.Count)
        {
            case < 1:
                return new PathsD(source);
        }

        customSizing = customSizing switch
        {
            0 => keyhole_sizing,
            _ => customSizing
        };

        // Limit the offset used in the removal otherwise we can cause self-intersections that
        // result in lots of artifacts and trouble.
        PathsD input = pRemoveFragments(source, 0.5 * customSizing, extension);
        input = pStripCollinear(input);

        switch (input.Count)
        {
            case < 2:
                return pClose(input);
        }

        PathsD[] decomp = pGetDecomposed(input);

        // No cutters found. Return.
        if (decomp[1].Count == 0)
        {
            return new PathsD(source);
        }

        // Cutters have no area (perhaps a collinear open path) - ignore them.
        if (Math.Abs(Clipper.Area(decomp[1])) <= Constants.tolerance)
        {
            return new PathsD(source);
        }

        PathsD[] odecomp = new PathsD[decomp.Length];
        for (int i = 0; i < decomp.Length; i++)
        {
            decomp[i] = pClose(decomp[i]);
            odecomp[i] = new PathsD(decomp[i]);
        }
        // Outer areas
        double maxOuterArea = double.MinValue;
        for (int i = 0; i < decomp[(int)Type.outer].Count; i++)
        {
            double a = Clipper.Area(decomp[(int)Type.outer][i]);
            if (a > maxOuterArea) maxOuterArea = a;
        }

        // So here, things get annoying. We can have nested donuts, which means that we have outers fully covered by cutters (from the larger donut).
        // Unless we massage things, these get killed as the cutters are applied en-masse to outers in the keyholer.
        
        // First, we'll run the keyholer as usual.
        PathsD ret = pMakeKeyHole(decomp[(int)Type.outer], decomp[(int)Type.cutter], reverseEval, biDirectionalEval, invert, customSizing, extension, angularTolerance);
        
        double origArea = 0.0;
        {
            PathsD tmp = sliverGapRemoval(source);
            for (int i = 0; i < tmp.Count; i++) origArea += Clipper.Area(tmp[i]);
            origArea = Math.Abs(origArea);
        }
        double newArea = 0.0;
        {
            PathsD tmp = sliverGapRemoval(ret);
            for (int i = 0; i < tmp.Count; i++) newArea += Clipper.Area(tmp[i]);
            newArea = Math.Abs(newArea);
        }

        // If we lost area, we probably had a cutter fully cover up one or more of our polygons.
        double lostArea = Math.Abs(origArea - newArea);
        
        if (lostArea > 1E-9) // Arbitrary due to the nature of floats
        {
            // We need to find out which cutters might have completely killed one or more outers and figure out a plan.
            for (int oIndex = 0; oIndex < odecomp[(int) Type.outer].Count; oIndex++)
            {
                PathD tOuter = odecomp[(int) Type.outer][oIndex];
                double outerArea = Clipper.Area(tOuter);
                if (outerArea > lostArea)
                {
                    if (Math.Abs(outerArea - maxOuterArea) <= Constants.tolerance)
                    {
                        continue;
                    }
                }
                PathsD tCutters = [];
                // Do any cutters cover up our outer?
                for (int cIndex = 0; cIndex < odecomp[(int) Type.cutter].Count; cIndex++)
                {
                    PathsD test = [];
                    c.Clear();
                    c.AddSubject(tOuter);
                    c.AddClip(odecomp[(int) Type.cutter][cIndex]);
                    c.Execute(ClipType.Difference, FillRule.EvenOdd, test);
                    test = pReorderXY(test);
                    double area = 0.0;
                    for (int ti = 0; ti < test.Count; ti++) area += Clipper.Area(test[ti]);

                    // If area is zero, the cutter fully covered the outer and we need to skip it.
                    if ((Math.Abs(area) > 0.0) && (Math.Abs(area) <= lostArea))
                    {
                        tCutters.Add(odecomp[(int) Type.cutter][cIndex]);
                    }

                }

                PathsD tOuters = [tOuter];

                PathsD tRet = pMakeKeyHole(tOuters, tCutters, reverseEval, biDirectionalEval, invert, customSizing, extension,
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
                return new PathsD(source); // something blew up. Send back the original geometry.
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
        PathsD clipped = [];

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
            rc = new RayCast(t, outers, 1000000, invert: invert, projectCorners: projectCorners);
            clipped = rc.getClippedRays();
        }

        if (!reverseEval || biDirectionalEval)
        {
            t.Reverse();
        }
        rc = new RayCast(t, outers, 1000000, invert: invert, projectCorners: projectCorners);
        clipped.AddRange(rc.getClippedRays());

        return clipped;
    }

    private static PathD pFindInsertionCandidate(PathsD edges, bool orthogonalInput)
    {
        PathD ret = null;
        // Need to find minimal length, ideally orthogonal.
        double minLengthSq = -1;
        bool minLength_ortho = false;
        int cutPathIndex = -1;
        double sqTol = Constants.tolerance * Constants.tolerance;

        for (int r = 0; r < edges.Count; r++)
        {
            if (orthogonalInput)
            {
                bool ray_isOrtho = Math.Abs(edges[r][0].x - edges[r][1].x) < Constants.tolerance || Math.Abs(edges[r][0].y - edges[r][1].y) < Constants.tolerance;
                if (minLength_ortho && !ray_isOrtho)
                {
                    // Don't replace an orthogonal ray with a non-orthogonal ray.
                    continue;
                }
                
                // Update our tracking with the ray ortho flag.
                minLength_ortho = ray_isOrtho;
            }

            double dx = edges[r][0].x - edges[r][1].x;
            double dy = edges[r][0].y - edges[r][1].y;
            double ray_length_sq = dx * dx + dy * dy;

            if (ray_length_sq <= sqTol) continue;

            if (r != 0 && !(minLengthSq < 0 || ray_length_sq < minLengthSq))
            {
                continue;
            }

            cutPathIndex = r;
            minLengthSq = ray_length_sq;
        }

        if (cutPathIndex != -1)
        {
            ret = new PathD(edges[cutPathIndex]);
        }

        return ret;
    }

    private static PathsD pCutKeyHole(PathsD outers, PathsD cutters, PathsD extraCutters)
    {
        ClipperD c = new(Constants.roundingDecimalPrecision);
        c.AddSubject(cutters);
        c.AddClip(extraCutters);

        PathsD mergedCutters = [];
        c.Execute(ClipType.Union, FillRule.EvenOdd, mergedCutters);

        mergedCutters = pReorderXY(mergedCutters);

        c.Clear();
        c.AddSubject(outers);
        c.AddClip(mergedCutters);

        // Reduce our geometry back to the simplest form.
        PathsD new_outers = [];
        c.Execute(ClipType.Difference, FillRule.EvenOdd, new_outers);//, PolyFillType.pftNonZero, PolyFillType.pftNegative);

        new_outers = pReorderXY(new_outers);
        
        return new_outers;
    }
    
    private static PathsD pMakeKeyHole(PathsD outers_, PathsD cutters_, bool reverseEval, bool biDirectionalEval, RayCast.inversionMode invert = RayCast.inversionMode.x, double customSizing = 0, double extension = 0, double angularTolerance = 0)
    {
        PathsD outers = new(outers_);
        PathsD cutters = new(cutters_);
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

                PathsD extraCutters = [];
                if (insertionCandidate != null)
                {
                    // Offset our cutter and assign to the clipping scenario.
                    PathsD sPaths = pInflateEdge(insertionCandidate, customSizing);
                    sPaths = Clipper.SimplifyPaths(sPaths, 0);

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
                new_outers = Clipper.SimplifyPaths(new_outers, 0);

                outers.Clear();
                for (int idx = 0; idx < new_outers.Count; idx++)
                {
                    PathD pth = new_outers[idx];
                    if (Clipper.IsPositive(pth) == outerOrient)
                    {
                        outers.Add(pth);
                    }
                }
            }
            
            return pClockwiseAndReorderXY(outers);
        }
        catch (Exception e)
        {
            throw new Exception("pMakeKeyHole error");
        }
    }

    private static PathsD pInflateEdge(PathD edge_, double customSizing)
    {
        customSizing = customSizing switch
        {
            0 => keyhole_sizing,
            _ => customSizing
        };
        // Force clockwise, which should get us something consistent to work with.
        PathD edge = pClockwise(edge_);

        // edge = pExtendEdge(edge, keyhole_sizing);

        // Need to workaround missing PathD support in ClipperOffset...
        // Note that the scalar has interaction with pRemoveFragments
        Path64 rescaledSource = _pPath64FromPathD(edge, Constants.scalar_1E2);

        ClipperOffset co = new() {PreserveCollinear = true};
        co.AddPath(rescaledSource, JoinType.Miter, EndType.Square);

        Paths64 tmp = [];
        co.Execute(2 * customSizing, tmp);

        // Size back down again.
        PathsD sPaths = _pPathsDFromPaths64(tmp, Constants.scalar_1E2_inv);

        return pReorderXY(sPaths);
    }

    public static PathsD sliverGapRemoval(PathD source, double customSizing = 0, double extension = 0, bool maySimplify = false)
    {
        return pSliverGapRemoval([source], customSizing, extension, maySimplify: maySimplify);
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
                new PathsD(source),
            _ => ret
        };
    }

    public static PathsD gapRemoval(PathD source, double customSizing = 0, double extension = 0, bool maySimplify = false)
    {
        return gapRemoval([source], customSizing, extension, maySimplify);
    }

    public static PathsD gapRemoval(PathsD source, double customSizing = 0, double extension = 0, bool maySimplify = false)
    {
        switch (source.Count)
        {
            case < 1:
                return new PathsD(source);
        }

        bool orig_orient_gw = isClockwise(source[0]);
        bool orig_orient_c = Clipper.IsPositive(source[0]);

        PathsD ret = pRemoveFragments(source, customSizing, extension, maySimplify);

        switch (ret.Count)
        {
            case 0:
                return new PathsD(source); // something blew up. Send back the original geometry.
        }

        // Clean-up the geometry.
        ClipperD c = new(Constants.roundingDecimalPrecision);
        c.AddSubject(ret);
        c.Execute(ClipType.Union, FillRule.EvenOdd, ret);

        ret = pReorderXY(ret);

        // ret = stripCollinear(ret);

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
        double oArea = 0.0;
        for (int i = 0; i < source.Count; i++) oArea += Clipper.Area(source[i]);
        PathsD ret = pRemoveFragments(source, -customSizing, extension, maySimplify: maySimplify);
        double nArea = 0.0;
        for (int i = 0; i < ret.Count; i++) nArea += Clipper.Area(ret[i]);

        return (Math.Abs(oArea) - Math.Abs(nArea)) switch
        {
            < 10000 => new PathsD(source),
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

        // This is an empirical 'fudge' factor in order to close keyholes reliably.
        // Testing landed on this as a minimal value which appears to work - lower values leave
        // gaps; higher values run the risk of crushing 'real' geometrical features.
        // I'm not really sure why this is now needed - earlier code worked fine.
        // Good tests for this are ILB16, ILB17.
        customSizing *= 1.05;

        extension = extension switch
        {
            0 => keyhole_extension_default,
            _ => extension
        };

        // Used to try and avoid residual fragments; empirically derived.
        customSizing *= extension;
        
        // Need to workaround missing PathD support in ClipperOffset...
        Paths64 rescaledSource = _pPaths64FromPathsD(source, Constants.scalar_1E4);

        ClipperOffset co = new() {PreserveCollinear = true};
        co.AddPaths(rescaledSource, joinType, EndType.Polygon);
        // The scalar below must be aligned with the usage in pInflateEdge
        Paths64 tmp = [];
        co.Execute(customSizing * Constants.scalar_1E2, tmp);
        co.Clear();
        co.AddPaths(new Paths64(tmp), joinType, EndType.Polygon);
        tmp.Clear();
        // The scalar below must be aligned with the usage in pInflateEdge
        co.Execute(-(customSizing  * Constants.scalar_1E2), tmp); // Size back to original dimensions

        PathsD cGeometry = _pPathsDFromPaths64(tmp, Constants.scalar_1E4_inv);
        
        cGeometry = pReorderXY(cGeometry);
        
        if (maySimplify)
        {
            cGeometry = stripCollinear(cGeometry);
        }

        double newArea = 0.0;
        for (int i = 0; i < cGeometry.Count; i++) newArea += Clipper.Area(cGeometry[i]);

        return Math.Abs(newArea) switch
        {
            <= Constants.tolerance =>
                // We crushed our geometry, it seems.
                // This was probably not the plan, so send back the original geometry instead.
                new PathsD(source),
            _ => pReorderXY(cGeometry)
        };
    }

    private static PathsD pRemoveFragments(PathD source_, double customSizing, bool maySimplify = false, JoinType joinType = JoinType.Miter)
    {
        PathD source = new(source_);
        customSizing = customSizing switch
        {
            0 => keyhole_sizing,
            _ => customSizing
        };
        double sourceArea = Clipper.Area(source);

        // Need to workaround missing PathD support in ClipperOffset...
        Path64 rescaledSource = _pPath64FromPathD(source, Constants.scalar_1E2);

        ClipperOffset co = new() {PreserveCollinear = true};
        co.AddPath(rescaledSource, joinType, EndType.Polygon);
        // The scalar below must be aligned with the usage in pInflateEdge
        Paths64 tmp = [];
        co.Execute(customSizing * Constants.scalar_1E2, tmp);
        co.Clear();
        co.AddPaths(new Paths64(tmp), joinType, EndType.Polygon);
        tmp.Clear();
        // The scalar below must be aligned with the usage in pInflateEdge
        co.Execute(-(customSizing * Constants.scalar_1E2), tmp); // Size back to original dimensions

        PathsD cGeometry = _pPathsDFromPaths64(tmp, Constants.scalar_1E2_inv);

        cGeometry = pReorderXY(cGeometry);
        
        if (maySimplify)
        {
            cGeometry = stripCollinear(cGeometry);
        }

        double newArea = 0.0;
        for (int i = 0; i < cGeometry.Count; i++) newArea += Clipper.Area(cGeometry[i]);

        switch (Math.Abs(newArea))
        {
            case <= Constants.tolerance:
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
                        // Multi-path handling gets interesting. The first path is assumed to be the outer. Let's compare that with the original geometry. If the orientation is different, reverse the necessary paths.
                        bool orientation = Clipper.IsPositive(source);
                        bool origCG0_o = Clipper.IsPositive(cGeometry[0]);
                        for (int ti = 0; ti < cGeometry.Count; ti++)
                        {
                            if (origCG0_o != orientation)
                            {
                                cGeometry[ti].Reverse();
                            }
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