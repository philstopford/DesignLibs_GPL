using Clipper2Lib;
using System;
using System.Linq;
using System.Threading.Tasks;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public enum outerCutterIndex { outer, cutter }
    
    public static PathsD[] getOutersAndCutters(PathsD source)
    {
        return pGetOutersAndCutters(source);
    }

    private static PathsD[] pGetOutersAndCutters(PathsD source)
    {
        PathsD[] ret = new PathsD[2];
        // Find cutters and outers.
        PathsD outers = new();
        PathsD cutters = new();
        foreach (PathD t in source)
        {
            if (pIsClockwise(t))
            {
                outers.Add(t);
            }
            else
            {
                cutters.Add(t);
            }
        }
        ret[(int)outerCutterIndex.outer] = outers;
        ret[(int)outerCutterIndex.cutter] = cutters;

        return ret;
    }

    public enum type { outer, cutter }

    public static PathsD decompose(PathsD source)
    {
        return pDecompose(source);
    }

    private static PathsD pDecompose(PathsD source)
    {
        switch (source.Count)
        {
            case < 1:
                return new(source);
        }

        PathsD ret = new();

        ClipperD c = new(constants.roundingDecimalPrecision);

        // Reconcile each path separately to get a clean representation.
        foreach (PathD t1 in source)
        {
            double a1 = Clipper.Area(t1);
            c.Clear();
            c.AddSubject(t1);
            PathsD t = new();
            c.Execute(ClipType.Union, FillRule.EvenOdd, t);
            t = pReorderXY(t);
            double a2 = t.Sum(Clipper.Area);
            
            switch (Math.Abs(Math.Abs(a1) - Math.Abs(a2)))
            {
                case <= constants.tolerance:
                    // shape didn't really change.
                    ret.Add(t1);
                    break;
                default:
                {
                    // Orientation tracking.
                    bool origOrient = Clipper.IsPositive(t1);

                    c.AddSubject(source);

                    PathsD cR = new();
                    // Non-zero here means that we also reconcile self-intersections without odd-even causing holes; positive only respects a certain orientation (unlike non-zero)
                    // Union is cheaper than finding the bounding box and using intersection; test-bed showed identical results.
                    c.Execute(ClipType.Union, FillRule.NonZero, cR);
                    cR = pReorderXY(cR);

                    int crCount = cR.Count;

                    // Review orientation. Fix if needed.
                    if (Clipper.IsPositive(cR[0]) != origOrient)
                    {
#if !GWSINGLETHREADED
                        Parallel.For(0, crCount, j =>
#else
                        for (int j = 0; j < crCount; j++)
#endif
                            {
                                cR[j].Reverse();
                            }
#if !GWSINGLETHREADED
                        );
#endif
                    }

                    ret.AddRange(cR);
                    c.Clear();
                    break;
                }
            }
        }
        switch (ret.Count)
        {
            case > 1:
            {
                // Need to reverse the orientations if Clipper indicates false here.
                bool reverse = !Clipper.IsPositive(ret[0]);

                switch (reverse)
                {
                    case true:
                    {
                        int rCount = ret.Count;
#if !GWSINGLETHREADED
                        Parallel.For(0, rCount, i =>
#else
                        for (int i = 0; i < rCount; i++)
#endif
                        {
                            ret[i].Reverse();
                        }
#if !GWSINGLETHREADED
                        );
#endif
                        break;
                    }
                }

                PathsD[] decomp = pGetDecomposed(ret);

                c.AddSubject(decomp[(int)type.outer]);
                c.AddClip(decomp[(int)type.cutter]);

                c.Execute(ClipType.Difference, FillRule.EvenOdd, ret);//, PolyFillType.pftPositive, PolyFillType.pftNegative);

                ret = pReorderXY(ret);
                
                switch (ret.Count)
                {
                    // Assume tone is wrong. We should not trigger this with the 'reverse' handling above.
                    case 0:
                        ret.Clear();
                        c.Execute(ClipType.Difference, FillRule.EvenOdd, ret);
                        ret = pReorderXY(ret);
                        break;
                }

                break;
            }
        }

        ret = pClose(ret);

        return ret;
    }

    public static PathsD[] getDecomposed(PathsD source)
    {
        return pGetDecomposed(source);
    }

    private static PathsD[] pGetDecomposed(PathsD source_)
    {
        PathsD source = new(source_);
        PathsD[] ret = new PathsD[2];
        ret[0] = new ();
        ret[1] = new ();
        
        // First path in source is always the outer orientation.
        bool outerOrient = Clipper.IsPositive(source[0]);
        
        foreach (PathD t in source)
        {
            // Outer was wrongly oriented, so fix up the current path for consistency.
            if (!outerOrient)
            {
                t.Reverse();
            }

            int r = (int)type.outer;
            if (!Clipper.IsPositive(t)) // cutter
            {
                r = (int)type.cutter;
            }

            ret[r].Add(new (t));
        }

        return ret;
    }


    // Scaling value below is because the incoming geometry is upsized to allow minor notches to be discarded in the conversion back to ints. Default value provided based on testing.
    public static PathsD rectangular_decomposition(ref bool abort, PathsD polys, long maxRayLength=-1, double angularTolerance = 0, bool vertical= true)
    {
        return pRectangular_decomposition(ref abort, polys, maxRayLength, angularTolerance, vertical);
    }

    private static PathsD pRectangular_decomposition(ref bool abort, PathsD polys_, long maxRayLength, double angularTolerance, bool vertical)
    {
        PathsD polys = new(polys_);
        PathsD ret = new();

        foreach (PathD t in polys)
        {
            if (abort)
            {
                ret.Clear();
                break;
            }
            ret.AddRange(pRectangular_decomposition(ref abort, t, maxRayLength, angularTolerance, vertical));
        }
        
        return ret;
    }
    public static PathsD rectangular_decomposition(ref bool abort, PathD _poly, long maxRayLength=-1, double angularTolerance = 0, bool vertical = true)
    {
        return pRectangular_decomposition(ref abort, _poly, maxRayLength, angularTolerance, vertical);
    }

    private static PathsD pRectangular_decomposition(ref bool abort, PathD _poly, long maxRayLength, double angularTolerance, bool vertical)
    {
        PathsD ret = new() {new(_poly)};

        bool changed = true;
        int startIndex = 0;
        while (changed)
        {
            if (abort)
            {
                ret.Clear();
                break;
            }
            // Set so that we break out of the loop if nothing changes (i.e. decomposition is no longer possible).
            changed = false;

            int retCount = ret.Count;
            for (int i = startIndex; i < retCount; i++)
            {
                if (abort)
                {
                    ret.Clear();
                    break;
                }
                PathsD decomp = decompose_poly_to_rectangles(ref abort, new (ret[i]), maxRayLength, angularTolerance, vertical);
                // If we got more than one polygon back, we decomposed across an internal edge.
                if (decomp.Count <= 1)
                {
                    continue;
                }

                // We decomposed something and need to store the new elements.
                // Remove original polygon from the list - we only want the decomposed entries in future.
                ret.RemoveAt(i);
                // Add our decomposed geometry to the end of the list for re-scan.
                ret.AddRange(decomp);
                // Set our start for the next pass to be this entry in the list. We avoid re-checking known-good polygons in this way.
                startIndex = i;
                // Flag that we changed something.
                changed = true;
                // Break out of this loop to start again.
                break;
            }
        }

        return ret;
    }

    private static PathsD decompose_poly_to_rectangles(ref bool abort, PathD poly, long maxRayLength, double angularTolerance, bool vertical)
    {
        PathD _poly = pClockwiseAndReorderXY(new PathD(poly));
        // Path64 lPoly = pathFromPoint(_poly, scaling);

        PathD lPoly = pClose(_poly);

        switch (_poly.Count)
        {
            case 5 when orthogonal(stripTerminators(_poly, false), angularTolerance):
                return new () { _poly };
        }

        // dirOverride switches from a horizontally-biased raycast to vertical in this case.
        RayCast rc = new(lPoly, lPoly, maxRayLength, projectCorners: true, invert: RayCast.inversionMode.x, runOuterLoopThreaded:true, runInnerLoopThreaded: true, dirOverride: vertical ? RayCast.forceSingleDirection.vertical : RayCast.forceSingleDirection.horizontal);

        PathsD rays = rc.getRays();

        // Contains edges from ray intersections that are not part of the original geometry.
        PathsD newEdges = new();

        ClipperD c = new(constants.roundingDecimalPrecision);

        foreach (PathD t in rays)
        {
            if (abort)
            {
                break;
            }
            c.AddOpenSubject(t);
            c.AddClip(lPoly);

            PolyTreeD pt = new();
            PathsD p = new();

            c.Execute(ClipType.Intersection, FillRule.EvenOdd, pt, p);
            c.Clear();
            
            int pCount = p.Count;

            switch (pCount)
            {
                case > 0:
                {
                    for (int path = pCount - 1; path >= 0; path--)
                    {
                        if (abort)
                        {
                            break;
                        }
                        double aDist = double.MaxValue;
                        double bDist = double.MaxValue;
                        // See whether the start or end point exists in the lPoly geometry. If not, we should drop this path from the list.
                        foreach (PointD t1 in lPoly)
                        {
                            if (abort)
                            {
                                break;
                            }
                            double aDist_t = distanceBetweenPoints(t1, p[path][0]);
                            double bDist_t = distanceBetweenPoints(t1, p[path][1]);

                            aDist = Math.Min(aDist_t, aDist);
                            bDist = Math.Min(bDist_t, bDist);
                        }

                        double minDist = Math.Min(aDist, bDist);

                        switch (minDist)
                        {
                            // Remove any path which has a min distance to existing points of more than 0.
                            //|| ((aDist == 0) && (bDist == 0)))
                            case > 0:
                                p.RemoveAt(path);
                                break;
                        }
                    }

                    break;
                }
            }

            if (p.Count <= 0)
            {
                continue;
            }

            {
                int pCount_ = p.Count;
                for (int p_ = pCount_ - 1; p_ >= 0; p_--)
                {
                    if (abort)
                    {
                        break;
                    }

                    switch (vertical)
                    {
                        case true:
                        {
                            if (Math.Abs(p[p_][0].x - p[p_][1].x) > constants.tolerance)
                            {
                                p.RemoveAt(p_);
                            }

                            break;
                        }
                        default:
                        {
                            if (Math.Abs(p[p_][0].y - p[p_][1].y) > constants.tolerance)
                            {
                                p.RemoveAt(p_);
                            }

                            break;
                        }
                    }
                }

                switch (p.Count)
                {
                    case 0:
                        continue;
                }

                // Should only have at least one path in the result, hopefully with desired direction. Could still have more than one, though.

                bool breakOut = false;
                foreach (PathD t1 in p)
                {
                    if (abort)
                    {
                        break;
                    }
                    bool edgeIsNew = true;
                    for (int e = 0; e < lPoly.Count - 1; e++)
                    {
                        if (abort)
                        {
                            break;
                        }
                        if (Math.Abs(lPoly[e].x - t1[0].x) < constants.tolerance && Math.Abs(lPoly[e].y - t1[0].y) < constants.tolerance)
                        {
                            int nextIndex = (e + 1) % lPoly.Count;
                            if (Math.Abs(lPoly[nextIndex].x - t1[1].x) < constants.tolerance && Math.Abs(lPoly[nextIndex].y - t1[1].y) < constants.tolerance)
                            {
                                edgeIsNew = false;
                            }
                        }

                        switch (edgeIsNew)
                        {
                            case true:
                            {
                                if (Math.Abs(lPoly[e].x - t1[1].x) < constants.tolerance && Math.Abs(lPoly[e].y - t1[1].y) < constants.tolerance)
                                {
                                    int nextIndex = (e + 1) % lPoly.Count;
                                    if (Math.Abs(lPoly[nextIndex].x - t1[0].x) < constants.tolerance && Math.Abs(lPoly[nextIndex].y - t1[0].y) < constants.tolerance)
                                    {
                                        edgeIsNew = false;
                                    }
                                }

                                break;
                            }
                        }
                    }

                    if (!edgeIsNew)
                    {
                        continue;
                    }

                    newEdges.Add(t1);// new Path(p[0]));
                    breakOut = true;
                    break;
                }
                if (breakOut)
                {
                    break;
                }
            }
        }

        PathsD final = new();
        switch (newEdges.Count)
        {
            case > 0 when !abort:
            {
                // Turn the new edges into cutters and slice. Not terribly elegant and we're relying on rounding to squash notches later.
                // Floating points cause trouble here - we need to snap the edges to integer intervals to avoid creating internal edges.
                Paths64 rescaledSource = _pPaths64FromPathsD(newEdges, constants.scalar_1E3);
                
                ClipperOffset co = new() {PreserveCollinear = true};
                co.AddPaths(rescaledSource, JoinType.Miter, EndType.Square);

                // Width is 2 for 1 unit each side (+/-), and the second value below is to balance the cut.
                Paths64 cutters = co.Execute(2.0);
                
                Clipper64 c1 = new();
                c1.AddSubject(_pPath64FromPathD(lPoly, constants.scalar_1E3));

                // Take first cutter only - we only cut once, no matter how many potential cutters we have.
                c1.AddClip(cutters[0]);
                Paths64 f = new();
                c1.Execute(ClipType.Difference, FillRule.EvenOdd, f);
                
                // Squash our notches by scaling the integers back down. The notches disappear when they can't
                // be represented by integer values.
                f = Clipper.ScalePaths(f, constants.scalar_1E3_inv);
                
                // Clean-up.
                f = simplify(f);
                
                // Conversion.
                final = _pPathsDFromPaths64(f, 1);
                
                break;
            }
        }

        return new(final);
    }
}