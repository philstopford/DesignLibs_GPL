using Clipper2Lib;
using System;
using System.Linq;
using System.Threading.Tasks;
using System.Collections.Generic;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public const double decomp_keyhole_sizing = 50;
    public enum outerCutterIndex { outer, cutter }
    
    public static PathsD[] getOutersAndCutters(PathsD source)
    {
        return pGetOutersAndCutters(source);
    }

    private static PathsD[] pGetOutersAndCutters(PathsD source)
    {
        PathsD[] ret = new PathsD[2];
        // Find cutters and outers.
        PathsD outers = [];
        PathsD cutters = [];
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

    public enum Type { outer, cutter }

    public static PathsD decompose(PathsD source)
    {
        return pDecompose(source);
    }

    private static PathsD pDecompose(PathsD source)
    {
        switch (source.Count)
        {
            case < 1:
                return new PathsD(source);
        }

        PathsD ret = [];

        ClipperD c = new(Constants.roundingDecimalPrecision);

        // Reconcile each path separately to get a clean representation.
        foreach (PathD t1 in source)
        {
            double a1 = Clipper.Area(t1);
            c.Clear();
            c.AddSubject(t1);
            PathsD t = [];
            c.Execute(ClipType.Union, FillRule.EvenOdd, t);
            t = pReorderXY(t);
            double a2 = t.Sum(Clipper.Area);
            
            switch (Math.Abs(Math.Abs(a1) - Math.Abs(a2)))
            {
                case <= Constants.tolerance:
                    // shape didn't really change.
                    ret.Add(t1);
                    break;
                default:
                {
                    // Orientation tracking.
                    bool origOrient = Clipper.IsPositive(t1);

                    c.AddSubject(source);

                    PathsD cR = [];
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

                c.AddSubject(decomp[(int)Type.outer]);
                c.AddClip(decomp[(int)Type.cutter]);

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
        ret[0] = [];
        ret[1] = [];
        
        // First path in source is always the outer orientation.
        bool outerOrient = Clipper.IsPositive(source[0]);
        
        foreach (PathD t in source)
        {
            // Outer was wrongly oriented, so fix up the current path for consistency.
            if (!outerOrient)
            {
                t.Reverse();
            }

            int r = (int)Type.outer;
            if (!Clipper.IsPositive(t)) // cutter
            {
                r = (int)Type.cutter;
            }

            ret[r].Add(new PathD(t));
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
        PathsD ret = [];

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
        PathsD ret = [new PathD(_poly)];

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

                PathsD decomp = decompose_poly_to_rectangles(ref abort, new PathD(ret[i]), maxRayLength, angularTolerance, vertical);
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

        PathD lPoly = pClose(_poly);

        switch (_poly.Count)
        {
            case 5 when orthogonal(stripTerminators(_poly, false), angularTolerance):
                return [_poly];
        }

        // dirOverride switches from a horizontally-biased raycast to vertical in this case.
        RayCast rc = new(lPoly, lPoly, maxRayLength, projectCorners: true, invert: RayCast.inversionMode.x, runOuterLoopThreaded:true, runInnerLoopThreaded: true, dirOverride: vertical ? RayCast.forceSingleDirection.vertical : RayCast.forceSingleDirection.horizontal);
        //***************************************************************************
        // Get rays (unchanged)
        PathsD rays = rc.getRays();

        // Contains edges from ray intersections that are not part of the original geometry.
        PathsD newEdges = [];

        ClipperD c = new(Constants.roundingDecimalPrecision);

        // ---------- BATCHED INTERSECTION + VERTEX LOOKUP ----------
        // If no rays, fallback to original behavior (nothing to do).
        if (rays == null || rays.Count == 0)
        {
            // nothing found
        }
        else
        {
            // Local function to create a stable key for a vertex using the tolerance magnitude.
            static string KeyFor(PointD p, double tolerance)
            {
                // Number of decimals to round to derived from tolerance (safety: clamp between 0 and 10).
                int digits = 6;
                if (tolerance > 0)
                {
                    double log = -Math.Floor(Math.Log10(tolerance));
                    digits = (int)Math.Max(0, Math.Min(10, log));
                }
                // Use invariant formatting to ensure stable keys.
                return $"{Math.Round(p.x, digits):F{digits}}:{Math.Round(p.y, digits):F{digits}}";
            }

            // Build vertex lookup for O(1) endpoint membership checks.
            var vertexLookup = new HashSet<string>();
            foreach (PointD v in lPoly)
            {
                vertexLookup.Add(KeyFor(v, Constants.tolerance));
            }

            // Batch all rays into a single Clipper execute - add all open subjects then the polygon as clip.
            c.Clear();
            foreach (PathD r in rays)
            {
                if (abort) break;
                c.AddOpenSubject(r);
            }
            c.AddClip(lPoly);

            PolyTreeD pt = new();
            PathsD p = new PathsD();

            c.Execute(ClipType.Intersection, FillRule.EvenOdd, pt, p);
            c.Clear();

            // Filter segments: only keep those that have at least one endpoint on an existing polygon vertex.
            for (int idx = p.Count - 1; idx >= 0; idx--)
            {
                if (abort) break;
                var seg = p[idx];
                if (seg.Count < 2)
                {
                    p.RemoveAt(idx);
                    continue;
                }

                var aKey = KeyFor(seg[0], Constants.tolerance);
                var bKey = KeyFor(seg[1], Constants.tolerance);

                if (!vertexLookup.Contains(aKey) && !vertexLookup.Contains(bKey))
                {
                    p.RemoveAt(idx);
                    continue;
                }
            }

            // Axis filter (vertical/horizontal) -- same logic as original but performed after the batch.
            for (int idx = p.Count - 1; idx >= 0; idx--)
            {
                if (abort) break;
                var seg = p[idx];
                if (vertical)
                {
                    if (Math.Abs(seg[0].x - seg[1].x) > Constants.tolerance)
                    {
                        p.RemoveAt(idx);
                    }
                }
                else
                {
                    if (Math.Abs(seg[0].y - seg[1].y) > Constants.tolerance)
                    {
                        p.RemoveAt(idx);
                    }
                }
            }

            // Use the remaining candidate segments for the original offset-area screening.
            if (p.Count > 0 && !abort)
            {
                bool breakOut = false;
                foreach (PathD t1 in p)
                {
                    if (abort) break;

                    bool edgeIsNew = true;
                    // Use an area check to see if our edge was somehow coincident with the original geometry.
                    ClipperOffset co = new();
                    co.AddPath(Clipper.ScalePath64(t1, 1.0), JoinType.Square, EndType.Butt);
                    Paths64 inflated = [];
                    co.Execute(1.0, inflated);

                    double orig_area = Clipper.Area(inflated);
                    Paths64 intersect = Clipper.Intersect(inflated, [Clipper.ScalePath64(lPoly, 1.0)], FillRule.EvenOdd);
                    double intersect_area = Clipper.Area(intersect);

                    // If we have a coincident edge, half of the offset will be inside the polygon and half outside, so we should get a measurable area difference.
                    if (Math.Abs((Math.Abs(orig_area) * 0.5) - Math.Abs(intersect_area)) < 0.0001 )
                    {
                        edgeIsNew = false;
                    }

                    if (!edgeIsNew)
                    {
                        continue;
                    }

                    newEdges.Add(t1);
                    breakOut = true;
                    break;
                }
                // NOTE: original code broke out of the outer ray loop once one new edge was added; we maintain that behavior by stopping here.
            }
        }
        // ---------- end batched section ----------
        //***************************************************************************

        PathsD final = [];
        switch (newEdges.Count)
        {
            case > 0 when !abort:
            {
                Paths64 f = [];
                Paths64 rescaledSources = _pPaths64FromPathsD(newEdges, Constants.scalar_1E2);

                // We tried to screen candidate edges earlier, but this is a defensive approach to avoid decomposition failures.
                // We may have some falsely detected 'new' edges, so we iterate our candidates to try and find one that increases the polygon count.
                foreach (Path64 rescaledSource in rescaledSources)
                {
                    f.Clear();
                    // Turn the new edges into cutters and slice. Not terribly elegant and we're relying on rounding to squash notches later.
                    // Floating points cause trouble here - we need to snap the edges to integer intervals to avoid creating internal edges.
                    ClipperOffset co = new() { PreserveCollinear = true };
                    co.AddPath(rescaledSource, JoinType.Miter, EndType.Square);

                    // Width is 2 for 1 unit each side (+/-), and the second value below is to balance the cut.
                    Paths64 cutters = [];
                    co.Execute(1.0, cutters);

                    Clipper64 c1 = new();
                    c1.AddSubject(_pPath64FromPathD(lPoly, Constants.scalar_1E2));

                    // Take first cutter only - we only cut once, no matter how many potential cutters we have.
                    c1.AddClip(cutters[0]);
                    c1.Execute(ClipType.Difference, FillRule.EvenOdd, f);

                    // Squash our notches by scaling the integers back down. The notches disappear when they can't
                    // be represented by integer values.
                    f = Clipper.ScalePaths(f, Constants.scalar_1E2_inv);

                    // Did we actually get a bisection? If not, try a different candidate edge.
                    if (f.Count > 1)
                    {
                        break;
                    }
                    // Ideally we don't get here. Use this to figure out if we do...
                    int x = 2;
                }
                
                // Clean-up.
                f = simplify(f);
                
                // Conversion.
                final = _pPathsDFromPaths64(f, 1);
                
                break;
            }
        }

        return new PathsD(final);
    }
}