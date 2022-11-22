using Clipper2Lib;
using geoLib;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public enum type { outer, cutter }

    public static Paths64 decompose(Paths64 source)
    {
        return pDecompose(source);
    }

    private static Paths64 pDecompose(Paths64 source)
    {
        switch (source.Count)
        {
            case < 1:
                return source;
        }

        Paths64 ret = new();

        Clipper64 c = new();

        // Reconcile each path separately to get a clean representation.
        foreach (Path64 t1 in source)
        {
            double a1 = Clipper.Area(t1);
            c.Clear();
            c.AddSubject(t1);
            Paths64 t = new();
            c.Execute(ClipType.Union, FillRule.EvenOdd, t);
            t = pReorderXY(t);
            double a2 = t.Sum(Clipper.Area);

            switch (Math.Abs(Math.Abs(a1) - Math.Abs(a2)))
            {
                case <= double.Epsilon:
                    // shape didn't really change.
                    ret.Add(t1);
                    break;
                default:
                {
                    // Orientation tracking.
                    bool origOrient = Clipper.IsPositive(t1);

                    c.AddSubject(source);

                    Paths64 cR = new();
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

                Paths64[] decomp = pGetDecomposed(ret);

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

    public static Paths64[] getDecomposed(Paths64 source)
    {
        return pGetDecomposed(source);
    }

    private static Paths64[] pGetDecomposed(Paths64 source)
    {
        Paths64[] ret = new Paths64[2];
        ret[0] = new ();
        ret[1] = new ();
        
        // First path in source is always the outer orientation.
        bool outerOrient = Clipper.IsPositive(source[0]);
        
        foreach (Path64 t in source)
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
    public static Paths64 rectangular_decomposition(ref bool abort, Paths64 polys, long maxRayLength=-1, double angularTolerance = 0, bool vertical= true)
    {
        return pRectangular_decomposition(ref abort, polys, maxRayLength, angularTolerance, vertical);
    }

    private static Paths64 pRectangular_decomposition(ref bool abort, Paths64 polys, long maxRayLength, double angularTolerance, bool vertical)
    {
        Paths64 ret = new();

        foreach (Path64 t in polys)
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
    public static Paths64 rectangular_decomposition(ref bool abort, Path64 _poly, long maxRayLength=-1, double angularTolerance = 0, bool vertical = true)
    {
        return pRectangular_decomposition(ref abort, _poly, maxRayLength, angularTolerance, vertical);
    }

    private static Paths64 pRectangular_decomposition(ref bool abort, Path64 _poly, long maxRayLength, double angularTolerance, bool vertical)
    {
        Paths64 ret = new() {new(_poly)};

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
                Paths64 decomp = decompose_poly_to_rectangles(ref abort, new (ret[i]), maxRayLength, angularTolerance, vertical);
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

    private static Paths64 decompose_poly_to_rectangles(ref bool abort, Path64 _poly, long maxRayLength, double angularTolerance, bool vertical)
    {
        _poly = pClockwiseAndReorderXY(_poly);
        // Path64 lPoly = pathFromPoint(_poly, scaling);

        Path64 lPoly = pClose(_poly);

        switch (_poly.Count)
        {
            case 5 when orthogonal(stripTerminators(_poly, false), angularTolerance):
                return new () { _poly };
        }

        // dirOverride switches from a horizontally-biased raycast to vertical in this case.
        RayCast rc = new(lPoly, lPoly, maxRayLength, projectCorners: true, invert: RayCast.inversionMode.x, runOuterLoopThreaded:true, runInnerLoopThreaded: true, dirOverride: vertical ? RayCast.forceSingleDirection.vertical : RayCast.forceSingleDirection.horizontal);

        Paths64 rays = rc.getRays();

        // Contains edges from ray intersections that are not part of the original geometry.
        Paths64 newEdges = new();

        Clipper64 c = new();

        foreach (Path64 t in rays)
        {
            if (abort)
            {
                break;
            }
            c.AddOpenSubject(t);
            c.AddClip(lPoly);

            PolyTree64 pt = new();
            Paths64 p = new();

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
                        foreach (Point64 t1 in lPoly)
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
                            if (p[p_][0].X != p[p_][1].X)
                            {
                                p.RemoveAt(p_);
                            }

                            break;
                        }
                        default:
                        {
                            if (p[p_][0].Y != p[p_][1].Y)
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
                foreach (Path64 t1 in p)
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
                        if (lPoly[e].X == t1[0].X && lPoly[e].Y == t1[0].Y)
                        {
                            int nextIndex = (e + 1) % lPoly.Count;
                            if (lPoly[nextIndex].X == t1[1].X && lPoly[nextIndex].Y == t1[1].Y)
                            {
                                edgeIsNew = false;
                            }
                        }

                        switch (edgeIsNew)
                        {
                            case true:
                            {
                                if (lPoly[e].X == t1[1].X && lPoly[e].Y == t1[1].Y)
                                {
                                    int nextIndex = (e + 1) % lPoly.Count;
                                    if (lPoly[nextIndex].X == t1[0].X && lPoly[nextIndex].Y == t1[0].Y)
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

        Paths64 final = new();
        switch (newEdges.Count)
        {
            case > 0 when !abort:
            {
                // Turn the new edges into cutters and slice. Not terribly elegant and we're relying on rounding to squash notches later.
                ClipperOffset co = new() {PreserveCollinear = true};
                co.AddPaths(newEdges, JoinType.Miter, EndType.Square);

                Paths64 cutters = co.Execute(2.0);

                c.Clear();

                c.AddSubject(lPoly);

                // Take first cutter only - we only cut once, no matter how many potential cutters we have.
                c.AddClip(cutters[0]);
                Paths64 f = new();
                c.Execute(ClipType.Difference, FillRule.EvenOdd, f);

                f = pReorderXY(f);
                
                final = pClose(f);

                final = simplify(final);

                final = clockwiseAndReorderXY(final);
                break;
            }
        }

        return new(final);
    }
}