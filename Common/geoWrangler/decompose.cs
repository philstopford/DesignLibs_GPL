using ClipperLib;
using geoLib;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace geoWrangler;

using Path = List<IntPoint>;
using Paths = List<List<IntPoint>>;

public static partial class GeoWrangler
{
    public enum type { outer, cutter }

    public static Paths decompose(Paths source)
    {
        return pDecompose(source);
    }

    private static Paths pDecompose(Paths source)
    {
        switch (source.Count)
        {
            case < 1:
                return source;
        }

        Paths ret = new();

        Clipper c = new() {PreserveCollinear = true};

        // Reconcile each path separately to get a clean representation.
        foreach (Path t1 in source)
        {
            double a1 = Clipper.Area(t1);
            double a2 = 0;
            c.Clear();
            c.AddPath(t1, PolyType.ptSubject, true);
            Paths t = new();
            c.Execute(ClipType.ctUnion, t);
            foreach (Path t2 in t)
            {
                a2 += Clipper.Area(t2);
            }

            switch (Math.Abs(Math.Abs(a1) - Math.Abs(a2)))
            {
                case <= double.Epsilon:
                    // shape didn't really change.
                    ret.Add(t1);
                    break;
                default:
                {
                    // Orientation tracking.
                    bool origOrient = Clipper.Orientation(t1);

                    c.AddPaths(source, PolyType.ptSubject, true);

                    Paths cR = new();
                    // Non-zero here means that we also reconcile self-intersections without odd-even causing holes; positive only respects a certain orientation (unlike non-zero)
                    // Union is cheaper than finding the bounding box and using intersection; test-bed showed identical results.
                    c.Execute(ClipType.ctUnion, cR, PolyFillType.pftNonZero);

                    int crCount = cR.Count;

                    // Review orientation. Fix if needed.
                    if (Clipper.Orientation(cR[0]) != origOrient)
                    {
#if !GWSINGLETHREADED
                        Parallel.For(0, crCount, (j) =>
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
                bool reverse = !Clipper.Orientation(ret[0]);

                switch (reverse)
                {
                    case true:
                    {
                        int rCount = ret.Count;
#if !GWSINGLETHREADED
                        Parallel.For(0, rCount, (i) =>
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

                Paths[] decomp = pGetDecomposed(ret);

                c.AddPaths(decomp[(int)type.outer], PolyType.ptSubject, true);
                c.AddPaths(decomp[(int)type.cutter], PolyType.ptClip, true);

                c.Execute(ClipType.ctDifference, ret, PolyFillType.pftPositive, PolyFillType.pftNegative);

                switch (ret.Count)
                {
                    // Assume tone is wrong. We should not trigger this with the 'reverse' handling above.
                    case 0:
                        ret.Clear();
                        c.Execute(ClipType.ctDifference, ret);
                        break;
                }

                break;
            }
        }

        ret = pClose(ret);

        return ret;
    }

    public static Paths[] getDecomposed(Paths source)
    {
        return pGetDecomposed(source);
    }

    private static Paths[] pGetDecomposed(Paths source)
    {
        Paths[] ret = new Paths[2];
        ret[0] = new Paths();
        ret[1] = new Paths();

        foreach (Path t in source)
        {
            int r = (int)type.outer;
            if (!Clipper.Orientation(t))
            {
                r = (int)type.cutter;
            }

            ret[r].Add(new Path(t));
        }

        return ret;
    }


    // Scaling value below is because the incoming geometry is upsized to allow minor notches to be discarded in the conversion back to ints. Default value provided based on testing.
    public static List<GeoLibPoint[]> rectangular_decomposition(ref bool abort, List<GeoLibPoint[]> polys, int scaling = 10000, long maxRayLength=-1, double angularTolerance = 0, bool vertical= true)
    {
        return pRectangular_decomposition(ref abort, polys, scaling, maxRayLength, angularTolerance, vertical);
    }

    private static List<GeoLibPoint[]> pRectangular_decomposition(ref bool abort, List<GeoLibPoint[]> polys, int scaling, long maxRayLength, double angularTolerance, bool vertical)
    {
        List<GeoLibPoint[]> ret = new();

        foreach (GeoLibPoint[] t in polys)
        {
            if (abort)
            {
                ret.Clear();
                break;
            }
            ret.AddRange(pRectangular_decomposition(ref abort, t, scaling, maxRayLength, angularTolerance, vertical));
        }

        return ret;
    }
    public static List<GeoLibPoint[]> rectangular_decomposition(ref bool abort, GeoLibPoint[] _poly, int scaling = 10000, long maxRayLength=-1, double angularTolerance = 0, bool vertical = true)
    {
        return pRectangular_decomposition(ref abort, _poly, scaling, maxRayLength, angularTolerance, vertical);
    }

    private static List<GeoLibPoint[]> pRectangular_decomposition(ref bool abort, GeoLibPoint[] _poly, int scaling, long maxRayLength, double angularTolerance, bool vertical)
    {
        List<GeoLibPoint[]> ret = new() {_poly.ToArray()};

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
                List<GeoLibPoint[]> decomp = decompose_poly_to_rectangles(ref abort, ret[i].ToArray(), scaling, maxRayLength, angularTolerance, vertical);
                // If we got more than one polygon back, we decomposed across an internal edge.
                if (decomp.Count > 1)
                {
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
        }



        return ret;
    }

    private static List<GeoLibPoint[]> decompose_poly_to_rectangles(ref bool abort, GeoLibPoint[] _poly, int scaling, long maxRayLength, double angularTolerance, bool vertical)
    {
        Path lPoly = pathFromPoint(pClockwiseAndReorder(_poly), scaling);

        switch (_poly.Length)
        {
            case 5 when orthogonal(stripTerminators(_poly, false), angularTolerance):
                return new List<GeoLibPoint[]> { _poly };
        }

        // dirOverride switches from a horizontally-biased raycast to vertical in this case.
        RayCast rc = new(lPoly, lPoly, maxRayLength * scaling, projectCorners: true, invert: true, runOuterLoopThreaded:true, runInnerLoopThreaded: true, dirOverride: vertical ? RayCast.forceSingleDirection.vertical : RayCast.forceSingleDirection.horizontal);

        Paths rays = rc.getRays();

        // Contains edges from ray intersections that are not part of the original geometry.
        Paths newEdges = new();

        Clipper c = new();

        foreach (Path t in rays)
        {
            if (abort)
            {
                break;
            }
            c.AddPath(t, PolyType.ptSubject, false);
            c.AddPath(lPoly, PolyType.ptClip, true);

            PolyTree pt = new();

            c.Execute(ClipType.ctIntersection, pt);
            c.Clear();

            Paths p = Clipper.OpenPathsFromPolyTree(pt);

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
                        foreach (IntPoint t1 in lPoly)
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

            if (p.Count > 0)
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
                foreach (Path t1 in p)
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

                    if (edgeIsNew)
                    {
                        newEdges.Add(t1);// new Path(p[0]));
                        breakOut = true;
                        break;
                    }
                }
                if (breakOut)
                {
                    break;
                }
            }
        }

        List<GeoLibPoint[]> final = new();
        switch (newEdges.Count)
        {
            case > 0 when !abort:
            {
                // Turn the new edges into cutters and slice. Not terribly elegant and we're relying on rounding to squash notches later.
                ClipperOffset co = new();
                co.AddPaths(newEdges, JoinType.jtMiter, EndType.etOpenSquare);
                PolyTree tp = new();
                co.Execute(ref tp, 1.0);

                Paths cutters = Clipper.ClosedPathsFromPolyTree(tp);

                c.Clear();

                c.AddPath(lPoly, PolyType.ptSubject, true);

                // Take first cutter only - we only cut once, no matter how many potential cutters we have.
                c.AddPath(cutters[0], PolyType.ptClip, true);
                Paths f = new();
                c.Execute(ClipType.ctDifference, f, PolyFillType.pftEvenOdd, PolyFillType.pftEvenOdd);

                final = pointsFromPaths(f, scaling);

                final = simplify(final);

                final = clockwiseAndReorder(final);
                break;
            }
        }

        return final;
    }
}