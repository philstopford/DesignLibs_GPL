using ClipperLib;
using geoLib;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace geoWrangler
{
    using Path = List<IntPoint>;
    using Paths = List<List<IntPoint>>;

    public static partial class GeoWrangler
    {
        public enum type { outer, cutter }

        public static Paths decompose(Paths source)
        {
            return pDecompose(source);
        }

        static Paths pDecompose(Paths source)
        {
            if (source.Count < 1)
            {
                return source;
            }

            Paths ret = new Paths();

            Clipper c = new Clipper();
            c.PreserveCollinear = true;

            // Reconcile each path separately to get a clean representation.
            for (int i = 0; i < source.Count; i++)
            {
                double a1 = Clipper.Area(source[i]);
                double a2 = 0;
                c.Clear();
                c.AddPath(source[i], PolyType.ptSubject, true);
                Paths t = new Paths();
                c.Execute(ClipType.ctUnion, t);
                for (int j = 0; j < t.Count; j++)
                {
                    a2 += Clipper.Area(t[j]);
                }

                if (Math.Abs(Math.Abs(a1) - Math.Abs(a2)) < double.Epsilon)
                {
                    // shape didn't really change.
                    ret.Add(source[i]);
                }
                else
                {
                    // Orientation tracking.
                    bool origOrient = Clipper.Orientation(source[i]);

                    c.AddPaths(source, PolyType.ptSubject, true);

                    Paths cR = new Paths();
                    // Non-zero here means that we also reconcile self-intersections without odd-even causing holes; positive only respects a certain orientation (unlike non-zero)
                    // Union is cheaper than finding the bounding box and using intersection; test-bed showed identical results.
                    c.Execute(ClipType.ctUnion, cR, PolyFillType.pftNonZero);

                    int crCount = cR.Count;

                    // Review orientation. Fix if needed.
                    if (Clipper.Orientation(cR[0]) != origOrient)
                    {
#if GWTHREADED
                        Parallel.For(0, crCount, (j) =>
#else
                        for (int j = 0; j < cR.Count; j++)
#endif
                        {
                            cR[j].Reverse();
                        }
#if GWTHREADED
                        );
#endif
                    }

                    ret.AddRange(cR);
                    c.Clear();
                }
            }
            if (ret.Count > 1)
            {
                // Need to reverse the orientations if Clipper indicates false here.
                bool reverse = !Clipper.Orientation(ret[0]);

                if (reverse)
                {
                    int rCount = ret.Count;
#if GWTHREADED
                    Parallel.For(0, rCount, (i) =>
#else
                    for (int i = 0; i < ret.Count; i++)
#endif
                    {
                        ret[i].Reverse();
                    }
#if GWTHREADED
                    );
#endif
                }

                Paths[] decomp = pGetDecomposed(ret);

                c.AddPaths(decomp[(int)type.outer], PolyType.ptSubject, true);
                c.AddPaths(decomp[(int)type.cutter], PolyType.ptClip, true);

                c.Execute(ClipType.ctDifference, ret, PolyFillType.pftPositive, PolyFillType.pftNegative);

                // Assume tone is wrong. We should not trigger this with the 'reverse' handling above.
                if (ret.Count == 0)
                {
                    ret.Clear();
                    c.Execute(ClipType.ctDifference, ret);
                }
            }

            ret = pClose(ret);

            return ret;
        }

        public static Paths[] getDecomposed(Paths source)
        {
            return pGetDecomposed(source);
        }

        static Paths[] pGetDecomposed(Paths source)
        {
            Paths[] ret = new Paths[2];
            ret[0] = new Paths();
            ret[1] = new Paths();

            for (int i = 0; i < source.Count; i++)
            {
                int r = (int)type.outer;
                if (!Clipper.Orientation(source[i]))
                {
                    r = (int)type.cutter;
                }

                ret[r].Add(new Path(source[i]));
            }

            return ret;
        }


        // Scaling value below is because the incoming geometry is upsized to allow minor notches to be discarded in the conversion back to ints. Default value provided based on testing.
        public static List<GeoLibPoint[]> rectangular_decomposition(List<GeoLibPoint[]> polys, Int32 scaling = 10000, Int64 maxRayLength=-1, bool vertical= true)
        {
            return pRectangular_decomposition(polys, scaling, maxRayLength, vertical);
        }
        static List<GeoLibPoint[]> pRectangular_decomposition(List<GeoLibPoint[]> polys, Int32 scaling, Int64 maxRayLength, bool vertical)
        {
            List<GeoLibPoint[]> ret = new List<GeoLibPoint[]>();

            for (int i = 0; i < polys.Count; i++)
            {
                ret.AddRange(pRectangular_decomposition(polys[i], scaling, maxRayLength, vertical));
            }

            return ret;
        }
        public static List<GeoLibPoint[]> rectangular_decomposition(GeoLibPoint[] _poly, Int32 scaling = 10000, Int64 maxRayLength=-1, bool vertical = true)
        {
            return pRectangular_decomposition(_poly, scaling, maxRayLength, vertical);
        }

        static List<GeoLibPoint[]> pRectangular_decomposition(GeoLibPoint[] _poly, Int32 scaling, Int64 maxRayLength, bool vertical)
        {
            List<GeoLibPoint[]> ret = new List<GeoLibPoint[]>();
            ret.Add(_poly.ToArray());

            bool changed = true;
            int startIndex = 0;
            while (changed)
            {
                // Set so that we break out of the loop if nothing changes (i.e. decomposition is no longer possible).
                changed = false;

                int retCount = ret.Count;
                for (int i = startIndex; i < retCount; i++)
                {
                    List<GeoLibPoint[]> decomp = decompose_poly_to_rectangles(ret[i].ToArray(), scaling, maxRayLength, vertical);
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

        static List<GeoLibPoint[]> decompose_poly_to_rectangles(GeoLibPoint[] _poly, Int32 scaling, Int64 maxRayLength, bool vertical)
        {
            Path lPoly = GeoWrangler.pathFromPoint(pClockwiseAndReorder(_poly), scaling);

            if ((_poly.Length == 5) && GeoWrangler.orthogonal(GeoWrangler.stripTerminators(_poly, false)))
            {
                return new List<GeoLibPoint[]>() { _poly };
            }

            // dirOverride switches from a horizontally-biased raycast to a vertical case in this case.
            RayCast rc = new RayCast(lPoly, lPoly, maxRayLength * scaling, projectCorners: true, invert: true, runOuterLoopThreaded:true, runInnerLoopThreaded: true, dirOverride: vertical ? RayCast.forceSingleDirection.vertical : RayCast.forceSingleDirection.horizontal);

            Paths rays = rc.getRays();

            Paths newEdges = new Paths();

            Clipper c = new Clipper();

            for (int r = 0; r < rays.Count; r++)
            {
                c.AddPath(rays[r], PolyType.ptSubject, false);
                c.AddPath(lPoly, PolyType.ptClip, true);

                PolyTree pt = new PolyTree();

                c.Execute(ClipType.ctIntersection, pt);
                c.Clear();

                Paths p = Clipper.OpenPathsFromPolyTree(pt);

                int pCount = p.Count;

                if (pCount > 0)
                {
                    for (int path = pCount - 1; path >= 0; path--)
                    {
                        double minDist = 100000;
                        // See whether the start or end point exists in the lPoly geometry. If not, we should drop this path from the list.
                        for (int lPolyPt = 0; lPolyPt < lPoly.Count; lPolyPt++)
                        {
                            double aDist = GeoWrangler.distanceBetweenPoints(lPoly[lPolyPt], p[path][0]);
                            double bDist = GeoWrangler.distanceBetweenPoints(lPoly[lPolyPt], p[path][1]);

                            minDist = Math.Min(minDist, aDist);
                            minDist = Math.Min(minDist, bDist);
                        }

                        if (minDist > 1000)
                        {
                            p.RemoveAt(path);
                        }
                    }

                }

                if (p.Count > 0)
                {
                    // Should only have one path in the result.
                    bool edgeIsNew = true;
                    for (int e = 0; e < lPoly.Count - 1; e++)
                    {
                        if ((lPoly[e].X == p[0][0].X) && (lPoly[e].Y == p[0][0].Y))
                        {
                            int nextIndex = (e + 1) % lPoly.Count;
                            if ((lPoly[nextIndex].X == p[0][1].X) && (lPoly[nextIndex].Y == p[0][1].Y))
                            {
                                edgeIsNew = false;
                            }
                        }

                        if (edgeIsNew)
                        {
                            if ((lPoly[e].X == p[0][1].X) && (lPoly[e].Y == p[0][1].Y))
                            {
                                int nextIndex = (e + 1) % lPoly.Count;
                                if ((lPoly[nextIndex].X == p[0][0].X) && (lPoly[nextIndex].Y == p[0][0].Y))
                                {
                                    edgeIsNew = false;
                                }
                            }
                        }
                    }

                    if (edgeIsNew)
                    {
                        newEdges.Add(new Path(p[0]));
                        break;
                    }
                    else
                    {
                    }
                }
            }

            List<GeoLibPoint[]> final = new List<GeoLibPoint[]>();
            if (newEdges.Count > 0)
            {
                // Turn the new edges into cutters and slice. Not terribly elegant and we're relying on rounding to squash notches later.
                ClipperOffset co = new ClipperOffset();
                co.AddPaths(newEdges, JoinType.jtMiter, EndType.etOpenSquare);
                PolyTree tp = new PolyTree();
                co.Execute(ref tp, 1.0);

                Paths cutters = Clipper.ClosedPathsFromPolyTree(tp);

                c.Clear();

                c.AddPath(lPoly, PolyType.ptSubject, true);
                c.AddPath(cutters[0], PolyType.ptClip, true);
                Paths f = new Paths();
                c.Execute(ClipType.ctDifference, f, PolyFillType.pftEvenOdd, PolyFillType.pftEvenOdd);

                final = GeoWrangler.pointsFromPaths(f, scaling);

                final = GeoWrangler.simplify(final);

                final = GeoWrangler.clockwiseAndReorder(final);
            }

            return final;
        }
    }
}
