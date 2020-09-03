using System;
using System.Collections.Generic;
using System.Linq;
using System.Net.Security;
using ClipperLib;
using geoLib;
using geoWrangler;

namespace partitionTest
{
    using Path = List<IntPoint>;
    using Paths = List<List<IntPoint>>;
    class Program
    {
        static void ZFillCallback(IntPoint bot1, IntPoint top1, IntPoint bot2, IntPoint top2, ref IntPoint pt)
        {
            pt.Z = -1; // Tag our intersection points.
        }

        static List<GeoLibPoint[]> processGeometry(GeoLibPoint[] _poly)
        {
            List<GeoLibPoint[]> ret = new List<GeoLibPoint[]>();
            ret.Add(_poly.ToArray());

            bool changed = true;
            while(changed)
            {
                changed = false;
                int retCount = ret.Count;
                for (int i = 0; i < retCount; i++)
                {
                    List<GeoLibPoint[]> decomp = rect_decomp(ret[i].ToArray());
                    if (decomp.Count > 1)
                    {
                        // We decomposed something and need to start over.
                        ret.RemoveAt(i);
                        ret.AddRange(decomp);
                        changed = true;
                        break;
                    }
                }
            }

            return ret;
        }


        static List<GeoLibPoint[]> rect_decomp(GeoLibPoint[] _poly)
        {
            Path lPoly = GeoWrangler.pathFromPoint(_poly, scaling);

            RayCast rc = new RayCast(lPoly, lPoly, 1000 * scaling, projectCorners: true, invert: true);

            Paths rays = rc.getRays();

            Paths newEdges = new Paths();

            Clipper c = new Clipper();
            c.ZFillFunction = ZFillCallback;

            for (int r = 0; r < rays.Count; r++)
            {
                c.AddPath(rays[r], PolyType.ptSubject, false);
                c.AddPath(lPoly, PolyType.ptClip, true);

                PolyTree pt = new PolyTree();

                c.Execute(ClipType.ctIntersection, pt);
                c.Clear();

                Paths p = Clipper.OpenPathsFromPolyTree(pt);

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

                /*
                c.Clear();
                c.AddPaths(cutters, PolyType.ptSubject, true);
                c.Execute(ClipType.ctUnion, cutters);
                */
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

        static int scaling = 10000;

        static void Main(string[] args)
        {
            // L
            GeoLibPoint[] L = new GeoLibPoint[] {

            new GeoLibPoint( 0, 0),
            new GeoLibPoint( 0, 50),
            new GeoLibPoint( 10, 50),
            new GeoLibPoint( 10, 20),
            new GeoLibPoint( 60, 20),
            new GeoLibPoint( 60, 0),
            new GeoLibPoint( 0, 0)

            };

            GeoLibPoint[] rL = new GeoLibPoint[] {

            new GeoLibPoint( 0, 0),
            new GeoLibPoint( 0, 20),
            new GeoLibPoint( 10, 20),
            new GeoLibPoint( 10, 50),
            new GeoLibPoint( 60, 50),
            new GeoLibPoint( 60, 0),
            new GeoLibPoint( 0, 0)

            };



            // U
            GeoLibPoint[] U = new GeoLibPoint[] {

            new GeoLibPoint( 0, 0),
            new GeoLibPoint( 0, 50),
            new GeoLibPoint( 10, 50),
            new GeoLibPoint( 10, 20),
            new GeoLibPoint( 60, 20),
            new GeoLibPoint( 60, 80),
            new GeoLibPoint( 120, 80),
            new GeoLibPoint( 120, 0),
            new GeoLibPoint( 0, 0)

            };

            // T
            GeoLibPoint[] T = new GeoLibPoint[] {

            new GeoLibPoint( 0, 50),
            new GeoLibPoint( 0, 80),
            new GeoLibPoint( 80, 80),
            new GeoLibPoint( 80, 50),
            new GeoLibPoint( 60, 50),
            new GeoLibPoint( 60, 0),
            new GeoLibPoint( 40, 0),
            new GeoLibPoint( 40, 50),
            new GeoLibPoint( 0, 50)

            };


            // X
            GeoLibPoint[] X = new GeoLibPoint[] {

            new GeoLibPoint( 0, 50),
            new GeoLibPoint( 0, 80),
            new GeoLibPoint( 60, 80),
            new GeoLibPoint( 60, 100),
            new GeoLibPoint( 80, 100),
            new GeoLibPoint( 80, 80),
            new GeoLibPoint( 100, 80),
            new GeoLibPoint( 100, 50),
            new GeoLibPoint( 80, 50),
            new GeoLibPoint( 80, 20),
            new GeoLibPoint( 60, 20),
            new GeoLibPoint( 60, 50),
            new GeoLibPoint( 0, 50)

            };


            // S
            GeoLibPoint[] S = new GeoLibPoint[] {

            new GeoLibPoint( 0, 0),
            new GeoLibPoint( 0, 20),
            new GeoLibPoint( 20, 20),
            new GeoLibPoint( 20, 50),
            new GeoLibPoint( 0, 50),
            new GeoLibPoint( 0, 110),
            new GeoLibPoint( 100, 110),
            new GeoLibPoint( 100, 80),
            new GeoLibPoint( 80, 80),
            new GeoLibPoint( 80, 60),
            new GeoLibPoint( 100, 60),
            new GeoLibPoint( 100, 0),
            new GeoLibPoint( 0, 0)

            };

            List<GeoLibPoint[]> l = processGeometry(L);

            List<GeoLibPoint[]> rl = processGeometry(rL);

            List<GeoLibPoint[]> u = processGeometry(U);

            List<GeoLibPoint[]> t = processGeometry(T);

            List<GeoLibPoint[]> x = processGeometry(X);

            List<GeoLibPoint[]> s = processGeometry(S);

            Console.WriteLine("Hello World!");
        }
    }
}
