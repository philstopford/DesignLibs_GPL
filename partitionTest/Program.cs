using System;
using System.Collections.Generic;
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

        static void Main(string[] args)
        {
            // L
            Path lPoly = new Path {
            
            new IntPoint( 0, 0),
            new IntPoint( 0, 50),
            new IntPoint( 10, 50),
            new IntPoint( 10, 20),
            new IntPoint( 60, 20),
            new IntPoint( 60, 0),
            new IntPoint( 0, 0)

            };

            RayCast rc = new RayCast(lPoly, lPoly, 1000, projectCorners: true, invert: true);

            Paths rays = rc.getRays();

            Clipper c = new Clipper();
            c.ZFillFunction = ZFillCallback;

            Paths newEdges = new Paths();

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
                    }
                    else
                    {
                    }
                }
            }

            // Turn the new edges into cutters and slice.
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
            c.Execute(ClipType.ctDifference, f);

            Console.WriteLine("Hello World!");
        }
    }
}
