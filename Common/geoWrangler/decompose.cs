using ClipperLib;
using System;
using System.Collections.Generic;
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
    }
}
