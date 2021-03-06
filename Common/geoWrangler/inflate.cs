﻿using ClipperLib;
using geoLib;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace geoWrangler
{
    using Path = List<IntPoint>;
    using Paths = List<List<IntPoint>>;

    public static partial class GeoWrangler
    {
        public static GeoLibPoint[] inflatePath(GeoLibPoint[] source, int width)
        {
            return pInflatePath(source, width);
        }

        static GeoLibPoint[] pInflatePath(GeoLibPoint[] source, int width)
        {
            if (width == 0)
            {
                return source;
            }

            Paths allSolutions = new Paths();

            for (int i = 0; i < source.Length - 1; i++)
            {
                ClipperOffset co = new ClipperOffset();
                Path o = new Path();
                o.Add(new IntPoint(source[i].X, source[i].Y));
                o.Add(new IntPoint(source[i + 1].X, source[i + 1].Y));
                co.AddPath(o, JoinType.jtMiter, EndType.etClosedLine);

                int offsetVal = width / 2;

                Paths solution = new Paths();

                co.Execute(ref solution, offsetVal);

                allSolutions.Add(new Path(solution[0]));

                // Need to add a patch polygon to link the segments.
                Path patchPoly = new Path();
                patchPoly.Add(new IntPoint(source[i + 1].X - offsetVal, source[i + 1].Y - offsetVal));
                patchPoly.Add(new IntPoint(source[i + 1].X - offsetVal, source[i + 1].Y + offsetVal));
                patchPoly.Add(new IntPoint(source[i + 1].X + offsetVal, source[i + 1].Y + offsetVal));
                patchPoly.Add(new IntPoint(source[i + 1].X + offsetVal, source[i + 1].Y - offsetVal));

                allSolutions.Add(new Path(patchPoly));
            }

            Clipper c = new Clipper();
            c.AddPaths(allSolutions, PolyType.ptSubject, true);

            IntRect b = Clipper.GetBounds(allSolutions);

            Path bPath = new Path();
            bPath.Add(new IntPoint(b.left, b.bottom));
            bPath.Add(new IntPoint(b.left, b.top));
            bPath.Add(new IntPoint(b.right, b.top));
            bPath.Add(new IntPoint(b.right, b.bottom));

            c.AddPaths(new Paths() { bPath }, PolyType.ptClip, true);

            Paths union = new Paths();
            c.Execute(ClipType.ctIntersection, union, PolyFillType.pftPositive);

            int test = union.Count();

            GeoLibPoint[] ret;
            if (union.Count() > 0)
            {
                // We should only have one result.
                union[0].Add(new IntPoint(union[0][0])); // force a close - it wasn't done in the Boolean.
                ret = pointFromPath(union[0], 1);
            }
            else
            {
                ret = new GeoLibPoint[1] { new GeoLibPoint(0, 0) };
            }

            return ret;
        }

        public static GeoLibPointF[] resize(GeoLibPointF[] source, double factor)
        {
            return pResize(source, factor);
        }
        static GeoLibPointF[] pResize(GeoLibPointF[] source, double factor)
        {
            int sLength = source.Length;
            GeoLibPointF[] ret = new GeoLibPointF[sLength];
#if GWTHREADED
            Parallel.For(0, sLength, (pt) =>
#else
            for (int pt = 0; pt < sLength; pt++)
#endif
            {
                ret[pt] = new GeoLibPointF(source[pt].X * factor, source[pt].Y * factor);
            }
#if GWTHREADED
            );
#endif
            return ret;
        }

        public static GeoLibPoint[] resize_to_int(GeoLibPointF[] source, double factor)
        {
            return pResize_to_int(source, factor);
        }
        static GeoLibPoint[] pResize_to_int(GeoLibPointF[] source, double factor)
        {
            int sLength = source.Length;
            GeoLibPoint[] ret = new GeoLibPoint[sLength];

#if GWTHREADED
            Parallel.For(0, sLength, (pt) =>
#else
            for (int pt = 0; pt < sLength; pt++)
#endif
            {
                ret[pt] = new GeoLibPoint((int)(source[pt].X * factor), (int)(source[pt].Y * factor));
            }
#if GWTHREADED
            );
#endif
            return ret;
        }

        public static GeoLibPoint[] resize(GeoLibPoint[] source, double factor)
        {
            return pResize(source, factor);
        }
        static GeoLibPoint[] pResize(GeoLibPoint[] source, double factor)
        {
            int sLength = source.Length;
            GeoLibPoint[] ret = new GeoLibPoint[sLength];
#if GWTHREADED
            Parallel.For(0, sLength, (pt) =>
#else
            for (int pt = 0; pt < sLength; pt++)
#endif
            {
                ret[pt] = new GeoLibPoint(source[pt].X * factor, source[pt].Y * factor);
            }
#if GWTHREADED
            );
#endif
            return ret;
        }

        public static GeoLibPoint[] resize(GeoLibPoint pivot, GeoLibPoint[] source, double factor)
        {
            return pResize(pivot, source, factor);
        }

        static GeoLibPoint[] pResize(GeoLibPoint pivot, GeoLibPoint[] source, double factor)
        {
            GeoLibPoint[] pointarray = new GeoLibPoint[source.Length];
#if GWTHREADED
            Parallel.For(0, pointarray.Length, (i) => 
#else
            for (int i = 0; i < pointarray.Length; i++)
#endif
            {
                pointarray[i] = new GeoLibPoint(pivot.X + ((source[i].X - pivot.X) * factor), pivot.Y + ((source[i].Y - pivot.Y) * factor));
            }
#if GWTHREADED
            );
#endif
            return pointarray;
        }
    }
}
