using ClipperLib;
using geoLib;
using System.Collections.Generic;
using System.Linq;

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
                ret = pointFromPath(union[0]);
            }
            else
            {
                ret = new GeoLibPoint[1] { new GeoLibPoint(0, 0) };
            }

            return ret;
        }
    }
}
