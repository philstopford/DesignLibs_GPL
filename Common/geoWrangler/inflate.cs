using ClipperLib2;
using geoLib;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace geoWrangler;

using Path = List<Point64>;
using Paths = List<List<Point64>>;
using PathsD = List<List<PointD>>;

public static partial class GeoWrangler
{
    public static GeoLibPoint[] inflatePath(GeoLibPoint[] source, int width)
    {
        return pInflatePath(source, width);
    }

    private static GeoLibPoint[] pInflatePath(GeoLibPoint[] source, int width)
    {
        switch (width)
        {
            case 0:
                return source;
        }

        Paths allSolutions = new();

        for (int i = 0; i < source.Length - 1; i++)
        {
            ClipperOffset co = new();
            Path o = new()
            {
                new Point64(source[i].X, source[i].Y), new Point64(source[i + 1].X, source[i + 1].Y)
            };
            co.AddPath(o, JoinType.Miter, EndType.Joined);

            int offsetVal = width / 2;

            Paths solution = ClipperFunc.Paths (co.Execute(offsetVal));

            allSolutions.Add(new Path(solution[0]));

            // Need to add a patch polygon to link the segments.
            Path patchPoly = new()
            {
                new Point64(source[i + 1].X - offsetVal, source[i + 1].Y - offsetVal),
                new Point64(source[i + 1].X - offsetVal, source[i + 1].Y + offsetVal),
                new Point64(source[i + 1].X + offsetVal, source[i + 1].Y + offsetVal),
                new Point64(source[i + 1].X + offsetVal, source[i + 1].Y - offsetVal)
            };

            allSolutions.Add(new Path(patchPoly));
        }

        Clipper c = new();
        c.AddSubject(allSolutions);

        Rect64 b = ClipperFunc.GetBounds(allSolutions);

        Path bPath = new()
        {
            new Point64(b.left, b.bottom),
            new Point64(b.left, b.top),
            new Point64(b.right, b.top),
            new Point64(b.right, b.bottom)
        };

        c.AddClip(new Paths { bPath });

        Paths union = new();
        PolyTree pt = new();
        c.Execute(ClipType.Intersection, FillRule.Positive, pt);
        union = ClipperFunc.PolyTreeToPaths(pt);

        GeoLibPoint[] ret;
        if (union.Any())
        {
            // We should only have one result.
            union[0].Add(new Point64(union[0][0])); // force a close - it wasn't done in the Boolean.
            ret = pointFromPath(union[0], 1);
        }
        else
        {
            ret = new [] { new GeoLibPoint(0, 0) };
        }

        return ret;
    }

    public static GeoLibPointF[] resize(GeoLibPointF[] source, double factor)
    {
        return pResize(source, factor);
    }

    private static GeoLibPointF[] pResize(GeoLibPointF[] source, double factor)
    {
        int sLength = source.Length;
        GeoLibPointF[] ret = new GeoLibPointF[sLength];
#if !GWSINGLETHREADED
        Parallel.For(0, sLength, pt =>
#else
            for (int pt = 0; pt < sLength; pt++)
#endif
            {
                ret[pt] = new GeoLibPointF(source[pt].X * factor, source[pt].Y * factor);
            }
#if !GWSINGLETHREADED
        );
#endif
        return ret;
    }

    public static GeoLibPoint[] resize_to_int(GeoLibPointF[] source, double factor)
    {
        return pResize_to_int(source, factor);
    }

    private static GeoLibPoint[] pResize_to_int(GeoLibPointF[] source, double factor)
    {
        int sLength = source.Length;
        GeoLibPoint[] ret = new GeoLibPoint[sLength];

#if !GWSINGLETHREADED
        Parallel.For(0, sLength, pt =>
#else
            for (int pt = 0; pt < sLength; pt++)
#endif
            {
                ret[pt] = new GeoLibPoint((int)(source[pt].X * factor), (int)(source[pt].Y * factor));
            }
#if !GWSINGLETHREADED
        );
#endif
        return ret;
    }

    public static GeoLibPoint[] resize(GeoLibPoint[] source, double factor)
    {
        return pResize(source, factor);
    }

    private static GeoLibPoint[] pResize(GeoLibPoint[] source, double factor)
    {
        int sLength = source.Length;
        GeoLibPoint[] ret = new GeoLibPoint[sLength];
#if !GWSINGLETHREADED
        Parallel.For(0, sLength, pt =>
#else
            for (int pt = 0; pt < sLength; pt++)
#endif
            {
                ret[pt] = new GeoLibPoint(source[pt].X * factor, source[pt].Y * factor);
            }
#if !GWSINGLETHREADED
        );
#endif
        return ret;
    }

    public static GeoLibPoint[] resize(GeoLibPoint pivot, GeoLibPoint[] source, double factor)
    {
        return pResize(pivot, source, factor);
    }

    private static GeoLibPoint[] pResize(GeoLibPoint pivot, GeoLibPoint[] source, double factor)
    {
        GeoLibPoint[] pointarray = new GeoLibPoint[source.Length];
#if !GWSINGLETHREADED
        Parallel.For(0, pointarray.Length, i => 
#else
            for (int i = 0; i < pointarray.Length; i++)
#endif
            {
                pointarray[i] = new GeoLibPoint(pivot.X + (source[i].X - pivot.X) * factor, pivot.Y + (source[i].Y - pivot.Y) * factor);
            }
#if !GWSINGLETHREADED
        );
#endif
        return pointarray;
    }
}