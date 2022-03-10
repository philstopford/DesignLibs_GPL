using System;
using ClipperLib2;
using geoLib;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace geoWrangler;

using Path = List<Point64>;
using Paths = List<List<Point64>>;

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
        
        ClipperOffset co = new();
        Path a = GeoWrangler.pathFromPoint(source, 1);
        // Path from Point auto-closes the input for historical reasons. We may not want this....
        if (pDistanceBetweenPoints(source[0], source[^1]) > Double.Epsilon)
        {
            pStripTerminators(a, false);
        }
        co.AddPath(a, JoinType.Miter, EndType.Square);
        Paths output = ClipperFunc.Paths(co.Execute(width));

        return pPointFromPath(pClose(output[0]), 1);

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