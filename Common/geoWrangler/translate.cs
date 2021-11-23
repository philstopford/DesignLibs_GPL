using System;
using System.Collections.Generic;
using geoLib;
using System.Linq;
using System.Threading.Tasks;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static List<GeoLibPointF[]> move(List<GeoLibPointF[]> source, decimal x, decimal y)
    {
        return pMove(source, Convert.ToDouble(x), Convert.ToDouble(y));
    }
    public static List<GeoLibPointF[]> move(List<GeoLibPointF[]> source, double x, double y)
    {
        return pMove(source, x, y);
    }

    private static List<GeoLibPointF[]> pMove(List<GeoLibPointF[]> source, double x, double y)
    {
        return source.Select(t => pMove(t, x, y)).ToList();
    }

    private static List<List<GeoLibPointF>> pMove(List<List<GeoLibPointF>> source, double x, double y)
    {
        return source.Select(t => pMove(t.ToArray(), x, y).ToList()).ToList();
    }

    public static GeoLibPointF[] move(GeoLibPointF[] source, decimal x, decimal y)
    {
        return pMove(source, Convert.ToDouble(x), Convert.ToDouble(y));
    }

    public static GeoLibPointF[] move(GeoLibPointF[] source, double x, double y)
    {
        return pMove(source, x, y);
    }

    private static GeoLibPointF[] pMove(GeoLibPointF[] source, double x, double y)
    {
        int sLength = source.Length;
        GeoLibPointF[] ret = new GeoLibPointF[sLength];
#if !GWSINGLETHREADED
        Parallel.For(0, sLength, i =>
#else
            for (int i = 0; i < source.Length; i++)
#endif
            {
                ret[i] = new GeoLibPointF(source[i].X + x, source[i].Y + y);
            }
#if !GWSINGLETHREADED
        );
#endif
        return ret;
    }

    public static List<GeoLibPointF> move(List<GeoLibPointF> source, double x, double y)
    {
        return pMove(source.ToArray(), x, y).ToList();
    }
}