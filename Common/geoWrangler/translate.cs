using System;
using System.Collections.Generic;
using geoLib;
using System.Linq;
using System.Threading.Tasks;

namespace geoWrangler
{
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

        static List<GeoLibPointF[]> pMove(List<GeoLibPointF[]> source, double x, double y)
        {
            List<GeoLibPointF[]> ret = new List<GeoLibPointF[]>();
            for (int i = 0; i < source.Count; i++)
            {
                ret.Add(pMove(source[i], x, y));
            }

            return ret;
        }

        static List<List<GeoLibPointF>> pMove(List<List<GeoLibPointF>> source, double x, double y)
        {
            List<List<GeoLibPointF>> ret = new List<List<GeoLibPointF>>();
            for (int i = 0; i < source.Count; i++)
            {
                ret.Add(pMove(source[i].ToArray(), x, y).ToList());
            }

            return ret;
        }

        public static GeoLibPointF[] move(GeoLibPointF[] source, decimal x, decimal y)
        {
            return pMove(source, Convert.ToDouble(x), Convert.ToDouble(y));
        }

        public static GeoLibPointF[] move(GeoLibPointF[] source, double x, double y)
        {
            return pMove(source, x, y);
        }

        static GeoLibPointF[] pMove(GeoLibPointF[] source, double x, double y)
        {
            int sLength = source.Length;
            GeoLibPointF[] ret = new GeoLibPointF[sLength];
#if GWTHREADED
            Parallel.For(0, sLength, (i) =>
#else
            for (int i = 0; i < source.Length; i++)
#endif
            {
                ret[i] = new GeoLibPointF(source[i].X + x, source[i].Y + y);
            }
#if GWTHREADED
            );
#endif
            return ret;
        }

        public static List<GeoLibPointF> move(List<GeoLibPointF> source, double x, double y)
        {
            return pMove(source.ToArray(), x, y).ToList();
        }
    }
}
