using System;
using System.Collections.Generic;
using System.Text;
using System.Threading.Tasks;
using geoLib;

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

        public static GeoLibPointF[] move(GeoLibPointF[] source, double x, double y)
        {
            return pMove(source, x, y);
        }

        static GeoLibPointF[] pMove(GeoLibPointF[] source, double x, double y)
        {
            int sLength = source.Length;
            GeoLibPointF[] ret = new GeoLibPointF[sLength];
            Parallel.For(0, sLength, (i) => // for (int i = 0; i < source.Length; i++)
            {
                ret[i] = new GeoLibPointF(source[i].X + x, source[i].Y + y);
            });

            return ret;
        }
    }
}
