using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using geoLib;

namespace geoWrangler
{
    public static partial class GeoWrangler
    {
        public static List<List<GeoLibPointF>> makeArray2(List<GeoLibPointF> source, int xCount, double xPitch, int yCount, double yPitch)
        {
            List<GeoLibPointF[]> r = pMakeArray(source.ToArray(), xCount, xPitch, yCount, yPitch);

            List<List<GeoLibPointF>> ret = new List<List<GeoLibPointF>>();
            for (int i = 0; i < r.Count; i++)
            {
                ret.Add(r[i].ToList());
            }

            return ret;
        }

        public static List<GeoLibPointF[]> makeArray(List<GeoLibPointF> source, int xCount, decimal xPitch, int yCount, decimal yPitch)
        {
            return pMakeArray(source.ToArray(), xCount, Convert.ToDouble(xPitch), yCount, Convert.ToDouble(yPitch));
        }

        public static List<GeoLibPointF[]> makeArray(List<GeoLibPointF> source, int xCount, double xPitch, int yCount, double yPitch)
        {
            return pMakeArray(source.ToArray(), xCount, Convert.ToDouble(xPitch), yCount, Convert.ToDouble(yPitch));
        }

        public static List<GeoLibPointF[]> makeArray(GeoLibPointF[] source, int xCount, decimal xPitch, int yCount, decimal yPitch)
        {
            return pMakeArray(source, xCount, Convert.ToDouble(xPitch), yCount, Convert.ToDouble(yPitch));
        }

        public static List<GeoLibPointF[]> makeArray(GeoLibPointF[] source, int xCount, double xPitch, int yCount, double yPitch)
        {
            return pMakeArray(source, xCount, xPitch, yCount, yPitch);
        }

        static List<GeoLibPointF[]> pMakeArray(GeoLibPointF[] source, int xCount, double xPitch, int yCount, double yPitch)
        {
            List<GeoLibPointF[]> ret = new List<GeoLibPointF[]>();
            for (int x = 0; x < xCount; x++)
            {
                for (int y = 0; y < yCount; y++)
                {
                    ret.Add((pMove(source, x * xPitch, y * yPitch)));
                }
            }

            return ret;
        }

        public static List<List<GeoLibPointF>> makeArray(List<List<GeoLibPointF>> source, int xCount, decimal xPitch, int yCount, decimal yPitch)
        {
            return pMakeArray(source, xCount, Convert.ToDouble(xPitch), yCount, Convert.ToDouble(yPitch));
        }

        public static List<List<GeoLibPointF>> makeArray(List<List<GeoLibPointF>> source, int xCount, double xPitch, int yCount, double yPitch)
        {
            return pMakeArray(source, xCount, xPitch, yCount, yPitch);
        }

        static List<List<GeoLibPointF>> pMakeArray(List<List<GeoLibPointF>> source, int xCount, double xPitch, int yCount, double yPitch)
        {
            List<List<GeoLibPointF>> ret = new List<List<GeoLibPointF>>();
            for (int x = 0; x < xCount; x++)
            {
                for (int y = 0; y < yCount; y++)
                {
                    ret.AddRange(pMove(source, x * xPitch, y * yPitch));
                }
            }

            return ret;
        }


        public static List<GeoLibPointF[]> makeArray(List<GeoLibPointF[]> source, int xCount, decimal xPitch, int yCount, decimal yPitch)
        {
            return pMakeArray(source, xCount, Convert.ToDouble(xPitch), yCount, Convert.ToDouble(yPitch));
        }

        public static List<GeoLibPointF[]> makeArray(List<GeoLibPointF[]> source, int xCount, double xPitch, int yCount, double yPitch)
        {
            return pMakeArray(source, xCount, xPitch, yCount, yPitch);
        }

        static List<GeoLibPointF[]> pMakeArray(List<GeoLibPointF[]> source, int xCount, double xPitch, int yCount, double yPitch)
        {
            List<GeoLibPointF[]> ret = new List<GeoLibPointF[]>();
            for (int x = 0; x < xCount; x++)
            {
                for (int y = 0; y < yCount; y++)
                {
                    ret.AddRange(pMove(source, x * xPitch, y * yPitch));
                }
            }

            return ret;
        }
    }
}
