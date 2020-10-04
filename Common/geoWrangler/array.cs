using System;
using System.Collections.Generic;
using System.Text;
using geoLib;

namespace geoWrangler
{
    public static partial class GeoWrangler
    {
        public static List<GeoLibPointF[]> makeArray(GeoLibPointF[] source, int xCount, decimal xPitch, int yCount, decimal yPitch)
        {
            return pMakeArray(source, xCount, Convert.ToDouble(xPitch), yCount, Convert.ToDouble(yPitch));
        }

        public static List<GeoLibPointF[]> makeArray(GeoLibPointF[] source, int xCount, double xPitch, int yCount, double yPitch)
        {
            return pMakeArray(source, xCount, xPitch, yCount, yPitch);
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
            for (int i = 0; i < source.Count; i++)
            {
                ret.AddRange(pMakeArray(source[i], xCount, xPitch, yCount, yPitch));
            }

            return ret;
        }


        static List<GeoLibPointF[]> pMakeArray(GeoLibPointF[] source, int xCount, double xPitch, int yCount, double yPitch)
        {
            List<GeoLibPointF[]> ret = new List<GeoLibPointF[]>();
            for (int x = 0; x < xCount; x++)
            {
                for (int y = 0; y < yCount; y++)
                {
                    ret.Add(pMove(source, x * xPitch, y * yPitch));
                }
            }

            return ret;
        }
    }
}
