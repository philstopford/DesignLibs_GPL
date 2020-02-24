using ClipperLib;
using geoLib;
using System;
using System.Collections.Generic;
using System.Linq;

namespace geoWrangler
{
    using Path = List<IntPoint>;
    using Paths = List<List<IntPoint>>;

    public static partial class GeoWrangler
    {
        public static GeoLibPoint[] clockwise(GeoLibPoint[] points)
        {
            return pClockwise(points);
        }

        static GeoLibPoint[] pClockwise(GeoLibPoint[] points)
        {
            if (!pIsClockwise(points))
            {
                Array.Reverse(points);
            }

            return points;
        }

        public static GeoLibPointF[] clockwise(GeoLibPointF[] iPoints)
        {
            return pClockwise(iPoints);
        }

        static GeoLibPointF[] pClockwise(GeoLibPointF[] iPoints)
        {
            if (!pIsClockwise(iPoints))
            {
                Array.Reverse(iPoints);
            }

            return iPoints;
        }

        public static Paths clockwise(Paths source)
        {
            Paths ret = new Paths();
            for (int i = 0; i < source.Count; i++)
            {
                ret.Add(pClockwise(source[i].ToList()));
            }

            return ret;
        }

        public static Path clockwise(Path iPoints)
        {
            return pClockwise(iPoints);
        }

        static Path pClockwise(Path iPoints)
        {
            if (!pIsClockwise(iPoints))
            {
                Array.Reverse(iPoints.ToArray());
            }

            return iPoints;
        }
    }
}
