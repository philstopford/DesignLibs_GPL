using ClipperLib;
using geoLib;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace geoWrangler
{
    using Path = List<IntPoint>;
    using Paths = List<List<IntPoint>>;

    public static partial class GeoWrangler
    {
        public static Paths pathsFromPoints(List<GeoLibPoint[]> source)
        {
            return pPathsFromPoints(source);
        }

        static Paths pPathsFromPoints(List<GeoLibPoint[]> source)
        {
            Paths ret = new Paths();
            try
            {
                for (int poly = 0; poly < source.Count; poly++)
                {
                    ret.Add(pathFromPoint(source[poly]));
                }
            }
            catch (Exception)
            {
            }
            return ret;
        }

        public static Path pathFromPoint(GeoLibPoint[] source)
        {
            return pPathFromPoint(source);
        }

        static Path pPathFromPoint(GeoLibPoint[] source)
        {
            int length = source.Length;
            if ((source[0].X != source[source.Length - 1].X) && (source[0].Y != source[source.Length - 1].Y))
            {
                length++; // close the geometry
            }
            Path returnPath = new Path();
            for (int pt = 0; pt < source.Count(); pt++)
            {
                try
                {
                    returnPath.Add(new IntPoint(source[pt].X, source[pt].Y));
                }
                catch (Exception)
                {
                }
            }

            // Close the shape
            if (length != source.Length)
            {
                returnPath.Add(new IntPoint(returnPath[0]));
            }
            return returnPath;
        }

        public static List<GeoLibPoint[]> pointsFromPaths(Paths source)
        {
            return pPointsFromPaths(source);
        }

        static List<GeoLibPoint[]> pPointsFromPaths(Paths source)
        {
            List<GeoLibPoint[]> ret = new List<GeoLibPoint[]>();
            for (int poly = 0; poly < source.Count; poly++)
            {
                ret.Add(pPointFromPath(source[poly]));
            }

            return ret;
        }

        public static GeoLibPoint[] pointFromPath(Path source)
        {
            return pPointFromPath(source);
        }

        static GeoLibPoint[] pPointFromPath(Path source)
        {
            int length = source.Count;
            int sCount = length;
            if ((source[0].X != source[sCount - 1].X) && (source[0].Y != source[sCount - 1].Y))
            {
                length++; // close the geometry
            }
            GeoLibPoint[] returnPoint = new GeoLibPoint[length];
            Parallel.For(0, sCount, (pt) => // (int pt = 0; pt < source.Count; pt++)
            {
                returnPoint[pt] = new GeoLibPoint(source[pt].X, source[pt].Y);
            });

            // Close the shape.
            if (length != sCount)
            {
                returnPoint[length - 1] = new GeoLibPoint(returnPoint[0]);
            }
            return returnPoint;
        }

        public static Paths pathsFromPointFs(List<GeoLibPointF[]> source, Int64 scaling)
        {
            return pPathsFromPointFs(source, scaling);
        }

        static Paths pPathsFromPointFs(List<GeoLibPointF[]> source, Int64 scaling)
        {
            Paths ret = new Paths();
            try
            {
                for (int poly = 0; poly < source.Count; poly++)
                {
                    ret.Add(pathFromPointF(source[poly], scaling));
                }
            }
            catch (Exception)
            {
            }
            return ret;
        }

        public static Path pathFromPointF(GeoLibPointF[] source, Int64 scaling)
        {
            return pPathFromPointF(source, scaling);
        }

        static Path pPathFromPointF(GeoLibPointF[] source, Int64 scaling)
        {
            int length = source.Length;
            if ((source[0].X != source[source.Length - 1].X) && (source[0].Y != source[source.Length - 1].Y))
            {
                length++; // close the geometry
            }
            Path returnPath = new Path();
            for (int pt = 0; pt < source.Count(); pt++)
            {
                try
                {
                    returnPath.Add(new IntPoint(Convert.ToInt64(source[pt].X * scaling),
                                                  Convert.ToInt64(source[pt].Y * scaling)));
                }
                catch (Exception)
                {
                }
            }

            // Close the shape
            if (length != source.Length)
            {
                returnPath.Add(new IntPoint(returnPath[0]));
            }
            return returnPath;
        }

        public static List<GeoLibPointF[]> pointFsFromPaths(Paths source, Int64 scaling)
        {
            return pPointFsFromPaths(source, scaling);
        }

        static List<GeoLibPointF[]> pPointFsFromPaths(Paths source, Int64 scaling)
        {
            List<GeoLibPointF[]> ret = new List<GeoLibPointF[]>();
            for (int poly = 0; poly < source.Count; poly++)
            {
                ret.Add(pPointFFromPath(source[poly], scaling));
            }

            return ret;
        }

        public static GeoLibPointF[] pointFFromPath(Path source, Int64 scaling)
        {
            return pPointFFromPath(source, scaling);
        }

        static GeoLibPointF[] pPointFFromPath(Path source, Int64 scaling)
        {
            int length = source.Count;
            int sourceCount = length;
            if ((source[0].X != source[length - 1].X) && (source[0].Y != source[length - 1].Y))
            {
                length++; // close the geometry
            }
            GeoLibPointF[] returnPointF = new GeoLibPointF[length];
            Parallel.For(0, sourceCount, (pt) => // for (int pt = 0; pt < source.Count(); pt++)
            {
                returnPointF[pt] = new GeoLibPointF((double)(source[pt].X) / scaling,
                                                    (double)(source[pt].Y) / scaling);
            });

            // Close the shape.
            if (length != sourceCount)
            {
                returnPointF[length - 1] = new GeoLibPointF(returnPointF[0]);
            }
            return returnPointF;
        }
    }
}
