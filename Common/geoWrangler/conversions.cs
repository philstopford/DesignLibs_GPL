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
        public static Paths pathsFromPoints(List<GeoLibPoint[]> source, Int64 scaling)
        {
            return pPathsFromPoints(source, scaling);
        }

        static Paths pPathsFromPoints(List<GeoLibPoint[]> source, Int64 scaling)
        {
            Paths ret = new Paths();
            try
            {
                foreach (GeoLibPoint[] t in source)
                {
                    ret.Add(pathFromPoint(t, scaling));
                }
            }
            catch (Exception)
            {
            }
            return ret;
        }

        public static Path pathFromPoint(GeoLibPoint[] source, Int64 scaling)
        {
            return pPathFromPoint(source, scaling);
        }

        static Path pPathFromPoint(GeoLibPoint[] source, Int64 scaling)
        {
            int length = source.Length;
            if ((source[0].X != source[^1].X) && (source[0].Y != source[^1].Y))
            {
                length++; // close the geometry
            }
            Path returnPath = new Path();
            try
            {
                foreach (GeoLibPoint t in source)
                {
                    returnPath.Add(new IntPoint(t.X * scaling, t.Y * scaling));
                }
            }
            catch (Exception)
            {
            }

            // Close the shape
            if (length != source.Length)
            {
                returnPath.Add(new IntPoint(returnPath[0]));
            }
            return returnPath;
        }

        public static List<GeoLibPoint[]> pointsFromPaths(Paths source, Int64 scaling)
        {
            return pPointsFromPaths(source, scaling);
        }

        static List<GeoLibPoint[]> pPointsFromPaths(Paths source, Int64 scaling)
        {
            List<GeoLibPoint[]> ret = new List<GeoLibPoint[]>();
            foreach (Path t in source)
            {
                ret.Add(pPointFromPath(t, scaling));
            }

            return ret;
        }

        public static GeoLibPoint[] pointFromPath(Path source, Int64 scaling)
        {
            return pPointFromPath(source, scaling);
        }

        static GeoLibPoint[] pPointFromPath(Path source, Int64 scaling)
        {
            int length = source.Count;
            int sCount = length;
            if (length > 1)
            {
                if ((source[0].X != source[sCount - 1].X) && (source[0].Y != source[sCount - 1].Y))
                {
                    length++; // close the geometry
                }
            }
            GeoLibPoint[] returnPoint = new GeoLibPoint[length];
#if !GWSINGLETHREADED
            Parallel.For(0, sCount, (pt) =>
#else
            for (int pt = 0; pt < source.Count; pt++)
#endif
            {
                
                returnPoint[pt] = new GeoLibPoint((Int64)Math.Round(Convert.ToDecimal(source[pt].X) / scaling), (Int64)Math.Round(Convert.ToDecimal(source[pt].Y) / scaling));
            }
#if !GWSINGLETHREADED
            );
#endif
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
                foreach (GeoLibPointF[] t in source)
                {
                    ret.Add(pathFromPointF(t, scaling));
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
            if ((Math.Abs(source[0].X - source[^1].X) > Double.Epsilon) && (Math.Abs(source[0].Y - source[^1].Y) > Double.Epsilon))
            {
                length++; // close the geometry
            }
            Path returnPath = new Path();
            try
            {
                foreach (GeoLibPointF t in source)
                {
                    returnPath.Add(new IntPoint(Convert.ToInt64(t.X * scaling),
                        Convert.ToInt64(t.Y * scaling)));
                }
            }
            catch (Exception)
            {
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
            foreach (Path t in source)
            {
                ret.Add(pPointFFromPath(t, scaling));
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
#if !GWSINGLETHREADED
            Parallel.For(0, sourceCount, (pt) =>
#else
            for (int pt = 0; pt < source.Count(); pt++)
#endif
            {
                returnPointF[pt] = new GeoLibPointF((double)(source[pt].X) / scaling,
                                                    (double)(source[pt].Y) / scaling);
            }
#if !GWSINGLETHREADED
            );
#endif
            // Close the shape.
            if (length != sourceCount)
            {
                returnPointF[length - 1] = new GeoLibPointF(returnPointF[0]);
            }
            return returnPointF;
        }

        public static List<GeoLibPoint[]> pointsFromPointFs(List<GeoLibPointF[]> source, Int64 scaling)
        {
            return pPointsFromPointFs(source, scaling);
        }

        static List<GeoLibPoint[]> pPointsFromPointFs(List<GeoLibPointF[]> source, Int64 scaling)
        {
            List<GeoLibPoint[]> ret = new List<GeoLibPoint[]>();
            try
            {
                foreach (GeoLibPointF[] t in source)
                {
                    ret.Add(pPointsFromPointF(t, scaling));
                }
            }
            catch (Exception)
            {
            }
            return ret;
        }

        public static GeoLibPoint[] pointsFromPointF(GeoLibPointF[] source, Int64 scaling)
        {
            return pPointsFromPointF(source, scaling);
        }

        static GeoLibPoint[] pPointsFromPointF(GeoLibPointF[] source, Int64 scaling)
        {
            GeoLibPoint[] ret = new GeoLibPoint[source.Length];
            for (int pt = 0; pt < source.Count(); pt++)
            {
                try
                {
                    ret[pt] = new GeoLibPoint(Convert.ToInt64(source[pt].X * scaling),
                                                  Convert.ToInt64(source[pt].Y * scaling));
                }
                catch (Exception)
                {
                }
            }

            return ret;
        }

        public static List<GeoLibPointF[]> pointFsFromPoints(List<GeoLibPoint[]> source, Int64 scaling)
        {
            return pPointFsFromPoints(source, scaling);
        }

        static List<GeoLibPointF[]> pPointFsFromPoints(List<GeoLibPoint[]> source, Int64 scaling)
        {
            List<GeoLibPointF[]> ret = new List<GeoLibPointF[]>();
            try
            {
                foreach (GeoLibPoint[] t in source)
                {
                    ret.Add(pPointFsFromPoint(t, scaling));
                }
            }
            catch (Exception)
            {
            }
            return ret;
        }

        public static GeoLibPointF[] pointFsFromPoint(GeoLibPoint[] source, Int64 scaling)
        {
            return pPointFsFromPoint(source, scaling);
        }

        static GeoLibPointF[] pPointFsFromPoint(GeoLibPoint[] source, Int64 scaling)
        {
            GeoLibPointF[] ret = new GeoLibPointF[source.Length];
            for (int pt = 0; pt < source.Count(); pt++)
            {
                try
                {
                    ret[pt] = new GeoLibPointF(Convert.ToDouble(source[pt].X) / scaling,
                                                  Convert.ToDouble(source[pt].Y) / scaling);
                }
                catch (Exception)
                {
                }
            }

            return ret;
        }

    }
}
