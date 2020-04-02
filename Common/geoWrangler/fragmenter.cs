using ClipperLib;
using geoLib;
using System;
using System.Collections.Generic;
using System.Linq;

namespace geoWrangler
{
    using Path = List<IntPoint>;
    using Paths = List<List<IntPoint>>;
    public class Fragmenter
    {
        double resolution;
        double scaleFactor;

        const double def_res = 1.0;
        const double def_scale = 1E4;

        public Fragmenter(double res = def_res, double scaling = def_scale)
        {
            pFragmenter(res, scaling);
        }

        void pFragmenter(double res = def_res, double scaling = def_scale)
        {
            resolution = res;
            scaleFactor = scaling;
        }

        public List<GeoLibPointF[]> fragmentPaths(List<GeoLibPointF[]> source)
        {
            return pFragmentPaths(source);
        }

        List<GeoLibPointF[]> pFragmentPaths(List<GeoLibPointF[]> source)
        {
            if (source.Count == 0)
            {
                return source.ToList();
            }

            List<GeoLibPointF[]> ret = new List<GeoLibPointF[]>();
            for (int i = 0; i < source.Count; i++)
            {
                ret.Add(pFragmentPath(source[i].ToArray()));
            }

            return ret.ToList();
        }

        public GeoLibPointF[] fragmentPath(GeoLibPointF[] pointList)
        {
            return pFragmentPath(pointList);
        }

        GeoLibPointF[] pFragmentPath(GeoLibPointF[] pointList)
        {
            List<GeoLibPointF> returnList = fragmentPath(pointList.ToList());

            return returnList.ToArray();
        }

        // Returns start and end point with fragmented in between.
        public List<GeoLibPointF> fragmentPath(List<GeoLibPointF> pointList)
        {
            return pFragmentPath(pointList);
        }

        public Paths fragmentPaths(Paths source)
        {
            return pFragmentPaths(source);
        }

        Paths pFragmentPaths(Paths source)
        {
            Paths ret = new Paths();
            for (int i = 0; i < source.Count; i++)
            {
                ret.Add(pFragmentPath(source[i]));
            }

            return ret;
        }

        public Path fragmentPath(Path source)
        {
            return pFragmentPath(source);
        }

        Path pFragmentPath(Path source)
        {
            Path ret = new Path();
            bool closed = false;
            if ((source[0].X == source[source.Count - 1].X) && (source[0].Y == source[source.Count - 1].Y))
            {
                closed = true;
            }
            // Decouple the geometry.
            Path t = new Path();
            for (int p = 0; p < source.Count; p++)
            {
                t.Add(new IntPoint(source[p].X, source[p].Y));
            }
            if (closed)
            {
                t = GeoWrangler.stripTerminators(t, false);
            }

            for (int p = 0; p < t.Count - 1; p++)
            {
                ret.Add(new IntPoint(t[p]));
                ret.AddRange(pFragmentPath(t[p], t[p + 1]));
            }

            if (closed)
            {
                ret.Add(new IntPoint(t[t.Count - 1]));
                ret.AddRange(pFragmentPath(t[t.Count - 1], t[0]));
                ret = GeoWrangler.close(ret);
            }

            return ret;
        }

        public Path fragmentPath(IntPoint pt1, IntPoint pt2)
        {
            return pFragmentPath(pt1, pt2);
        }

        Path pFragmentPath(IntPoint pt1, IntPoint pt2, bool startAndEndPoints = false)
        {
            Path returnPath = new Path();
            Int64 x_Distance = pt2.X - pt1.X;
            Int64 y_Distance = pt2.Y - pt1.Y;
            // We do some manipulation here because the points in the call are scaled by CentralProperties.scaleFactorForOperation.
            // That would lead to too many points being injected. We therefore scale resolution by 1E4
            Int32 fragmentCount = Convert.ToInt32(Math.Floor(GeoWrangler.distanceBetweenPoints(pt1, pt2) / (resolution * (scaleFactor / 1E4))));

            if (startAndEndPoints)
            {
                returnPath.Add(new IntPoint(pt1));
            }
            if (fragmentCount > 0)
            {
                double x_Step = x_Distance / fragmentCount;
                double y_Step = y_Distance / fragmentCount;

                for (Int32 i = 1; i < fragmentCount; i++)
                {
                    returnPath.Add(new IntPoint(pt1.X + (i * x_Step), pt1.Y + (i * y_Step)));
                }
            }
            if (startAndEndPoints)
            {
                returnPath.Add(new IntPoint(pt2));
            }
            return returnPath;
        }

        public List<GeoLibPointF> fragmentPath(List<GeoLibPointF> pointList, double res, double scaling)
        {
            resolution = res;
            scaleFactor = scaling;
            return pFragmentPath(pointList);
        }

        List<GeoLibPointF> pFragmentPath(List<GeoLibPointF> pointList)
        {
            List<GeoLibPointF> returnList = new List<GeoLibPointF>();
            for (Int32 pt = 0; pt < pointList.Count(); pt++)
            {
                returnList.Add(pointList[pt]);
                if (pt != (pointList.Count() - 1))
                {
                    // Fragment path doesn't return start and end points - just points between.
                    List<GeoLibPointF> newPtsList = fragmentPath(pointList[pt], pointList[pt + 1]);
                    returnList.AddRange(newPtsList);
                }
            }
            return returnList;
        }

        // Internals that don't return start and end points.

        List<GeoLibPointF> fragmentPath(GeoLibPointF pt1, GeoLibPointF pt2)
        {
            List<GeoLibPointF> returnList = new List<GeoLibPointF>();
            double x_Distance = pt2.X - pt1.X;
            double y_Distance = pt2.Y - pt1.Y;
            double fragmentCount = Math.Floor(GeoWrangler.distanceBetweenPoints(pt1, pt2) / (resolution));

            double x_Step = x_Distance / fragmentCount;
            double y_Step = y_Distance / fragmentCount;

            // Avoid sending back the first and last vertices of the path.
            for (Int32 i = 1; i < fragmentCount; i++)
            {
                returnList.Add(new GeoLibPointF((pt1.X + (i * x_Step)), (pt1.Y + (i * y_Step))));
            }
            return returnList;
        }
    }
}
