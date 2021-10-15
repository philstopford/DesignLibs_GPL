using ClipperLib;
using geoLib;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using utility;

namespace geoWrangler
{
    public static partial class GeoWrangler
    {
        public static GeoLibPoint[] Rotate(GeoLibPoint pivot, GeoLibPoint[] pointList, double angleDegree)
        {
            return pRotate(pivot, pointList, angleDegree);
        }

        static GeoLibPoint[] pRotate(GeoLibPoint pivot, GeoLibPoint[] pointList, double angleDegree)
        {
            if (Math.Abs(angleDegree) < 1E-2) // essentially a zero rotation clamp at 0.01 degrees
            {
                return pointList;
            }

            int pointListLength = pointList.Length;
            GeoLibPoint[] returnList = new GeoLibPoint[pointListLength];

#if !GWSINGLETHREADED
            Parallel.For(0, pointListLength, (i) =>
#else
            for (Int32 i = 0; i < pointList.Count(); i++)
#endif
            {
                returnList[i] = Rotate(pivot, pointList[i], angleDegree);
            }
#if !GWSINGLETHREADED
            );
#endif
            return returnList;
        }

        public static GeoLibPointF[] Rotate(GeoLibPointF pivot, GeoLibPointF[] pointList, double angleDegree)
        {
            return pRotate(pivot, pointList, angleDegree);
        }

        static GeoLibPointF[] pRotate(GeoLibPointF pivot, GeoLibPointF[] pointList, double angleDegree)
        {
            if (Math.Abs(angleDegree) < 1E-2) // essentially a zero rotation clamp at 0.01 degrees
            {
                return pointList;
            }

            int pointListLength = pointList.Length;
            GeoLibPointF[] returnList = new GeoLibPointF[pointListLength];

#if !GWSINGLETHREADED
            Parallel.For(0, pointListLength, (i) =>
#else
            for (Int32 i = 0; i < pointList.Count(); i++)
#endif
            {
                returnList[i] = Rotate(pivot, pointList[i], angleDegree);
            }
#if !GWSINGLETHREADED
            );
#endif
            return returnList;
        }

        public static List<GeoLibPointF> Rotate(GeoLibPointF pivot, List<GeoLibPointF> pointList, double angleDegree)
        {
            return pRotate(pivot, pointList, angleDegree);
        }

        static List<GeoLibPointF> pRotate(GeoLibPointF pivot, List<GeoLibPointF> pointList, double angleDegree)
        {
            if (Math.Abs(angleDegree) < 1E-2) // essentially a zero rotation clamp at 0.01 degrees
            {
                return pointList.ToList();
            }

            List<GeoLibPointF> returnList = new List<GeoLibPointF>();

            for (Int32 i = 0; i < pointList.Count(); i++)
            {
                returnList.Add(Rotate(pivot, pointList[i], angleDegree));
            }

            return returnList;
        }

        public static GeoLibPointF Rotate(GeoLibPointF pivot, GeoLibPointF point, double angleDegree)
        {
            return pRotate(pivot, point, angleDegree);
        }

        static GeoLibPointF pRotate(GeoLibPointF pivot, GeoLibPointF point, double angleDegree)
        {
            if (Math.Abs(angleDegree) < 1E-2) // essentially a zero rotation clamp at 0.01 degrees
            {
                return point;
            }

            double angle = Utils.toRadians(angleDegree);
            double x = pivot.X + ((point.X - pivot.X) * Math.Cos(angle) -
                (point.Y - pivot.Y) * Math.Sin(angle));
            double y = pivot.Y + ((point.X - pivot.X) * Math.Sin(angle) +
                (point.Y - pivot.Y) * Math.Cos(angle));

            GeoLibPointF rotated = new GeoLibPointF(x, y);
            return rotated;
        }

        public static GeoLibPoint Rotate(GeoLibPoint pivot, GeoLibPoint point, double angleDegree)
        {
            return pRotate(pivot, point, angleDegree);
        }

        static GeoLibPoint pRotate(GeoLibPoint pivot, GeoLibPoint point, double angleDegree)
        {
            if (Math.Abs(angleDegree) < 1E-2) // essentially a zero rotation clamp at 0.01 degrees
            {
                return point;
            }

            double angle = Utils.toRadians(angleDegree);
            double x = pivot.X + ((point.X - pivot.X) * Math.Cos(angle) -
                (point.Y - pivot.Y) * Math.Sin(angle));
            double y = pivot.Y + ((point.X - pivot.X) * Math.Sin(angle) +
                (point.Y - pivot.Y) * Math.Cos(angle));

            GeoLibPoint rotated = new GeoLibPoint(x, y);
            return rotated;
        }

        public static IntPoint Rotate(IntPoint pivot, IntPoint point, double angleDegree)
        {
            return pRotate(pivot, point, angleDegree);
        }

        static IntPoint pRotate(IntPoint pivot, IntPoint point, double angleDegree)
        {
            if (Math.Abs(angleDegree) < 1E-2) // essentially a zero rotation clamp at 0.01 degrees
            {
                return point;
            }

            double angle = Utils.toRadians(angleDegree);
            double x = pivot.X + ((point.X - pivot.X) * Math.Cos(angle) -
                (point.Y - pivot.Y) * Math.Sin(angle));
            double y = pivot.Y + ((point.X - pivot.X) * Math.Sin(angle) +
                (point.Y - pivot.Y) * Math.Cos(angle));

            IntPoint rotated = new IntPoint(x, y, point.Z);
            return rotated;
        }

        public static GeoLibPointF[] flip(bool H, bool V, bool alignX, bool alignY, GeoLibPointF pivot, GeoLibPointF[] source)
        {
            return pFlip(H, V, alignX, alignY, pivot, source);
        }

        public static List<GeoLibPointF> flip(bool H, bool V, bool alignX, bool alignY, GeoLibPointF pivot, List<GeoLibPointF> source)
        {
            return pFlip(H, V, alignX, alignY, pivot, source.ToArray()).ToList();
        }

        static GeoLibPointF[] pFlip(bool H, bool V, bool alignX, bool alignY, GeoLibPointF pivot, GeoLibPointF[] source)
        {
            int sLength = source.Length;
            GeoLibPointF[] ret = new GeoLibPointF[sLength];
#if !GWSINGLETHREADED
            Parallel.For(0, sLength, (pt) =>
#else
            for (Int32 pt = 0; pt < source.Length; pt++)
#endif
            {
                double newX = source[pt].X;
                if (H)
                {
                    newX = -newX;
                    if (alignX)
                    {
                        newX += (2 * pivot.X);
                    }
                }
                double newY = source[pt].Y;
                if (V)
                {
                    newY = -newY;
                    if (alignY)
                    {
                        newY += (2 * pivot.Y);
                    }
                }
                ret[pt] = new GeoLibPointF(newX, newY);
            }
#if !GWSINGLETHREADED
            );
#endif
            if (H ^ V)
            {
                ret.Reverse(); // preserve ordering.
            }

            return ret;
        }

    }
}
