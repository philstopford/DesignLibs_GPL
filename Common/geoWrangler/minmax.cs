using ClipperLib;
using geoLib;
using System;
using System.Collections.Generic;

namespace geoWrangler
{
    using Path = List<IntPoint>;

    public static partial class GeoWrangler
    {
        public static Int32 MinX(Path sourcePath)
        {
            return pMinX(sourcePath);
        }

        static Int32 pMinX(Path sourcePath)
        {
            Int64 min = sourcePath[0].X;
            Int32 minIndex = 0;

            for (Int32 i = 1; i < sourcePath.Count; ++i)
            {
                if (sourcePath[i].X < min)
                {
                    min = sourcePath[i].X;
                    minIndex = i;
                }
            }

            return minIndex;
        }

        public static Int32 MaxX(Path sourcePath)
        {
            return pMaxX(sourcePath);
        }

        static Int32 pMaxX(Path sourcePath)
        {
            Int64 max = sourcePath[0].X;
            Int32 maxIndex = 0;

            for (Int32 i = 1; i < sourcePath.Count; ++i)
            {
                if (sourcePath[i].X > max)
                {
                    max = sourcePath[i].X;
                    maxIndex = i;
                }
            }

            return maxIndex;
        }

        public static Int32 MinY(Path sourcePath)
        {
            return pMinY(sourcePath);
        }

        static Int32 pMinY(Path sourcePath)
        {
            Int64 min = sourcePath[0].Y;
            Int32 minIndex = 0;

            for (Int32 i = 1; i < sourcePath.Count; ++i)
            {
                if (sourcePath[i].Y < min)
                {
                    min = sourcePath[i].Y;
                    minIndex = i;
                }
            }

            return minIndex;
        }

        public static Int32 MaxY(Path sourcePath)
        {
            return pMaxY(sourcePath);
        }

        static Int32 pMaxY(Path sourcePath)
        {
            Int64 max = sourcePath[0].Y;
            Int32 maxIndex = 0;

            for (Int32 i = 1; i < sourcePath.Count; ++i)
            {
                if (sourcePath[i].Y > max)
                {
                    max = sourcePath[i].Y;
                    maxIndex = i;
                }
            }

            return maxIndex;
        }

        public static Int32 MinX(GeoLibPoint[] iPoints)
        {
            return pMinX(iPoints);
        }

        static Int32 pMinX(GeoLibPoint[] iPoints)
        {
            Int64 min = iPoints[0].X;
            Int32 minIndex = 0;

            for (Int32 i = 1; i < iPoints.Length; ++i)
            {
                if (iPoints[i].X < min)
                {
                    min = iPoints[i].X;
                    minIndex = i;
                }
            }

            return minIndex;
        }

        public static Int32 MinY(GeoLibPoint[] iPoints)
        {
            return pMinY(iPoints);
        }

        static Int32 pMinY(GeoLibPoint[] iPoints)
        {
            Int64 min = iPoints[0].Y;
            Int32 minIndex = 0;

            for (Int32 i = 1; i < iPoints.Length; ++i)
            {
                if (iPoints[i].Y < min)
                {
                    min = iPoints[i].Y;
                    minIndex = i;
                }
            }

            return minIndex;
        }

        public static Int32 MaxX(GeoLibPoint[] iPoints)
        {
            return pMaxX(iPoints);
        }

        static Int32 pMaxX(GeoLibPoint[] iPoints)
        {
            Int64 max = iPoints[0].X;
            Int32 maxIndex = 0;

            for (Int32 i = 1; i < iPoints.Length; ++i)
            {
                if (iPoints[i].X > max)
                {
                    max = iPoints[i].X;
                    maxIndex = i;
                }
            }

            return maxIndex;
        }

        public static Int32 MaxY(GeoLibPoint[] iPoints)
        {
            return pMaxY(iPoints);
        }

        static Int32 pMaxY(GeoLibPoint[] iPoints)
        {
            Int64 max = iPoints[0].Y;
            Int32 maxIndex = 0;

            for (Int32 i = 1; i < iPoints.Length; ++i)
            {
                if (iPoints[i].Y > max)
                {
                    max = iPoints[i].Y;
                    maxIndex = i;
                }
            }

            return maxIndex;
        }

        public static Int32 MinX(GeoLibPointF[] iPoints)
        {
            return pMinX(iPoints);
        }

        static Int32 pMinX(GeoLibPointF[] iPoints)
        {
            double min = iPoints[0].X;
            Int32 minIndex = 0;

            for (Int32 i = 1; i < iPoints.Length; ++i)
            {
                if (iPoints[i].X < min)
                {
                    min = iPoints[i].X;
                    minIndex = i;
                }
            }

            return minIndex;
        }

        public static Int32 MinY(GeoLibPointF[] iPoints)
        {
            return pMinY(iPoints);
        }

        static Int32 pMinY(GeoLibPointF[] iPoints)
        {
            double min = iPoints[0].Y;
            Int32 minIndex = 0;

            for (Int32 i = 1; i < iPoints.Length; ++i)
            {
                if (iPoints[i].Y < min)
                {
                    min = iPoints[i].Y;
                    minIndex = i;
                }
            }

            return minIndex;
        }

        public static Int32 MaxX(GeoLibPointF[] iPoints)
        {
            return pMaxX(iPoints);
        }

        static Int32 pMaxX(GeoLibPointF[] iPoints)
        {
            double max = iPoints[0].X;
            Int32 maxIndex = 0;

            for (Int32 i = 1; i < iPoints.Length; ++i)
            {
                if (iPoints[i].X > max)
                {
                    max = iPoints[i].X;
                    maxIndex = i;
                }
            }

            return maxIndex;
        }

        public static Int32 MaxY(GeoLibPointF[] iPoints)
        {
            return pMaxY(iPoints);
        }

        static Int32 pMaxY(GeoLibPointF[] iPoints)
        {
            double max = iPoints[0].Y;
            Int32 maxIndex = 0;

            for (Int32 i = 1; i < iPoints.Length; ++i)
            {
                if (iPoints[i].Y > max)
                {
                    max = iPoints[i].Y;
                    maxIndex = i;
                }
            }

            return maxIndex;
        }

    }
}
