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
        public static List<GeoLibPointF[]> xySequence(List<GeoLibPointF[]> source)
        {
            return pXYSequence(source);
        }

        static List<GeoLibPointF[]> pXYSequence(List<GeoLibPointF[]> source)
        {
            int sourceCount = source.Count;

            List<GeoLibPointF> midPoints = new List<GeoLibPointF>();

            for (int i = 0; i < sourceCount; i++)
            {
                midPoints.Add(midPoint(source[i]));
                midPoints[i].tag = i; // track our original poly for this midpoint through the re-order
            }

            midPoints = midPoints.OrderBy(p => p.X).ThenBy(p => p.Y).ToList();

            List<GeoLibPointF[]> ret = new List<GeoLibPointF[]>();

            for (int i = 0; i < sourceCount; i++)
            {
                ret.Add(source[midPoints[i].tag].ToArray());
            }

            return ret;
        }

        public static List<GeoLibPoint[]> xySequence(List<GeoLibPoint[]> source)
        {
            return pXYSequence(source);
        }

        static List<GeoLibPoint[]> pXYSequence(List<GeoLibPoint[]> source)
        {
            int sourceCount = source.Count;

            List<GeoLibPointF> midPoints = new List<GeoLibPointF>();

            for (int i = 0; i < sourceCount; i++)
            {
                midPoints.Add(midPoint(source[i]));
                midPoints[i].tag = i; // track our original poly for this midpoint through the re-order
            }

            midPoints = midPoints.OrderBy(p => p.X).ThenBy(p => p.Y).ToList();

            List<GeoLibPoint[]> ret = new List<GeoLibPoint[]>();

            for (int i = 0; i < sourceCount; i++)
            {
                ret.Add(source[midPoints[i].tag].ToArray());
            }

            return ret;
        }

        public static List<GeoLibPointF[]> clockwiseSequence(List<GeoLibPointF[]> source)
        {
            return pClockwiseSequence(source);
        }

        static List<GeoLibPointF[]> pClockwiseSequence(List<GeoLibPointF[]> source)
        {
            int sourceCount = source.Count;
            GeoLibPointF[] midPoints = new GeoLibPointF[sourceCount];

            for (int i = 0; i < sourceCount; i++)
            {
                midPoints[i] = midPoint(source[i]);
            }

            GeoLibPointF[] reorderedMidPoints = pClockwise(midPoints);

            // Position in array is new polygon; value is old index.
            int[] newOrder = new int[source.Count];
            for (int i = 0; i < sourceCount; i++)
            {
                for (int j = 0; j < sourceCount; j++)
                {
                    if ((reorderedMidPoints[i].X == midPoints[j].X) && (reorderedMidPoints[i].Y == midPoints[j].Y))
                    {
                        newOrder[i] = j;
                        break;
                    }
                }
            }

            List<GeoLibPointF[]> ret = new List<GeoLibPointF[]>();
            for (int i = 0; i < sourceCount; i++)
            {
                ret.Add(source[newOrder[i]].ToArray());
            }

            return ret;
        }

        public static List<GeoLibPoint[]> clockwiseSequence(List<GeoLibPoint[]> source)
        {
            return pClockwiseSequence(source);
        }

        static List<GeoLibPoint[]> pClockwiseSequence(List<GeoLibPoint[]> source)
        {
            int sourceCount = source.Count;
            GeoLibPointF[] midPoints = new GeoLibPointF[sourceCount];

            for (int i = 0; i < sourceCount; i++)
            {
                midPoints[i] = midPoint(source[i]);
            }

            GeoLibPointF[] reorderedMidPoints = pClockwise(midPoints);

            // Position in array is new polygon; value is old index.
            int[] newOrder = new int[source.Count];
            for (int i = 0; i < sourceCount; i++)
            {
                for (int j = 0; j < sourceCount; j++)
                {
                    if ((reorderedMidPoints[i].X == midPoints[j].X) && (reorderedMidPoints[i].Y == midPoints[j].Y))
                    {
                        newOrder[i] = j;
                        break;
                    }
                }
            }

            List<GeoLibPoint[]> ret = new List<GeoLibPoint[]>();
            for (int i = 0; i < sourceCount; i++)
            {
                ret.Add(source[newOrder[i]].ToArray());
            }

            return ret;
        }

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
