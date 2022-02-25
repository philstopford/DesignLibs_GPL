using ClipperLib2;
using geoLib;
using System;
using System.Collections.Generic;
using System.Linq;

namespace geoWrangler;

using Path = List<Point64>;
using Paths = List<List<Point64>>;

public static partial class GeoWrangler
{
    public static List<GeoLibPointF[]> xySequence(List<GeoLibPointF[]> source, bool useMidPoint = false)
    {
        return pXYSequence(source, useMidPoint);
    }

    private static List<GeoLibPointF[]> pXYSequence(List<GeoLibPointF[]> source, bool useMidPoint)
    {
        int sourceCount = source.Count;

        List<GeoLibPointF> sortPoints = new();

        for (int i = 0; i < sourceCount; i++)
        {
            switch (useMidPoint)
            {
                case true:
                    sortPoints.Add(midPoint(source[i]));
                    break;
                default:
                    sortPoints.Add(new GeoLibPointF(source[i][0]));
                    break;
            }

            sortPoints[i].tag = i; // track our original poly for this midpoint through the re-order
        }

        sortPoints = sortPoints.OrderBy(p => p.X).ThenBy(p => p.Y).ToList();

        List<GeoLibPointF[]> ret = new();

        for (int i = 0; i < sourceCount; i++)
        {
            ret.Add(source[sortPoints[i].tag].ToArray());
        }

        return ret;
    }

    public static List<GeoLibPoint[]> xySequence(List<GeoLibPoint[]> source, bool useMidPoint = false)
    {
        return pXYSequence(source, useMidPoint);
    }

    private static List<GeoLibPoint[]> pXYSequence(List<GeoLibPoint[]> source, bool useMidPoint)
    {
        int sourceCount = source.Count;

        List<GeoLibPointF> sortPoints = new();

        for (int i = 0; i < sourceCount; i++)
        {
            switch (useMidPoint)
            {
                case true:
                    sortPoints.Add(midPoint(source[i]));
                    break;
                default:
                    sortPoints.Add(new GeoLibPointF(source[i][0]));
                    break;
            }

            sortPoints[i].tag = i; // track our original poly for this midpoint through the re-order
        }

        sortPoints = sortPoints.OrderBy(p => p.X).ThenBy(p => p.Y).ToList();

        List<GeoLibPoint[]> ret = new();

        for (int i = 0; i < sourceCount; i++)
        {
            ret.Add(source[sortPoints[i].tag].ToArray());
        }

        return ret;
    }

    public static List<GeoLibPointF[]> clockwiseSequence(List<GeoLibPointF[]> source)
    {
        return pClockwiseSequence(source);
    }

    private static List<GeoLibPointF[]> pClockwiseSequence(List<GeoLibPointF[]> source)
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
                if (!(Math.Abs(reorderedMidPoints[i].X - midPoints[j].X) <= double.Epsilon) ||
                    !(Math.Abs(reorderedMidPoints[i].Y - midPoints[j].Y) <= double.Epsilon))
                {
                    continue;
                }

                newOrder[i] = j;
                break;
            }
        }

        List<GeoLibPointF[]> ret = new();
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

    private static List<GeoLibPoint[]> pClockwiseSequence(List<GeoLibPoint[]> source)
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
                if (!(Math.Abs(reorderedMidPoints[i].X - midPoints[j].X) <= double.Epsilon) ||
                    !(Math.Abs(reorderedMidPoints[i].Y - midPoints[j].Y) <= double.Epsilon))
                {
                    continue;
                }

                newOrder[i] = j;
                break;
            }
        }

        List<GeoLibPoint[]> ret = new();
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

    private static GeoLibPoint[] pClockwise(GeoLibPoint[] points)
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

    private static GeoLibPointF[] pClockwise(GeoLibPointF[] iPoints)
    {
        if (!pIsClockwise(iPoints))
        {
            Array.Reverse(iPoints);
        }

        return iPoints;
    }

    public static Paths clockwise(Paths source)
    {
        return source.Select(t => pClockwise(t.ToList())).ToList();
    }

    public static Path clockwise(Path iPoints)
    {
        return pClockwise(iPoints);
    }

    private static Path pClockwise(Path iPoints)
    {
        if (!pIsClockwise(iPoints))
        {
            Array.Reverse(iPoints.ToArray());
        }

        return iPoints;
    }
}