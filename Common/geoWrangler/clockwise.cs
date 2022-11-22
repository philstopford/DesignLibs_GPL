using Clipper2Lib;
using geoLib;
using System;
using System.Collections.Generic;
using System.Linq;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static PathsD xySequence(PathsD source, bool useMidPoint = false)
    {
        return pXYSequence(source, useMidPoint);
    }

    private static PathsD pXYSequence(PathsD source, bool useMidPoint)
    {
        int sourceCount = source.Count;

        PathD sortPoints = new();

        for (int i = 0; i < sourceCount; i++)
        {
            switch (useMidPoint)
            {
                case true:
                    sortPoints.Add(midPoint(source[i]));
                    break;
                default:
                    sortPoints.Add(new (source[i][0]));
                    break;
            }

            sortPoints[i] = new (sortPoints[i].x, sortPoints[i].y, i); // track our original poly for this midpoint through the re-order
        }

        IOrderedEnumerable<PointD> tmp_sortPoints = sortPoints.OrderBy(p => p.x).ThenBy(p => p.y);
        sortPoints.Clear();
        foreach (PointD t in tmp_sortPoints)
        {
            sortPoints.Add(new PointD(t));
        }

        PathsD ret = new();

        for (int i = 0; i < sourceCount; i++)
        {
            ret.Add(source[(int)sortPoints[i].z]);
        }

        return ret;
    }

    public static Paths64 xySequence(Paths64 source, bool useMidPoint = false)
    {
        return pXYSequence(source, useMidPoint);
    }

    private static Paths64 pXYSequence(Paths64 source, bool useMidPoint)
    {
        int sourceCount = source.Count;

        PathD sortPoints = new();

        for (int i = 0; i < sourceCount; i++)
        {
            switch (useMidPoint)
            {
                case true:
                    sortPoints.Add(midPoint(source[i]));
                    break;
                default:
                    sortPoints.Add(new (source[i][0]));
                    break;
            }

            sortPoints[i] = new (sortPoints[i].x, sortPoints[i].y, i); // track our original poly for this midpoint through the re-order
        }

        IOrderedEnumerable<PointD> tmp_sortPoints = sortPoints.OrderBy(p => p.x).ThenBy(p => p.y);
        sortPoints.Clear();
        foreach (PointD t in tmp_sortPoints)
        {
            sortPoints.Add(new PointD(t));
        }

        Paths64 ret = new();

        for (int i = 0; i < sourceCount; i++)
        {
            ret.Add(source[(int)sortPoints[i].z]);
        }

        return ret;
    }

    public static PathsD clockwiseSequence(PathsD source)
    {
        return pClockwiseSequence(source);
    }

    private static PathsD pClockwiseSequence(PathsD source)
    {
        int sourceCount = source.Count;
        PathD midPoints = Helper.initedPathD(sourceCount);

        for (int i = 0; i < sourceCount; i++)
        {
            midPoints[i] = midPoint(source[i]);
        }

        PathD reorderedMidPoints = pClockwise(midPoints);

        // Position in array is new polygon; value is old index.
        int[] newOrder = new int[source.Count];
        for (int i = 0; i < sourceCount; i++)
        {
            for (int j = 0; j < sourceCount; j++)
            {
                if (!(Math.Abs(reorderedMidPoints[i].x - midPoints[j].x) <= double.Epsilon) ||
                    !(Math.Abs(reorderedMidPoints[i].y - midPoints[j].y) <= double.Epsilon))
                {
                    continue;
                }

                newOrder[i] = j;
                break;
            }
        }

        PathsD ret = new();
        for (int i = 0; i < sourceCount; i++)
        {
            ret.Add(new(source[newOrder[i]]));
        }

        return ret;
    }

    public static Paths64 clockwiseSequence(Paths64 source)
    {
        return pClockwiseSequence(source);
    }

    private static Paths64 pClockwiseSequence(Paths64 source)
    {
        int sourceCount = source.Count;
        PathD midPoints = Helper.initedPathD(sourceCount);

        for (int i = 0; i < sourceCount; i++)
        {
            midPoints[i] = midPoint(source[i]);
        }

        PathD reorderedMidPoints = pClockwise(midPoints);

        // Position in array is new polygon; value is old index.
        int[] newOrder = new int[source.Count];
        for (int i = 0; i < sourceCount; i++)
        {
            for (int j = 0; j < sourceCount; j++)
            {
                if (!(Math.Abs(reorderedMidPoints[i].x - midPoints[j].x) <= double.Epsilon) ||
                    !(Math.Abs(reorderedMidPoints[i].y - midPoints[j].y) <= double.Epsilon))
                {
                    continue;
                }

                newOrder[i] = j;
                break;
            }
        }

        Paths64 ret = new();
        for (int i = 0; i < sourceCount; i++)
        {
            ret.Add(new(source[newOrder[i]]));
        }

        return ret;
    }

    public static Path64 clockwise(Path64 points)
    {
        return pClockwise(points);
    }

    private static Path64 pClockwise(Path64 points)
    {
        if (!pIsClockwise(points))
        {
            points.Reverse();
        }

        return points;
    }

    public static PathD clockwise(PathD iPoints)
    {
        return pClockwise(iPoints);
    }

    private static PathD pClockwise(PathD iPoints)
    {
        if (!pIsClockwise(iPoints))
        {
            iPoints.Reverse();
        }

        return iPoints;
    }

    public static Paths64 clockwise(Paths64 source)
    {
        return new (source.Select(t => pClockwise(new Path64 (t))));
    }
}