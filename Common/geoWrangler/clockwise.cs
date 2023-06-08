using Clipper2Lib;
using System;
using System.Linq;

namespace geoWrangler;

public static partial class GeoWrangler
{
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
                if (!(Math.Abs(reorderedMidPoints[i].x - midPoints[j].x) <= Constants.tolerance) ||
                    !(Math.Abs(reorderedMidPoints[i].y - midPoints[j].y) <= Constants.tolerance))
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
    
    public static PathD clockwise(PathD iPoints)
    {
        return pClockwise(iPoints);
    }

    private static PathD pClockwise(PathD iPoints)
    {
        PathD ret = new(iPoints);
        if (!pIsClockwise(ret))
        {
            ret.Reverse();
        }

        return ret;
    }

    public static PathsD clockwise(PathsD source)
    {
        return new (source.Select(t => pClockwise(new PathD (t))));
    }

    public static PathsD clockwiseAndReorderXY(PathsD iPoints)
    {
        return pClockwiseAndReorderXY(iPoints);
    }

    private static PathsD pClockwiseAndReorderXY(PathsD iPoints)
    {
        PathsD ret = new();
        foreach (PathD t in iPoints)
        {
            ret.Add(pClockwiseAndReorderXY(t));
        }
        return ret;
    }

    public static PathD clockwiseAndReorderXY(PathD iPoints)
    {
        return pClockwiseAndReorderXY(iPoints);
    }

    private static PathD pClockwiseAndReorderXY(PathD iPoints)
    {
        PathD ret = pClockwise(iPoints);
        ret = pReorderXY(ret);
        return ret;
    }
    
    public static PathsD clockwiseAndReorderYX(PathsD iPoints)
    {
        return pClockwiseAndReorderYX(iPoints);
    }

    private static PathsD pClockwiseAndReorderYX(PathsD iPoints)
    {
        PathsD ret = new();
        foreach (PathD t in iPoints)
        {
            ret.Add(pClockwiseAndReorderYX(t));
        }
        return ret;
    }
    
    public static PathD clockwiseAndReorderYX(PathD iPoints)
    {
        return pClockwiseAndReorderYX(iPoints);
    }

    private static PathD pClockwiseAndReorderYX(PathD iPoints)
    {
        PathD ret = pClockwise(iPoints);
        ret = pReorderYX(ret);
        return ret;
    }
}