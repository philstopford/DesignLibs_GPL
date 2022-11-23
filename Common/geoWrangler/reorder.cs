using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Clipper2Lib;

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
            sortPoints.Add(new (t));
        }

        PathsD ret = new();

        for (int i = 0; i < sourceCount; i++)
        {
            ret.Add(new(source[(int)sortPoints[i].z]));
        }

        return ret;
    }

    public static PathsD reOrderXY(PathsD iPoints)
    {
        return pReorderXY(iPoints);
    }
    
    private static PathsD pReorderXY(PathsD iPoints)
    {
        int sLength = iPoints.Count;
        PathsD ret = Helper.initedPathsD(sLength);
#if !GWSINGLETHREADED
        Parallel.For(0, sLength, pt =>
#else
            for (int pt = 0; pt < sLength; pt++)
#endif
            {
                ret[pt] = pReorderXY(iPoints[pt]);
            }
#if !GWSINGLETHREADED
        );
#endif
        return ret;
    }

    public static PathD reOrderXY(PathD iPoints)
    {
        return pReorderXY(iPoints);
    }
    
    public static PathsD reOrderYX(PathsD iPoints)
    {
        return pReorderYX(iPoints);
    }

    private static PathsD pReorderYX(PathsD iPoints)
    {
        int sLength = iPoints.Count;
        PathsD ret = Helper.initedPathsD(sLength);
#if !GWSINGLETHREADED
        Parallel.For(0, sLength, pt =>
#else
            for (int pt = 0; pt < sLength; pt++)
#endif
            {
                ret[pt] = pReorderYX(iPoints[pt]);
            }
#if !GWSINGLETHREADED
        );
#endif
        return ret;
    }
    
    private static PathD pReorderXY(PathD iPoints)
    {
        int minX_index = MinX(iPoints);
        double minX = iPoints[minX_index].x;
        // This will reorder the point index so that the 0-indexed point is at the minimum X value, and, in the case of multiple points at min X, at the lowest Y of all of those.
        List<int> minXPoints = new();
        for (int pt = 0; pt < iPoints.Count; pt++)
        {
            switch (Math.Abs(iPoints[pt].x - minX))
            {
                case <= constants.tolerance:
                    minXPoints.Add(pt);
                    break;
            }
        }
        // Now we need to query our minXPoints to find the point with the lowest Y value.
        double minY = iPoints[minXPoints[0]].y;
        int reIndexStart = minXPoints[0];
        for (int index = 1; index < minXPoints.Count; index++)
        {
            if (!(iPoints[minXPoints[index]].y < minY))
            {
                continue;
            }

            minY = iPoints[minXPoints[index]].y;
            reIndexStart = minXPoints[index];
        }

        if (reIndexStart == 0)
        {
            return iPoints;
        }

        {
            PathD tempList = new();
            // Now to start the re-indexing.
            for (int pt = reIndexStart; pt < iPoints.Count; pt++)
            {
                switch (tempList.Count)
                {
                    // Avoid adding duplicate vertices
                    case > 1 when Math.Abs(tempList[^1].x - iPoints[pt].x) <= constants.tolerance && Math.Abs(tempList[^1].y - iPoints[pt].y) <= constants.tolerance:
                        continue;
                    default:
                        tempList.Add(new (iPoints[pt]));
                        break;
                }
            }
            // Ensure we close the shape by hitting the reIndexStart point again, since we will possibly have pushed it to the beginning of the shape.
            for (int pt = 0; pt <= reIndexStart; pt++)
            {
                switch (tempList.Count)
                {
                    // Avoid adding duplicate vertices
                    case > 1 when Math.Abs(tempList[^1].x - iPoints[pt].x) <= constants.tolerance && Math.Abs(tempList[^1].y - iPoints[pt].y) <= constants.tolerance:
                        continue;
                    default:
                        tempList.Add(new (iPoints[pt]));
                        break;
                }
            }

            iPoints = new (tempList);
        }

        return iPoints;
    }
    
    public static PathD reOrderYX(PathD iPoints)
    {
        return pReorderYX(iPoints);
    }

    private static PathD pReorderYX(PathD iPoints)
    {
        int minY_index = MinY(iPoints);
        double minY = iPoints[minY_index].y;
        // This will reorder the point index so that the 0-indexed point is at the minimum Y value, and, in the case of multiple points at min Y, at the lowest X of all of those.
        List<int> minYPoints = new();
        for (int pt = 0; pt < iPoints.Count; pt++)
        {
            switch (Math.Abs(iPoints[pt].y - minY))
            {
                case <= constants.tolerance:
                    minYPoints.Add(pt);
                    break;
            }
        }
        // Now we need to query our minYPoints to find the point with the lowest X value.
        double minX = iPoints[minYPoints[0]].x;
        int reIndexStart = minYPoints[0];
        for (int index = 1; index < minYPoints.Count; index++)
        {
            if (!(iPoints[minYPoints[index]].x < minX))
            {
                continue;
            }

            minX = iPoints[minYPoints[index]].x;
            reIndexStart = minYPoints[index];
        }

        if (reIndexStart == 0)
        {
            return new(iPoints);
        }

        PathD tempList = new();
        // Now to start the re-indexing.
        for (int pt = reIndexStart; pt < iPoints.Count; pt++)
        {
            switch (tempList.Count)
            {
                // Avoid adding duplicate vertices
                case > 1 when Math.Abs(tempList[^1].x - iPoints[pt].x) <= constants.tolerance && Math.Abs(tempList[^1].y - iPoints[pt].y) <= constants.tolerance:
                    continue;
                default:
                    tempList.Add(new (iPoints[pt]));
                    break;
            }
        }
        // Ensure we close the shape by hitting the reIndexStart point again, since we will possibly have pushed it to the beginning of the shape.
        for (int pt = 0; pt <= reIndexStart; pt++)
        {
            switch (tempList.Count)
            {
                // Avoid adding duplicate vertices
                case > 1 when Math.Abs(tempList[^1].x - iPoints[pt].x) <= constants.tolerance && Math.Abs(tempList[^1].y - iPoints[pt].y) <= constants.tolerance:
                    continue;
                default:
                    tempList.Add(new (iPoints[pt]));
                    break;
            }
        }

        return new(tempList);
    }
}