using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Clipper2Lib;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static Paths64 xySequence(Paths64 source, bool useMidPoint = false)
    {
        return pXYSequence(source, useMidPoint);
    }

    private static Paths64 pXYSequence(Paths64 source, bool useMidPoint)
    {
        int sourceCount = source.Count;

        PathD sortPoints = [];

        for (int i = 0; i < sourceCount; i++)
        {
            switch (useMidPoint)
            {
                case true:
                    sortPoints.Add(midPoint(source[i]));
                    break;
                default:
                    sortPoints.Add(new PointD(source[i][0]));
                    break;
            }

            sortPoints[i] = new PointD(sortPoints[i].x, sortPoints[i].y, i); // track our original poly for this midpoint through the re-order
        }

        IOrderedEnumerable<PointD> tmp_sortPoints = sortPoints.OrderBy(p => p.x).ThenBy(p => p.y);
        sortPoints.Clear();
        sortPoints.AddRange(tmp_sortPoints.Select(t => new PointD(t)));

        Paths64 ret = [];

        for (int i = 0; i < sourceCount; i++)
        {
            ret.Add(new Path64(source[(int)sortPoints[i].z]));
        }

        return ret;
    }

    public static Paths64 reOrderXY(Paths64 iPoints)
    {
        return pReorderXY(iPoints);
    }
    
    private static Paths64 pReorderXY(Paths64 iPoints)
    {
        int sLength = iPoints.Count;
        Paths64 ret = Helper.initedPaths64(sLength);
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

    public static Path64 reOrderXY(Path64 iPoints)
    {
        return pReorderXY(iPoints);
    }
    
    public static Paths64 reOrderYX(Paths64 iPoints)
    {
        return pReorderYX(iPoints);
    }

    private static Paths64 pReorderYX(Paths64 iPoints)
    {
        int sLength = iPoints.Count;
        Paths64 ret = Helper.initedPaths64(sLength);
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
    
    private static Path64 pReorderXY(Path64 iPoints)
    {
        int minX_index = MinX(iPoints);
        double minX = iPoints[minX_index].X;
        // This will reorder the point index so that the 0-indexed point is at the minimum X value, and, in the case of multiple points at min X, at the lowest Y of all of those.
        List<int> minXPoints = [];
        for (int pt = 0; pt < iPoints.Count; pt++)
        {
            switch (Math.Abs(iPoints[pt].X - minX))
            {
                case <= Constants.tolerance:
                    minXPoints.Add(pt);
                    break;
            }
        }
        // Now we need to query our minXPoints to find the point with the lowest Y value.
        double minY = iPoints[minXPoints[0]].Y;
        int reIndexStart = minXPoints[0];
        for (int index = 1; index < minXPoints.Count; index++)
        {
            if (!(iPoints[minXPoints[index]].Y < minY))
            {
                continue;
            }

            minY = iPoints[minXPoints[index]].Y;
            reIndexStart = minXPoints[index];
        }

        if (reIndexStart == 0)
        {
            return new Path64(iPoints);
        }

        Path64 tempList = [];
        // Now to start the re-indexing.
        for (int pt = reIndexStart; pt < iPoints.Count; pt++)
        {
            switch (tempList.Count)
            {
                // Avoid adding duplicate vertices
                case > 1 when Math.Abs(tempList[^1].X - iPoints[pt].X) <= Constants.tolerance && Math.Abs(tempList[^1].Y - iPoints[pt].Y) <= Constants.tolerance:
                    continue;
                default:
                    tempList.Add(new Point64(iPoints[pt]));
                    break;
            }
        }
        // Ensure we close the shape by hitting the reIndexStart point again, since we will possibly have pushed it to the beginning of the shape.
        for (int pt = 0; pt <= reIndexStart; pt++)
        {
            switch (tempList.Count)
            {
                // Avoid adding duplicate vertices
                case > 1 when Math.Abs(tempList[^1].X - iPoints[pt].X) <= Constants.tolerance && Math.Abs(tempList[^1].Y - iPoints[pt].Y) <= Constants.tolerance:
                    continue;
                default:
                    tempList.Add(new Point64(iPoints[pt]));
                    break;
            }
        }

        return new Path64(tempList);
    }
    
    public static Path64 reOrderYX(Path64 iPoints)
    {
        return pReorderYX(iPoints);
    }

    private static Path64 pReorderYX(Path64 iPoints)
    {
        int minY_index = MinY(iPoints);
        double minY = iPoints[minY_index].Y;
        // This will reorder the point index so that the 0-indexed point is at the minimum Y value, and, in the case of multiple points at min Y, at the lowest X of all of those.
        List<int> minYPoints = [];
        for (int pt = 0; pt < iPoints.Count; pt++)
        {
            switch (Math.Abs(iPoints[pt].Y - minY))
            {
                case <= Constants.tolerance:
                    minYPoints.Add(pt);
                    break;
            }
        }
        // Now we need to query our minYPoints to find the point with the lowest X value.
        double minX = iPoints[minYPoints[0]].X;
        int reIndexStart = minYPoints[0];
        for (int index = 1; index < minYPoints.Count; index++)
        {
            if (!(iPoints[minYPoints[index]].X < minX))
            {
                continue;
            }

            minX = iPoints[minYPoints[index]].X;
            reIndexStart = minYPoints[index];
        }

        if (reIndexStart == 0)
        {
            return new Path64(iPoints);
        }

        Path64 tempList = [];
        // Now to start the re-indexing.
        for (int pt = reIndexStart; pt < iPoints.Count; pt++)
        {
            switch (tempList.Count)
            {
                // Avoid adding duplicate vertices
                case > 1 when Math.Abs(tempList[^1].X - iPoints[pt].X) <= Constants.tolerance && Math.Abs(tempList[^1].Y - iPoints[pt].Y) <= Constants.tolerance:
                    continue;
                default:
                    tempList.Add(new Point64(iPoints[pt]));
                    break;
            }
        }
        // Ensure we close the shape by hitting the reIndexStart point again, since we will possibly have pushed it to the beginning of the shape.
        for (int pt = 0; pt <= reIndexStart; pt++)
        {
            switch (tempList.Count)
            {
                // Avoid adding duplicate vertices
                case > 1 when Math.Abs(tempList[^1].X - iPoints[pt].X) <= Constants.tolerance && Math.Abs(tempList[^1].Y - iPoints[pt].Y) <= Constants.tolerance:
                    continue;
                default:
                    tempList.Add(new Point64(iPoints[pt]));
                    break;
            }
        }

        return new Path64(tempList);
    }
}