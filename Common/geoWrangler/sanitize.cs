using geoLib;
using LibTessDotNet.Double;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using Clipper2Lib;
using utility;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public enum outerCutterIndex { outer, cutter }

    public static Paths64 clockwiseAndReorderXY(Paths64 iPoints)
    {
        return pClockwiseAndReorderXY(iPoints);
    }

    private static Paths64 pClockwiseAndReorderXY(Paths64 iPoints)
    {
        int sLength = iPoints.Count;
        Paths64 ret = new (sLength);
#if !GWSINGLETHREADED
        Parallel.For(0, sLength, pt =>
#else
            for (int pt = 0; pt < sLength; pt++)
#endif
            {
                ret[pt] = pClockwiseAndReorderXY(iPoints[pt]);
            }
#if !GWSINGLETHREADED
        );
#endif
        return ret;
    }

    public static Path64 clockwiseAndReorderXY(Path64 iPoints)
    {
        return pClockwiseAndReorderXY(iPoints);
    }

    private static Path64 pClockwiseAndReorderXY(Path64 iPoints)
    {
        iPoints = pClockwise(iPoints);
        iPoints = pReorderXY(iPoints);
        return iPoints;
    }

    public static Paths64 reorderXY(Paths64 iPoints)
    {
        return pReorderXY(iPoints);
    }

    private static Paths64 pReorderXY(Paths64 iPoints)
    {
        int sLength = iPoints.Count;
        Paths64 ret = new (sLength);
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

    private static Path64 pReorderXY(Path64 iPoints)
    {
        int minX_index = MinX(iPoints);
        long minX = iPoints[minX_index].X;
        // This will reorder the point index so that the 0-indexed point is at the minimum X value, and, in the case of multiple points at min X, at the lowest Y of all of those.
        List<int> minXPoints = new();
        for (int pt = 0; pt < iPoints.Count; pt++)
        {
            if (iPoints[pt].X == minX)
            {
                minXPoints.Add(pt);
            }
        }
        // Now we need to query our minXPoints to find the point with the lowest Y value.
        long minY = iPoints[minXPoints[0]].Y;
        int reIndexStart = minXPoints[0];
        for (int index = 1; index < minXPoints.Count; index++)
        {
            if (iPoints[minXPoints[index]].Y >= minY)
            {
                continue;
            }

            minY = iPoints[minXPoints[index]].Y;
            reIndexStart = minXPoints[index];
        }

        if (reIndexStart == 0)
        {
            return iPoints;
        }

        {
            Path64 tempList = new();
            // Now to start the re-indexing.
            for (int pt = reIndexStart; pt < iPoints.Count; pt++)
            {
                tempList.Add(new (iPoints[pt]));
            }
            // Ensure we close the shape by hitting the reIndexStart point again, since we will possibly have pushed it to the beginning of the shape.
            for (int pt = 0; pt <= reIndexStart; pt++)
            {
                tempList.Add(new (iPoints[pt]));
            }

            iPoints = new (tempList);
        }

        return iPoints;
    }


    public static Paths64 clockwiseAndReorderYX(Paths64 iPoints)
    {
        return pClockwiseAndReorderYX(iPoints);
    }

    private static Paths64 pClockwiseAndReorderYX(Paths64 iPoints)
    {
        int sLength = iPoints.Count;
        Paths64 ret = new (sLength);
#if !GWSINGLETHREADED
        Parallel.For(0, sLength, pt =>
#else
            for (int pt = 0; pt < sLength; pt++)
#endif
            {
                ret[pt] = pClockwiseAndReorderYX(iPoints[pt]);
            }
#if !GWSINGLETHREADED
        );
#endif
        return ret;
    }

    public static Path64 clockwiseAndReorderYX(Path64 iPoints)
    {
        return pClockwiseAndReorderYX(iPoints);
    }

    private static Path64 pClockwiseAndReorderYX(Path64 iPoints)
    {
        iPoints = pClockwise(iPoints);
        iPoints = pReorderYX(iPoints);
        return iPoints;
    }

    public static Path64 reOrderYX(Path64 iPoints)
    {
        return pReorderYX(iPoints);
    }

    private static Path64 pReorderYX(Path64 iPoints)
    {
        int minY_index = MinY(iPoints);
        long minY = iPoints[minY_index].Y;
        // This will reorder the point index so that the 0-indexed point is at the minimum Y value, and, in the case of multiple points at min Y, at the lowest X of all of those.
        List<int> minYPoints = new();
        for (int pt = 0; pt < iPoints.Count; pt++)
        {
            if (iPoints[pt].Y == minY)
            {
                minYPoints.Add(pt);
            }
        }
        // Now we need to query our minYPoints to find the point with the lowest X value.
        long minX = iPoints[minYPoints[0]].X;
        int reIndexStart = minYPoints[0];
        for (int index = 1; index < minYPoints.Count; index++)
        {
            if (iPoints[minYPoints[index]].X >= minX)
            {
                continue;
            }

            minX = iPoints[minYPoints[index]].X;
            reIndexStart = minYPoints[index];
        }

        if (reIndexStart == 0)
        {
            return iPoints;
        }

        Path64 tempList = new();
        // Now to start the re-indexing.
        for (int pt = reIndexStart; pt < iPoints.Count; pt++)
        {
            tempList.Add(new (iPoints[pt]));
        }
        // Ensure we close the shape by hitting the reIndexStart point again, since we will possibly have pushed it to the beginning of the shape.
        for (int pt = 0; pt <= reIndexStart; pt++)
        {
            tempList.Add(new (iPoints[pt]));
        }

        iPoints = new(tempList);

        return iPoints;
    }
    
    public static Paths64 reOrderXY(Paths64 iPoints)
    {
        return pReorderXY(iPoints);
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
        Paths64 ret = new (sLength);
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
        iPoints = pClockwise(iPoints);
        iPoints = pReorderXY(iPoints);
        return iPoints;
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
                case <= double.Epsilon:
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
                    case > 1 when Math.Abs(tempList[^1].x - iPoints[pt].x) <= double.Epsilon && Math.Abs(tempList[^1].y - iPoints[pt].y) <= double.Epsilon:
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
                    case > 1 when Math.Abs(tempList[^1].x - iPoints[pt].x) <= double.Epsilon && Math.Abs(tempList[^1].y - iPoints[pt].y) <= double.Epsilon:
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
        iPoints = pClockwise(iPoints);
        iPoints = pReorderYX(iPoints);
        return iPoints;
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
                case <= double.Epsilon:
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
            return iPoints;
        }

        PathD tempList = new();
        // Now to start the re-indexing.
        for (int pt = reIndexStart; pt < iPoints.Count; pt++)
        {
            switch (tempList.Count)
            {
                // Avoid adding duplicate vertices
                case > 1 when Math.Abs(tempList[^1].x - iPoints[pt].x) <= double.Epsilon && Math.Abs(tempList[^1].y - iPoints[pt].y) <= double.Epsilon:
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
                case > 1 when Math.Abs(tempList[^1].x - iPoints[pt].x) <= double.Epsilon && Math.Abs(tempList[^1].y - iPoints[pt].y) <= double.Epsilon:
                    continue;
                default:
                    tempList.Add(new (iPoints[pt]));
                    break;
            }
        }

        iPoints = new(tempList);

        return iPoints;
    }
    
    public static Paths64 simplify(Paths64 source)
    {
        return pSimplify(source);
    }

    private static Paths64 pSimplify(Paths64 source)
    {
        int sLength = source.Count;
        Paths64 ret = new (sLength);
#if !GWSINGLETHREADED
        Parallel.For(0, sLength, pt =>
#else
            for (int pt = 0; pt < sLength; pt++)
#endif
            {
                ret[pt] = pSimplify(source[pt]);
            }
#if !GWSINGLETHREADED
        );
#endif
        return ret;
    }
   
    public static Path64 simplify(Path64 iPoints)
    {
        return pSimplify(iPoints);
    }

    private static Path64 pSimplify(Path64 iPoints)
    {
        Path64 iPoly = pathFromPoint(iPoints, 1);
        Clipper64 c = new();
        c.PreserveCollinear = false;
        c.AddSubject(iPoly);
        Paths64 oPoly = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, oPoly);

        oPoly = pReorderXY(oPoly);

        Path64 working = pointFromPath(oPoly[0], 1);

        return working;
    }

    public static Paths64[] getOutersAndCutters(Paths64 source)
    {
        return pGetOutersAndCutters(source);
    }

    private static Paths64[] pGetOutersAndCutters(Paths64 source)
    {
        Paths64[] ret = new Paths64[2];
        // Find cutters and outers.
        Paths64 outers = new();
        Paths64 cutters = new();
        foreach (Path64 t in source)
        {
            if (pIsClockwise(t))
            {
                outers.Add(t);
            }
            else
            {
                cutters.Add(t);
            }
        }
        ret[(int)outerCutterIndex.outer] = outers;
        ret[(int)outerCutterIndex.cutter] = cutters;

        return ret;
    }
    
    public static Path64 stripColinear(Path64 source, double angularTolerance = 0.0f)
    {
        return pStripColinear(source, angularTolerance);
    }

    private static Path64 pStripColinear(Path64 source, double angularTolerance = 0.0f)
    {
        switch (source.Count)
        {
            case < 3:
                return source;
        }

        Path64 ret = new();

        for (int pt = 0; pt < source.Count; pt++)
        {
            Point64 interSection_A, interSection_B, interSection_C;
            switch (pt)
            {
                // Assess angle.
                case 0:
                    interSection_B = source[^1]; // map to last point
                    interSection_C = source[pt];
                    interSection_A = source[pt + 1];
                    break;
                default:
                {
                    if (pt == source.Count - 1) // last point in the list
                    {
                        interSection_B = source[pt - 1];
                        interSection_C = source[pt];
                        interSection_A = source[0]; // map to the first point
                    }
                    else
                    {
                        interSection_B = source[pt - 1];
                        interSection_C = source[pt];
                        interSection_A = source[pt + 1];
                    }

                    break;
                }
            }

            double theta = pAngleBetweenPoints(interSection_A, interSection_B, interSection_C, false);

            bool addPoint = true;
            if (pt != 0 && pt != source.Count - 1)
            {
                if (Math.Abs(theta - 180) <= angularTolerance)
                {
                    addPoint = false;
                }
            }

            switch (addPoint)
            {
                case true:
                    ret.Add(new Point64(source[pt]));
                    break;
            }
            
        }
        return ret;
    }

    public static Paths64 stripColinear(Paths64 source, double angularTolerance = 0.0f)
    {
        return pStripColinear(source, angularTolerance);
    }

    private static Paths64 pStripColinear(Paths64 source, double angularTolerance = 0.0f)
    {
        Paths64 ret = new();
        foreach (Path64 t in source)
        {
            ret.Add(pStripColinear(t, angularTolerance));
        }

        return ret;
    }
    
    public static Path64 removeDuplicates(Path64 source)
    {
        return pRemoveDuplicates(source);
    }

    private static Path64 pRemoveDuplicates(Path64 source)
    {
        Path64 ret = new();
        switch (source.Count)
        {
            case > 0:
            {
                ret.Add(new (source[0]));
                int retIndex = 1;
                for (int i = 1; i < source.Count - 1; i++)
                {
                    if (source[i].X == ret[retIndex - 1].X && source[i].Y == ret[retIndex - 1].Y)
                    {
                        continue;
                    }

                    ret.Add(new (source[i]));
                    retIndex++;
                }

                break;
            }
        }

        return ret;
    }

    public static Path64 removeDuplicates(Path64 source, double threshold = Double.Epsilon)
    {
        return pRemoveDuplicates(source, threshold);
    }

    private static Path64 pRemoveDuplicates(Path64 source, double threshold = Double.Epsilon)
    {
        Path64 ret = new();
        switch (source.Count)
        {
            case > 0:
            {
                ret.Add(new Point64(source[0]));
                int retIndex = 1;
                for (int i = 1; i < source.Count - 1; i++)
                {
                    if (!(Math.Abs(source[i].X - ret[retIndex - 1].X) > threshold) &&
                        !(Math.Abs(source[i].Y - ret[retIndex - 1].Y) > threshold))
                    {
                        continue;
                    }

                    ret.Add(new Point64(source[i]));
                    retIndex++;
                }

                break;
            }
        }

        return ret;
    }
    
    public static PathD removeDuplicates(PathD source, double threshold = Double.Epsilon)
    {
        return pRemoveDuplicates(source, threshold);
    }
    
    private static PathsD pRemoveDuplicates(PathsD source, double threshold = Double.Epsilon)
    {
        PathsD ret = new PathsD();
        foreach (var t in source)
        {
            ret.Add(removeDuplicates(t, threshold));
        }

        return ret;
    }
    
    private static PathD pRemoveDuplicates(PathD source, double threshold = Double.Epsilon)
    {
        PathD ret = new();
        switch (source.Count)
        {
            case > 0:
            {
                ret.Add(new (source[0]));
                int retIndex = 1;
                for (int i = 1; i < source.Count - 1; i++)
                {
                    if (!(Math.Abs(source[i].x - ret[retIndex - 1].x) > threshold) &&
                        !(Math.Abs(source[i].y - ret[retIndex - 1].y) > threshold))
                    {
                        continue;
                    }

                    ret.Add(new (source[i]));
                    retIndex++;
                }

                break;
            }
        }

        return ret;
    }
    
    public static PathsD stripColinear(PathsD source, double angularTolerance = 0.0f)
    {
        return pStripColinear(source, angularTolerance);
    }

    private static PathsD pStripColinear(PathsD source, double angularTolerance = 0.0f)
    {
        int sLength = source.Count;
        PathsD ret = new (sLength);
#if !GWSINGLETHREADED
        Parallel.For(0, sLength, pt =>
#else
            for (int pt = 0; pt < sLength; pt++)
#endif
            {
                ret[pt] = pStripColinear(source[pt], angularTolerance);
            }
#if !GWSINGLETHREADED
        );
#endif

        return ret;
    }

    public static PathD stripColinear(PathD source, double angularTolerance = 0.0f)
    {
        return pStripColinear(source, angularTolerance);
    }

    private static PathD pStripColinear(PathD source, double angularTolerance = 0.0f)
    {
        switch (source.Count)
        {
            case < 3:
                return source;
        }

        PathD ret = new();

        for (int pt = 0; pt < source.Count; pt++)
        {
            PointD interSection_A, interSection_B, interSection_C;
            switch (pt)
            {
                // Assess angle.
                case 0:
                    interSection_B = source[^1]; // map to last point
                    interSection_C = source[pt];
                    interSection_A = source[pt + 1];
                    break;
                default:
                {
                    if (pt == source.Count - 1) // last point in the list
                    {
                        interSection_B = source[pt - 1];
                        interSection_C = source[pt];
                        interSection_A = source[0]; // map to the first point
                    }
                    else
                    {
                        interSection_B = source[pt - 1];
                        interSection_C = source[pt];
                        interSection_A = source[pt + 1];
                    }

                    break;
                }
            }

            double theta = pAngleBetweenPoints(interSection_A, interSection_B, interSection_C, false);

            bool addPoint = true;
            if (pt != 0 && pt != source.Count - 1)
            {
                if (Math.Abs(theta - 180) <= angularTolerance)
                {
                    addPoint = false;
                }
            }

            switch (addPoint)
            {
                case true:
                    ret.Add(new (source[pt]));
                    break;
            }
        }
        return ret;
    }

    public static PathD stripTerminators(PathD source, bool keepLast)
    {
        return pStripTerminators(source, keepLast);
    }
    
    private static PathD pStripTerminators(PathD source, bool keepLast)
    {
        bool firstLast_same = false;
        int pt_Check = source.Count - 1;
        if (distanceBetweenPoints(source[pt_Check], source[0]) < 0.01)
        {
            firstLast_same = true; // remove duplicated points. The shape will be closed later.
        }
        while (firstLast_same)
        {
            source.RemoveAt(pt_Check); // remove duplicated points. The shape will be closed later
            pt_Check--;
            if (distanceBetweenPoints(source[pt_Check], source[0]) > 0.01)
            {
                firstLast_same = false; // stop at the first unmatched point.
            }
        }

        source = keepLast switch
        {
            true => pClose(source),
            _ => source
        };

        return source;
    }

    public static Paths64 stripTerminators(Paths64 source, bool keepLast)
    {
        return pStripTerminators(source, keepLast);
    }

    private static Paths64 pStripTerminators(Paths64 source, bool keepLast)
    {
        Paths64 ret = new();
        foreach (Path64 t in source)
        {
            ret.Add(pStripTerminators(t, keepLast));
        }

        return ret;
    }

    public static Path64 stripTerminators(Path64 source, bool keepLast)
    {
        return pStripTerminators(source, keepLast);
    }

    private static Path64 pStripTerminators(Path64 source, bool keepLast)
    {
        switch (source.Count)
        {
            case <= 1:
                return source;
        }

        bool firstLast_same = false;
        int pt_Check = source.Count - 1;
        if (distanceBetweenPoints(source[pt_Check], source[0]) < 10)
        {
            firstLast_same = true; // remove duplicated points. The shape will be closed later.
        }
        while (firstLast_same)
        {
            source.RemoveAt(pt_Check); // remove duplicated points. The shape will be closed later
            pt_Check--;
            switch (source.Count)
            {
                case < 1:
                    return source;
            }
            if (distanceBetweenPoints(source[pt_Check], source[0]) > 10)
            {
                firstLast_same = false; // stop at the first unmatched point.
            }
        }

        source = keepLast switch
        {
            true => pClose(source),
            _ => source
        };

        return source;
    }

    public static Paths64 close(Paths64 source)
    {
        return pClose(source);
    }

    private static Paths64 pClose(Paths64 source)
    {
        Paths64 ret = new();
        foreach (Path64 t in source)
        {
            ret.Add(pClose(t));
        }
        return ret;
    }
    public static Path64 close(Path64 source)
    {
        return pClose(source);
    }

    private static Path64 pClose(Path64 source)
    {
        switch (source.Count)
        {
            case < 1:
                return source;
        }
        if (source[0].X != source[^1].X || source[0].Y != source[^1].Y)
        {
            source.Add(new Point64(source[0]));
        }
        return source;
    }

    public static PathsD close(PathsD source)
    {
        return pClose(source);
    }
    
    public static PathD close(PathD source)
    {
        return pClose(source);
    }

    private static PathD pClose(PathD source)
    {
        switch (source.Count)
        {
            case < 1:
                return source;
        }
        if (Math.Abs(source[0].x - source[^1].x) > double.Epsilon || Math.Abs(source[0].y - source[^1].y) > double.Epsilon)
        {
            source.Add(new (source[0]));
        }
        return source;
    }
    
    private static PathsD pClose(PathsD source)
    {
        PathsD ret = new();
        foreach (PathD p in source)
        {
            ret.Add(pClose(p));
        }

        return ret;
    }
    
    public static Paths64 fromSoup(Paths64 source)
    {
        return pFromSoup(source, WindingRule.NonZero);
    }

    private static Paths64 pFromSoup(Paths64 source, WindingRule wr)
    {
        Paths64 outers = new();
        Paths64 cutters = new();

        foreach (Path64 t in source)
        {
            if (Clipper.IsPositive(t) == Clipper.IsPositive(source[0]))
            {
                outers.Add(new (t));
            }
            else
            {
                cutters.Add(new (t));
            }
        }

        // Set up our contours. Since Clipper sets up the subject as the first item, we'll make that clockwise and force the rest to counterclockwise.	
        Tess tess = new();

        foreach (Path64 t in outers)
        {
            // Skip tiny fragments. The tessellator has trouble with them.	
            /*
            if (Math.Abs(Clipper.Area(outers[poly])) < 1E-16)
            {
                continue;
            }
            */
            ContourVertex[] contour = new ContourVertex[t.Count];
            for (int i = 0; i < contour.Length; i++)
            {
                contour[i].Position = new Vec3 { X = t[i].X, Y = t[i].Y, Z = 0 };
            }

            tess.AddContour(contour);
        }

        foreach (Path64 t in cutters)
        {
            // Skip tiny fragments. The tessellator has trouble with them.
            /*
            if (Math.Abs(Clipper.Area(cutters[poly])) < 1E-16)
            {
                continue;
            }
            */
            ContourVertex[] contour = new ContourVertex[t.Count];
            for (int i = 0; i < contour.Length; i++)
            {
                contour[i].Position = new Vec3 { X = t[i].X, Y = t[i].Y, Z = 0 };
            }

            tess.AddContour(contour, ContourOrientation.CounterClockwise);
        }

        try
        {
            // Triangulate. This gives us a triangle soup, which may not be entirely helpful for the case where we're punching multiple holes into a mesh. For now, this works, but the limitation will need to be reviewd.	
            // It might be that a PolyTree is needed as the input to the keyhole for these complex cases.... or clockwise paths may need to be evaluated piecewise using all counterclockwise paths for the tessellation.	
            const int polysize = 3;
            tess.Tessellate(wr, ElementType.Polygons, polysize);

            // Iterate triangles and create output geometry. We'll use clipper to simplify the output geometry.	
            Clipper64 c = new() {PreserveCollinear = true};
            Paths64 retPaths = new();

            Paths64 cPaths = new();
            Paths64 aPaths = new();

            for (int i = 0; i < tess.ElementCount; i++)
            {
                Path64 trianglePath = new();
                for (int p = 0; p < polysize; p++)
                {
                    Point64 tmpPt = new((long)tess.Vertices[tess.Elements[i * polysize + p]].Position.X, (long)tess.Vertices[tess.Elements[i * polysize + p]].Position.Y);
                    trianglePath.Add(tmpPt);
                }

                if (Clipper.IsPositive(trianglePath))
                {
                    cPaths.Add(trianglePath);
                }
                else
                {
                    aPaths.Add(trianglePath);
                }
            }

            // Add paths to the clipper.	
            c.AddSubject(cPaths);
            c.AddClip(aPaths);

            c.Execute(ClipType.Union, FillRule.NonZero, retPaths);

            retPaths = pReorderXY(retPaths);

            retPaths = pClose(retPaths);

            return retPaths;
        }
        catch (Exception)
        {
            return source;
        }
    }

    public static PathsD clean_and_flatten(PathsD source, long scaling, double customSizing = 0, double extension = 0)
    {
        return pClean_and_flatten(source, scaling, customSizing, extension);
    }

    private static PathsD pClean_and_flatten(PathsD source, long scaling, double customSizing = 0, double extension = 0)
    {
        Paths64 sourcePaths = pPathsFromPointFs(source, scaling);
        Clipper64 c = new();
        c.AddSubject(sourcePaths);
        Paths64 solution = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, solution);

        solution = pReorderXY(solution);

        Paths64 keyHoled = pMakeKeyHole(solution, reverseEval:false, biDirectionalEval:true, customSizing: customSizing, extension: extension);

        return pPointFsFromPaths(pClockwiseAndReorderXY(keyHoled), scaling);
    }

    public static Paths64 removeDuplicatePaths(Paths64 source)
    {
        return pRemoveDuplicatePaths(source);
    }
    
    private static Paths64 pRemoveDuplicatePaths(Paths64 source)
    {
        Paths64 ret = new();
        List<string> polyHashCodes = new();

        foreach (Path64 p in source)
        {
            string polyHash = Utils.GetMD5Hash(p);
            if (polyHashCodes.IndexOf(polyHash) != -1)
            {
                continue;
            }
            polyHashCodes.Add(polyHash);
            ret.Add(new(p));
        }

        return ret;
    }
    
}