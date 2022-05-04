using geoLib;
using LibTessDotNet.Double;
using System;
using System.Collections.Generic;
using System.Linq;
using Clipper2Lib;
using utility;

namespace geoWrangler;

using Path = List<Point64>;
using Paths = List<List<Point64>>;

public static partial class GeoWrangler
{
    public enum outerCutterIndex { outer, cutter }

    public static List<GeoLibPoint[]> clockwiseAndReorder(List<GeoLibPoint[]> iPoints)
    {
        return pClockwiseAndReorder(iPoints);
    }

    private static List<GeoLibPoint[]> pClockwiseAndReorder(List<GeoLibPoint[]> iPoints)
    {
        return iPoints.Select(t => pClockwiseAndReorder(t)).ToList();
    }

    public static GeoLibPoint[] clockwiseAndReorder(GeoLibPoint[] iPoints)
    {
        return pClockwiseAndReorder(iPoints);
    }

    private static GeoLibPoint[] pClockwiseAndReorder(GeoLibPoint[] iPoints)
    {
        iPoints = pClockwise(iPoints);
        iPoints = pReorder(iPoints);
        return iPoints;
    }

    private static GeoLibPoint[] pReorder(GeoLibPoint[] iPoints)
    {
        int minX_index = MinX(iPoints);
        long minX = iPoints[minX_index].X;
        // This will reorder the point index so that the 0-indexed point is at the minimum X value, and, in the case of multiple points at min X, at the lowest Y of all of those.
        List<int> minXPoints = new();
        for (int pt = 0; pt < iPoints.Length; pt++)
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
            List<GeoLibPoint> tempList = new();
            // Now to start the re-indexing.
            for (int pt = reIndexStart; pt < iPoints.Length; pt++)
            {
                tempList.Add(new GeoLibPoint(iPoints[pt].X, iPoints[pt].Y));
            }
            // Ensure we close the shape by hitting the reIndexStart point again, since we will possibly have pushed it to the beginning of the shape.
            for (int pt = 0; pt <= reIndexStart; pt++)
            {
                tempList.Add(new GeoLibPoint(iPoints[pt].X, iPoints[pt].Y));
            }

            iPoints = tempList.ToArray();
        }

        return iPoints;
    }

    public static Paths clockwiseAndReorder(Paths iPoints)
    {
        return pClockwiseAndReorder(iPoints);
    }

    private static Paths pClockwiseAndReorder(Paths iPoints)
    {
        Paths retPaths = new();
        foreach (Path t in iPoints.Select(t1 => pClockwiseAndReorder(t1)))
        {
            t.Reverse(); // Getting a reversed path from the above, not sure why.
            retPaths.Add(pClose(t));
        }

        return retPaths;
    }

    public static Path clockwiseAndReorder(Path iPoints)
    {
        return pClockwiseAndReorder(iPoints);
    }

    private static Path pClockwiseAndReorder(Path iPoints)
    {
        iPoints = pClockwise(iPoints);
        iPoints = pReorder(iPoints);

        return iPoints;
    }

    public static Paths reOrder(Paths iPoints)
    {
        return pReorder(iPoints);
    }

    private static Paths pReorder(Paths iPoints)
    {
        return iPoints.Select(t => pReorder(t)).ToList();
    }

    public static Path reOrder(Path iPoints)
    {
        return pReorder(iPoints);
    }

    private static Path pReorder(Path iPoints)
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
            Path tempList = new();
            // Now to start the re-indexing.
            for (int pt = reIndexStart; pt < iPoints.Count; pt++)
            {
                tempList.Add(new Point64(iPoints[pt].X, iPoints[pt].Y, iPoints[pt].Z));
            }
            // Ensure we close the shape by hitting the reIndexStart point again, since we will possibly have pushed it to the beginning of the shape.
            for (int pt = 0; pt <= reIndexStart; pt++)
            {
                tempList.Add(new Point64(iPoints[pt].X, iPoints[pt].Y, iPoints[pt].Z));
            }

            iPoints = tempList.ToList();
        }

        return iPoints;
    }

    public static List<GeoLibPointF[]> clockwiseAndReorder(List<GeoLibPointF[]> iPoints)
    {
        return pClockwiseAndReorder(iPoints);
    }

    private static List<GeoLibPointF[]> pClockwiseAndReorder(List<GeoLibPointF[]> iPoints)
    {
        List<GeoLibPointF[]> ret = new();
        foreach (GeoLibPointF[] t in iPoints)
        {
            ret.Add(pClockwiseAndReorder(t));
        }
        return ret;
    }

    public static GeoLibPointF[] clockwiseAndReorder(GeoLibPointF[] iPoints)
    {
        return pClockwiseAndReorder(iPoints);
    }

    private static GeoLibPointF[] pClockwiseAndReorder(GeoLibPointF[] iPoints)
    {
        iPoints = pClockwise(iPoints);
        iPoints = pReorder(iPoints);
        return iPoints;
    }

    private static GeoLibPointF[] pReorder(GeoLibPointF[] iPoints)
    {
        int minX_index = MinX(iPoints);
        double minX = iPoints[minX_index].X;
        // This will reorder the point index so that the 0-indexed point is at the minimum X value, and, in the case of multiple points at min X, at the lowest Y of all of those.
        List<int> minXPoints = new();
        for (int pt = 0; pt < iPoints.Length; pt++)
        {
            switch (Math.Abs(iPoints[pt].X - minX))
            {
                case <= double.Epsilon:
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
            return iPoints;
        }

        {
            List<GeoLibPointF> tempList = new();
            // Now to start the re-indexing.
            for (int pt = reIndexStart; pt < iPoints.Length; pt++)
            {
                switch (tempList.Count)
                {
                    // Avoid adding duplicate vertices
                    case > 1 when Math.Abs(tempList[^1].X - iPoints[pt].X) <= double.Epsilon && Math.Abs(tempList[^1].Y - iPoints[pt].Y) <= double.Epsilon:
                        continue;
                    default:
                        tempList.Add(new GeoLibPointF(iPoints[pt].X, iPoints[pt].Y));
                        break;
                }
            }
            // Ensure we close the shape by hitting the reIndexStart point again, since we will possibly have pushed it to the beginning of the shape.
            for (int pt = 0; pt <= reIndexStart; pt++)
            {
                switch (tempList.Count)
                {
                    // Avoid adding duplicate vertices
                    case > 1 when Math.Abs(tempList[^1].X - iPoints[pt].X) <= double.Epsilon && Math.Abs(tempList[^1].Y - iPoints[pt].Y) <= double.Epsilon:
                        continue;
                    default:
                        tempList.Add(new GeoLibPointF(iPoints[pt].X, iPoints[pt].Y));
                        break;
                }
            }

            iPoints = tempList.ToArray();
        }

        return iPoints;
    }

    public static List<GeoLibPoint[]> simplify(List<GeoLibPoint[]> source)
    {
        return pSimplify(source);
    }

    private static List<GeoLibPoint[]> pSimplify(List<GeoLibPoint[]> source)
    {
        return source.Select(t => pSimplify(t)).ToList();
    }
    public static GeoLibPoint[] simplify(GeoLibPoint[] iPoints)
    {
        return pSimplify(iPoints);
    }

    private static GeoLibPoint[] pSimplify(GeoLibPoint[] iPoints)
    {
        List<Point64> iPoly = pathFromPoint(iPoints, 1);
        Clipper c = new();
        c.PreserveCollinear = false;
        c.AddSubject(iPoly);
        List<List<Point64>> oPoly = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, oPoly);

        oPoly = pReorder(oPoly);

        GeoLibPoint[] working = pointFromPath(oPoly[0], 1);

        return working;
    }

    public static Paths[] getOutersAndCutters(Paths source)
    {
        return pGetOutersAndCutters(source);
    }

    private static Paths[] pGetOutersAndCutters(Paths source)
    {
        Paths[] ret = new Paths[2];
        // Find cutters and outers.
        Paths outers = new();
        Paths cutters = new();
        foreach (Path t in source)
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

    public static Paths stripColinear(Paths source, double angularTolerance = 0.0f)
    {
        return pStripColinear(source, angularTolerance);
    }
    
    private static Paths pStripColinear(Paths source, double angularTolerance = 0.0f)
    {
        return source.Select(t => pStripColinear(t, angularTolerance)).ToList();
    }

    public static Path stripColinear(Path source, double angularTolerance = 0.0f)
    {
        return pStripColinear(source, angularTolerance);
    }

    private static Path pStripColinear(Path source, double angularTolerance = 0.0f)
    {
        switch (source.Count)
        {
            case < 3:
                return source;
        }

        Path ret = new();

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

            if (pt == 0 || Math.Abs(theta - 180) > angularTolerance)
            {
                ret.Add(new Point64(source[pt]));
            }
        }
        return ret;
    }

    public static List<GeoLibPoint[]> stripColinear(List<GeoLibPoint[]> source, double angularTolerance = 0.0f)
    {
        return pStripColinear(source, angularTolerance);
    }

    private static List<GeoLibPoint[]> pStripColinear(List<GeoLibPoint[]> source, double angularTolerance = 0.0f)
    {
        List<GeoLibPoint[]> ret = new();
        foreach (GeoLibPoint[] t in source)
        {
            ret.Add(pStripColinear(t, angularTolerance));
        }

        return ret;
    }

    public static GeoLibPoint[] stripColinear(GeoLibPoint[] source, double angularTolerance = 0.0f)
    {
        return pStripColinear(source, angularTolerance);
    }

    private static GeoLibPoint[] pStripColinear(GeoLibPoint[] source, double angularTolerance = 0.0f)
    {
        switch (source.Length)
        {
            case < 3:
                return source;
        }

        List<GeoLibPoint> ret = new();

        for (int pt = 0; pt < source.Length; pt++)
        {
            GeoLibPoint interSection_A, interSection_B, interSection_C;
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
                    if (pt == source.Length - 1) // last point in the list
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
            if (pt != 0 && pt != source.Length - 1)
            {
                if (Math.Abs(theta - 180) < angularTolerance)
                {
                    addPoint = false;
                }
            }

            switch (addPoint)
            {
                case true:
                    ret.Add(new GeoLibPoint(source[pt]));
                    break;
            }
        }
        return ret.ToArray();
    }

    public static List<GeoLibPoint> stripColinear(List<GeoLibPoint> source, double angularTolerance = 0.0f)
    {
        return pStripColinear(source, angularTolerance);
    }

    private static List<GeoLibPoint> pStripColinear(List<GeoLibPoint> source, double angularTolerance = 0.0f)
    {
        switch (source.Count)
        {
            case < 3:
                return source;
        }

        List<GeoLibPoint> ret = new();

        for (int pt = 0; pt < source.Count; pt++)
        {
            GeoLibPoint interSection_A, interSection_B, interSection_C;
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
                if (Math.Abs(theta - 180) < angularTolerance)
                {
                    addPoint = false;
                }
            }

            switch (addPoint)
            {
                case true:
                    ret.Add(new GeoLibPoint(source[pt]));
                    break;
            }
        }
        return ret;
    }


    public static GeoLibPoint[]  removeDuplicates(GeoLibPoint[] source)
    {
        return pRemoveDuplicates(source).ToArray();
    }

    public static List<GeoLibPoint> pRemoveDuplicates(GeoLibPoint[] source)
    {
        return pRemoveDuplicates(source.ToList());
    }

    private static List<GeoLibPoint> pRemoveDuplicates(List<GeoLibPoint> source)
    {
        List<GeoLibPoint> ret = new();
        switch (source.Count)
        {
            case > 0:
            {
                ret.Add(new GeoLibPoint(source[0]));
                int retIndex = 1;
                for (int i = 1; i < source.Count - 1; i++)
                {
                    if (source[i].X == ret[retIndex - 1].X && source[i].Y == ret[retIndex - 1].Y)
                    {
                        continue;
                    }

                    ret.Add(new GeoLibPoint(source[i]));
                    retIndex++;
                }

                break;
            }
        }

        return ret;
    }


    public static GeoLibPointF[] removeDuplicates(GeoLibPointF[] source)
    {
        return pRemoveDuplicates(source).ToArray();
    }

    public static List<GeoLibPointF> pRemoveDuplicates(GeoLibPointF[] source)
    {
        return pRemoveDuplicates(source.ToList());
    }

    private static List<GeoLibPointF> pRemoveDuplicates(List<GeoLibPointF> source)
    {
        List<GeoLibPointF> ret = new();
        switch (source.Count)
        {
            case > 0:
            {
                ret.Add(new GeoLibPointF(source[0]));
                int retIndex = 1;
                for (int i = 1; i < source.Count - 1; i++)
                {
                    if (!(Math.Abs(source[i].X - ret[retIndex - 1].X) > double.Epsilon) &&
                        !(Math.Abs(source[i].Y - ret[retIndex - 1].Y) > double.Epsilon))
                    {
                        continue;
                    }

                    ret.Add(new GeoLibPointF(source[i]));
                    retIndex++;
                }

                break;
            }
        }

        return ret;
    }

    public static List<GeoLibPoint> stripTerminators(List<GeoLibPoint> source, bool keepLast)
    {
        return pStripTerminators(source, keepLast);
    }

    public static GeoLibPoint[] stripTerminators(GeoLibPoint[] source, bool keepLast)
    {
        return pStripTerminators(source, keepLast);
    }

    private static GeoLibPoint[] pStripTerminators(GeoLibPoint[] source, bool keepLast)
    {
        return pStripTerminators(source.ToList(), keepLast).ToArray();
    }

    private static List<GeoLibPoint> pStripTerminators(List<GeoLibPoint> source, bool keepLast)
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

    public static List<GeoLibPointF[]> stripColinear(List<GeoLibPointF[]> source, double angularTolerance = 0.0f)
    {
        return pStripColinear(source, angularTolerance);
    }

    private static List<GeoLibPointF[]> pStripColinear(List<GeoLibPointF[]> source, double angularTolerance = 0.0f)
    {
        return source.Select(t => pStripColinear(t, angularTolerance)).ToList();
    }

    public static GeoLibPointF[] stripColinear(GeoLibPointF[] source, double angularTolerance = 0.0f)
    {
        return pStripColinear(source, angularTolerance);
    }

    private static GeoLibPointF[] pStripColinear(GeoLibPointF[] source, double angularTolerance = 0.0f)
    {
        switch (source.Length)
        {
            case < 3:
                return source;
        }

        List<GeoLibPointF> ret = new();

        for (int pt = 0; pt < source.Length; pt++)
        {
            GeoLibPointF interSection_A, interSection_B, interSection_C;
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
                    if (pt == source.Length - 1) // last point in the list
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
            if (pt != 0 && pt != source.Length - 1)
            {
                if (Math.Abs(theta - 180) < angularTolerance)
                {
                    addPoint = false;
                }
            }

            switch (addPoint)
            {
                case true:
                    ret.Add(new GeoLibPointF(source[pt]));
                    break;
            }
        }
        return ret.ToArray();
    }

    public static List<GeoLibPointF> stripColinear(List<GeoLibPointF> source, double angularTolerance = 0.0f)
    {
        return pStripColinear(source, angularTolerance);
    }

    private static List<GeoLibPointF> pStripColinear(List<GeoLibPointF> source, double angularTolerance = 0.0f)
    {
        switch (source.Count)
        {
            case < 3:
                return source;
        }

        List<GeoLibPointF> ret = new();

        for (int pt = 0; pt < source.Count; pt++)
        {
            GeoLibPointF interSection_A, interSection_B, interSection_C;
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
                if (Math.Abs(theta - 180) < angularTolerance)
                {
                    addPoint = false;
                }
            }

            switch (addPoint)
            {
                case true:
                    ret.Add(new GeoLibPointF(source[pt]));
                    break;
            }
        }
        return ret;
    }

    public static List<GeoLibPointF> stripTerminators(List<GeoLibPointF> source, bool keepLast)
    {
        return pStripTerminators(source, keepLast);
    }

    public static GeoLibPointF[] stripTerminators(GeoLibPointF[] source, bool keepLast)
    {
        return pStripTerminators(source, keepLast);
    }

    private static GeoLibPointF[] pStripTerminators(GeoLibPointF[] source, bool keepLast)
    {
        return pStripTerminators(source.ToList(), keepLast).ToArray();
    }

    private static List<GeoLibPointF> pStripTerminators(List<GeoLibPointF> source, bool keepLast)
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

    public static Paths stripTerminators(Paths source, bool keepLast)
    {
        return pStripTerminators(source, keepLast);
    }

    private static Paths pStripTerminators(Paths source, bool keepLast)
    {
        Paths ret = new();
        foreach (Path t in source)
        {
            ret.Add(pStripTerminators(t, keepLast));
        }

        return ret;
    }

    public static Path stripTerminators(Path source, bool keepLast)
    {
        return pStripTerminators(source, keepLast);
    }

    private static Path pStripTerminators(Path source, bool keepLast)
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

    public static Paths close(Paths source)
    {
        return pClose(source);
    }

    private static Paths pClose(Paths source)
    {
        Paths ret = new();
        foreach (Path t in source)
        {
            ret.Add(pClose(t));
        }
        return ret;
    }
    public static Path close(Path source)
    {
        return pClose(source);
    }

    private static Path pClose(Path source)
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

    public static List<GeoLibPoint> close(List<GeoLibPoint> source)
    {
        return pClose(source);
    }

    private static List<GeoLibPoint> pClose(List<GeoLibPoint> source)
    {
        switch (source.Count)
        {
            case < 1:
                return source;
        }
        if (source[0].X != source[^1].X || source[0].Y != source[^1].Y)
        {
            source.Add(new GeoLibPoint(source[0]));
        }
        return source;
    }

    public static List<GeoLibPoint[]> close(List<GeoLibPoint[]> source)
    {
        return pClose(source);
    }

    static List<GeoLibPoint[]> pClose(List<GeoLibPoint[]> source)
    {
        List<GeoLibPoint[]> ret = new();
        foreach (GeoLibPoint[] t in source)
        {
            ret.Add(pClose(t));
        }

        return ret;
    }

    public static GeoLibPoint[] close(GeoLibPoint[] source)
    {
        return pClose(source);
    }

    private static GeoLibPoint[] pClose(GeoLibPoint[] source)
    {
        switch (source.Length)
        {
            case < 1:
                return source;
        }
        List<GeoLibPoint> n = source.ToList();
        n = close(n);
        return n.ToArray();
    }

    public static List<GeoLibPointF> close(List<GeoLibPointF> source)
    {
        return pClose(source);
    }

    private static List<GeoLibPointF> pClose(List<GeoLibPointF> source)
    {
        switch (source.Count)
        {
            case < 1:
                return source;
        }
        if (Math.Abs(source[0].X - source[^1].X) > double.Epsilon || Math.Abs(source[0].Y - source[^1].Y) > double.Epsilon)
        {
            source.Add(new GeoLibPointF(source[0]));
        }
        return source;
    }
    
    public static List<GeoLibPointF[]> close(List<GeoLibPointF[]> source)
    {
        return pClose(source);
    }

    private static List<GeoLibPointF[]> pClose(List<GeoLibPointF[]> source)
    {
        List<GeoLibPointF[]> ret = new();
        foreach (GeoLibPointF[] p in source)
        {
            ret.Add(pClose(p));
        }

        return ret;
    }

    public static GeoLibPointF[] close(GeoLibPointF[] source)
    {
        return pClose(source);
    }

    private static GeoLibPointF[] pClose(GeoLibPointF[] source)
    {
        switch (source.Length)
        {
            case < 1:
                return source;
        }
        List<GeoLibPointF> n = source.ToList();
        n = close(n);
        return n.ToArray();
    }

    public static Paths fromSoup(Paths source)
    {
        return pFromSoup(source, WindingRule.NonZero);
    }

    private static Paths pFromSoup(Paths source, WindingRule wr)
    {
        Paths outers = new();
        Paths cutters = new();

        foreach (Path t in source)
        {
            if (ClipperFunc.IsClockwise(t) == ClipperFunc.IsClockwise(source[0]))
            {
                outers.Add(new Path(t));
            }
            else
            {
                cutters.Add(new Path(t));
            }
        }

        // Set up our contours. Since Clipper sets up the subject as the first item, we'll make that clockwise and force the rest to counterclockwise.	
        Tess tess = new();

        foreach (Path t in outers)
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

        foreach (Path t in cutters)
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
            Clipper c = new() {PreserveCollinear = true};
            Paths retPaths = new();

            Paths cPaths = new();
            Paths aPaths = new();

            for (int i = 0; i < tess.ElementCount; i++)
            {
                Path trianglePath = new();
                for (int p = 0; p < polysize; p++)
                {
                    Point64 tmpPt = new((long)tess.Vertices[tess.Elements[i * polysize + p]].Position.X, (long)tess.Vertices[tess.Elements[i * polysize + p]].Position.Y);
                    trianglePath.Add(tmpPt);
                }

                if (ClipperFunc.IsClockwise(trianglePath))
                {
                    cPaths.Add(trianglePath.ToList());
                }
                else
                {
                    aPaths.Add(trianglePath.ToList());
                }
            }

            // Add paths to the clipper.	
            c.AddSubject(cPaths);
            c.AddClip(aPaths);

            c.Execute(ClipType.Union, FillRule.NonZero, retPaths);

            retPaths = pReorder(retPaths);

            retPaths = pClose(retPaths);

            return retPaths;
        }
        catch (Exception)
        {
            return source;
        }
    }

    public static List<GeoLibPointF[]> clean_and_flatten(List<GeoLibPointF[]> source, long scaling, double customSizing = 0, double extension = 0)
    {
        return pClean_and_flatten(source, scaling, customSizing, extension);
    }

    private static List<GeoLibPointF[]> pClean_and_flatten(List<GeoLibPointF[]> source, long scaling, double customSizing = 0, double extension = 0)
    {
        Paths sourcePaths = pPathsFromPointFs(source, scaling);
        Clipper c = new();
        c.AddSubject(sourcePaths);
        Paths solution = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, solution);

        solution = pReorder(solution);

        Paths keyHoled = pMakeKeyHole(solution, customSizing: customSizing, extension: extension);

        return pPointFsFromPaths(pClockwiseAndReorder(keyHoled), scaling);
    }

    public static Paths removeDuplicatePaths(Paths source)
    {
        return pRemoveDuplicatePaths(source);
    }
    
    private static Paths pRemoveDuplicatePaths(Paths source)
    {
        Paths ret = new();
        List<string> polyHashCodes = new();

        foreach (Path p in source)
        {
            string polyHash = Utils.GetMD5Hash(p);
            if (polyHashCodes.IndexOf(polyHash) != -1)
            {
                continue;
            }
            polyHashCodes.Add(polyHash);
            ret.Add(p.ToList());
        }

        return ret;
    }
    
}