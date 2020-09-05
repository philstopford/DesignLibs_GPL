using ClipperLib;
using geoLib;
using LibTessDotNet.Double;
using System;
using System.Collections.Generic;
using System.Linq;

namespace geoWrangler
{
    using Path = List<IntPoint>;
    using Paths = List<List<IntPoint>>;

    public static partial class GeoWrangler
    {
        public enum outerCutterIndex { outer, cutter }

        public static List<GeoLibPoint[]> clockwiseAndReorder(List<GeoLibPoint[]> iPoints)
        {
            return pClockwiseAndReorder(iPoints);
        }

        static List<GeoLibPoint[]> pClockwiseAndReorder(List<GeoLibPoint[]> iPoints)
        {
            List<GeoLibPoint[]> ret = new List<GeoLibPoint[]>();
            for (int i = 0; i < iPoints.Count; i++)
            {
                ret.Add(pClockwiseAndReorder(iPoints[i]));
            }
            return ret;
        }

        public static GeoLibPoint[] clockwiseAndReorder(GeoLibPoint[] iPoints)
        {
            return pClockwiseAndReorder(iPoints);
        }

        static GeoLibPoint[] pClockwiseAndReorder(GeoLibPoint[] iPoints)
        {
            iPoints = pClockwise(iPoints);
            iPoints = pReorder(iPoints);
            return iPoints;
        }

        static GeoLibPoint[] pReorder(GeoLibPoint[] iPoints)
        {
            Int32 reIndexStart = 0;
            Int32 minX_index = MinX(iPoints);
            Int64 minX = iPoints[minX_index].X;
            // This will reorder the point index so that the 0-indexed point is at the minimum X value, and, in the case of multiple points at min X, at the lowest Y of all of those.
            List<Int32> minXPoints = new List<Int32>();
            for (Int32 pt = 0; pt < iPoints.Length; pt++)
            {
                if (iPoints[pt].X == minX)
                {
                    minXPoints.Add(pt);
                }
            }
            // Now we need to query our minXPoints to find the point with the lowest Y value.
            Int64 minY = iPoints[minXPoints[0]].Y;
            reIndexStart = minXPoints[0];
            for (Int32 index = 1; index < minXPoints.Count; index++)
            {
                if (iPoints[minXPoints[index]].Y < minY)
                {
                    minY = iPoints[minXPoints[index]].Y;
                    reIndexStart = minXPoints[index];
                }
            }
            if (reIndexStart != 0)
            {
                List<GeoLibPoint> tempList = new List<GeoLibPoint>();
                // Now to start the re-indexing.
                for (Int32 pt = reIndexStart; pt < iPoints.Length; pt++)
                {
                    tempList.Add(new GeoLibPoint(iPoints[pt].X, iPoints[pt].Y));
                }
                // Ensure we close the shape by hitting the reIndexStart point again, since we will possibly have pushed it to the beginning of the shape.
                for (Int32 pt = 0; pt <= reIndexStart; pt++)
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

        static Paths pClockwiseAndReorder(Paths iPoints)
        {
            Paths retPaths = new Paths();
            for (int i = 0; i < iPoints.Count; i++)
            {
                Path t = pClockwiseAndReorder(iPoints[i]);
                t.Reverse(); // Getting a reversed path from the above, not sure why.
                t = pClose(t);
                retPaths.Add(t);

            }

            return retPaths;
        }

        public static Path clockwiseAndReorder(Path iPoints)
        {
            return pClockwiseAndReorder(iPoints);
        }

        static Path pClockwiseAndReorder(Path iPoints)
        {
            iPoints = pClockwise(iPoints);
            iPoints = pReorder(iPoints);

            return iPoints;
        }

        public static Path reOrder(Path iPoints)
        {
            return pReorder(iPoints);
        }

        static Path pReorder(Path iPoints)
        {
            Int32 reIndexStart = 0;
            Int32 minX_index = MinX(iPoints);
            Int64 minX = iPoints[minX_index].X;
            // This will reorder the point index so that the 0-indexed point is at the minimum X value, and, in the case of multiple points at min X, at the lowest Y of all of those.
            List<Int32> minXPoints = new List<Int32>();
            for (Int32 pt = 0; pt < iPoints.Count; pt++)
            {
                if (iPoints[pt].X == minX)
                {
                    minXPoints.Add(pt);
                }
            }
            // Now we need to query our minXPoints to find the point with the lowest Y value.
            Int64 minY = iPoints[minXPoints[0]].Y;
            reIndexStart = minXPoints[0];
            for (Int32 index = 1; index < minXPoints.Count; index++)
            {
                if (iPoints[minXPoints[index]].Y < minY)
                {
                    minY = iPoints[minXPoints[index]].Y;
                    reIndexStart = minXPoints[index];
                }
            }
            if (reIndexStart != 0)
            {
                Path tempList = new Path();
                // Now to start the re-indexing.
                for (Int32 pt = reIndexStart; pt < iPoints.Count; pt++)
                {
                    tempList.Add(new IntPoint(iPoints[pt].X, iPoints[pt].Y));
                }
                // Ensure we close the shape by hitting the reIndexStart point again, since we will possibly have pushed it to the beginning of the shape.
                for (Int32 pt = 0; pt <= reIndexStart; pt++)
                {
                    tempList.Add(new IntPoint(iPoints[pt].X, iPoints[pt].Y));
                }

                iPoints = tempList.ToList();
            }

            return iPoints;
        }

        public static List<GeoLibPointF[]> clockwiseAndReorder(List<GeoLibPointF[]> iPoints)
        {
            return pClockwiseAndReorder(iPoints);
        }

        static List<GeoLibPointF[]> pClockwiseAndReorder(List<GeoLibPointF[]> iPoints)
        {
            List<GeoLibPointF[]> ret = new List<GeoLibPointF[]>();
            for (int i = 0; i < iPoints.Count; i++)
            {
                ret.Add(pClockwiseAndReorder(iPoints[i]));
            }
            return ret;
        }

        public static GeoLibPointF[] clockwiseAndReorder(GeoLibPointF[] iPoints)
        {
            return pClockwiseAndReorder(iPoints);
        }

        static GeoLibPointF[] pClockwiseAndReorder(GeoLibPointF[] iPoints)
        {
            iPoints = pClockwise(iPoints);
            iPoints = pReorder(iPoints);
            return iPoints;
        }

        static GeoLibPointF[] pReorder(GeoLibPointF[] iPoints)
        {
            Int32 reIndexStart = 0;
            Int32 minX_index = MinX(iPoints);
            double minX = iPoints[minX_index].X;
            // This will reorder the point index so that the 0-indexed point is at the minimum X value, and, in the case of multiple points at min X, at the lowest Y of all of those.
            List<Int32> minXPoints = new List<Int32>();
            for (Int32 pt = 0; pt < iPoints.Length; pt++)
            {
                if (iPoints[pt].X == minX)
                {
                    minXPoints.Add(pt);
                }
            }
            // Now we need to query our minXPoints to find the point with the lowest Y value.
            double minY = iPoints[minXPoints[0]].Y;
            reIndexStart = minXPoints[0];
            for (Int32 index = 1; index < minXPoints.Count; index++)
            {
                if (iPoints[minXPoints[index]].Y < minY)
                {
                    minY = iPoints[minXPoints[index]].Y;
                    reIndexStart = minXPoints[index];
                }
            }
            if (reIndexStart != 0)
            {
                List<GeoLibPointF> tempList = new List<GeoLibPointF>();
                // Now to start the re-indexing.
                for (Int32 pt = reIndexStart; pt < iPoints.Length; pt++)
                {
                    // Avoid adding duplicate vertices
                    if (tempList.Count > 1)
                    {
                        if ((tempList[tempList.Count - 1].X == iPoints[pt].X) && (tempList[tempList.Count - 1].Y == iPoints[pt].Y))
                        {
                            continue;
                        }
                    }
                    tempList.Add(new GeoLibPointF(iPoints[pt].X, iPoints[pt].Y));
                }
                // Ensure we close the shape by hitting the reIndexStart point again, since we will possibly have pushed it to the beginning of the shape.
                for (Int32 pt = 0; pt <= reIndexStart; pt++)
                {
                    // Avoid adding duplicate vertices
                    if (tempList.Count > 1)
                    {
                        if ((tempList[tempList.Count - 1].X == iPoints[pt].X) && (tempList[tempList.Count - 1].Y == iPoints[pt].Y))
                        {
                            continue;
                        }
                    }
                    tempList.Add(new GeoLibPointF(iPoints[pt].X, iPoints[pt].Y));
                }

                iPoints = tempList.ToArray();
            }

            return iPoints;
        }

        public static List<GeoLibPoint[]> simplify(List<GeoLibPoint[]> source)
        {
            return pSimplify(source);
        }

        static List<GeoLibPoint[]> pSimplify(List<GeoLibPoint[]> source)
        {
            List<GeoLibPoint[]> ret = new List<GeoLibPoint[]>();
            for (int i = 0; i < source.Count; i++)
            {
                ret.Add(pSimplify(source[i]));
            }

            return ret;
        }
        public static GeoLibPoint[] simplify(GeoLibPoint[] iPoints)
        {
            return pSimplify(iPoints);
        }

        static GeoLibPoint[] pSimplify(GeoLibPoint[] iPoints)
        {
            List<IntPoint> ePoly = new List<IntPoint>();
            ePoly = Clipper.SimplifyPolygon(pathFromPoint(iPoints, 1))[0]; // assuming only one element
            ePoly = pStripTerminators(ePoly, false);

            List<IntPoint> iPoly = pathFromPoint(iPoints, 1);
            Clipper c = new Clipper();
            c.AddPath(iPoly, PolyType.ptClip, true);
            c.AddPath(iPoly, PolyType.ptSubject, true);
            List<List<IntPoint>> oPoly = new List<List<IntPoint>>();
            c.Execute(ClipType.ctIntersection, oPoly, PolyFillType.pftEvenOdd, PolyFillType.pftEvenOdd);

            GeoLibPoint[] working = pointFromPath(oPoly[0], 1);
            GeoLibPoint[] notWorking = pointFromPath(ePoly, 1);

            return working;
        }

        public static Paths[] getOutersAndCutters(Paths source)
        {
            return pGetOutersAndCutters(source);
        }

        static Paths[] pGetOutersAndCutters(Paths source)
        {
            Paths[] ret = new Paths[2];
            // Find cutters and outers.
            Paths outers = new Paths();
            Paths cutters = new Paths();
            for (int i = 0; i < source.Count; i++)
            {
                if (pIsClockwise(source[i]))
                {
                    outers.Add(source[i]);
                }
                else
                {
                    cutters.Add(source[i]);
                }
            }
            ret[(int)outerCutterIndex.outer] = outers;
            ret[(int)outerCutterIndex.cutter] = cutters;

            return ret;
        }

        public static Path stripColinear(Path source)
        {
            return pStripColinear(source);
        }

        static Path pStripColinear(Path source)
        {
            if (source.Count < 3)
            {
                return source;
            }

            Path ret = new Path();

            for (int pt = 0; pt < source.Count; pt++)
            {
                IntPoint interSection_A, interSection_B, interSection_C;
                // Assess angle.
                if (pt == 0)
                {
                    interSection_B = source[source.Count - 1]; // map to last point
                    interSection_C = source[pt];
                    interSection_A = source[pt + 1];
                }
                else if (pt == source.Count - 1) // last point in the list
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

                double theta = pAngleBetweenPoints(interSection_A, interSection_B, interSection_C, false);

                if ((pt == 0) || (theta != 180))
                {
                    ret.Add(new IntPoint(source[pt]));
                }
            }
            return ret;
        }

        public static List<GeoLibPoint[]> stripColinear(List<GeoLibPoint[]> source)
        {
            return pStripColinear(source);
        }

        static List<GeoLibPoint[]> pStripColinear(List<GeoLibPoint[]> source)
        {
            List<GeoLibPoint[]> ret = new List<GeoLibPoint[]>();
            for (int i = 0; i < source.Count; i++)
            {
                ret.Add(pStripColinear(source[i]));
            }

            return ret;
        }

        public static GeoLibPoint[] stripColinear(GeoLibPoint[] source)
        {
            return pStripColinear(source);
        }

        static GeoLibPoint[] pStripColinear(GeoLibPoint[] source)
        {
            if (source.Length < 3)
            {
                return source;
            }

            List<GeoLibPoint> ret = new List<GeoLibPoint>();

            for (int pt = 0; pt < source.Length; pt++)
            {
                GeoLibPoint interSection_A, interSection_B, interSection_C;
                // Assess angle.
                if (pt == 0)
                {
                    interSection_B = source[source.Length - 1]; // map to last point
                    interSection_C = source[pt];
                    interSection_A = source[pt + 1];
                }
                else if (pt == source.Length - 1) // last point in the list
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

                double theta = pAngleBetweenPoints(interSection_A, interSection_B, interSection_C, false);

                bool addPoint = true;
                if ((pt != 0) && (pt != source.Length - 1))
                {
                    if (theta == 180)
                    {
                        addPoint = false;
                    }
                }

                if (addPoint)
                {
                    ret.Add(new GeoLibPoint(source[pt]));
                }
            }
            return ret.ToArray();
        }

        public static List<GeoLibPoint> stripColinear(List<GeoLibPoint> source)
        {
            return pStripColinear(source);
        }

        static List<GeoLibPoint> pStripColinear(List<GeoLibPoint> source)
        {
            if (source.Count < 3)
            {
                return source;
            }

            List<GeoLibPoint> ret = new List<GeoLibPoint>();

            for (int pt = 0; pt < source.Count; pt++)
            {
                GeoLibPoint interSection_A, interSection_B, interSection_C;
                // Assess angle.
                if (pt == 0)
                {
                    interSection_B = source[source.Count - 1]; // map to last point
                    interSection_C = source[pt];
                    interSection_A = source[pt + 1];
                }
                else if (pt == source.Count - 1) // last point in the list
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

                double theta = pAngleBetweenPoints(interSection_A, interSection_B, interSection_C, false);

                bool addPoint = true;
                if ((pt != 0) && (pt != source.Count - 1))
                {
                    if (theta == 180)
                    {
                        addPoint = false;
                    }
                }

                if (addPoint)
                {
                    ret.Add(new GeoLibPoint(source[pt]));
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

        static List<GeoLibPoint> pRemoveDuplicates(List<GeoLibPoint> source)
        {
            List<GeoLibPoint> ret = new List<GeoLibPoint>();
            if (source.Count > 0)
            {
                ret.Add(new GeoLibPoint(source[0]));
                int retIndex = 1;
                for (int i = 1; i < source.Count - 1; i++)
                {
                    if ((source[i].X != ret[retIndex-1].X) || (source[i].Y != ret[retIndex - 1].Y))
                    {
                        ret.Add(new GeoLibPoint(source[i]));
                        retIndex++;
                    }
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

        static List<GeoLibPointF> pRemoveDuplicates(List<GeoLibPointF> source)
        {
            List<GeoLibPointF> ret = new List<GeoLibPointF>();
            if (source.Count > 0)
            {
                ret.Add(new GeoLibPointF(source[0]));
                int retIndex = 1;
                for (int i = 1; i < source.Count - 1; i++)
                {
                    if ((Math.Abs(source[i].X - ret[retIndex - 1].X) > double.Epsilon) || (Math.Abs(source[i].Y - ret[retIndex - 1].Y) > double.Epsilon))
                    {
                        ret.Add(new GeoLibPointF(source[i]));
                        retIndex++;
                    }
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

        static GeoLibPoint[] pStripTerminators(GeoLibPoint[] source, bool keepLast)
        {
            return pStripTerminators(source.ToList(), keepLast).ToArray();
        }
        static List<GeoLibPoint> pStripTerminators(List<GeoLibPoint> source, bool keepLast)
        {
            bool firstLast_same = false;
            int pt_Check = source.Count - 1;
            if (GeoWrangler.distanceBetweenPoints(source[pt_Check], source[0]) < 0.01)
            {
                firstLast_same = true; // remove duplicated points. The shape will be closed later.
            }
            while (firstLast_same)
            {
                source.RemoveAt(pt_Check); // remove duplicated points. The shape will be closed later
                pt_Check--;
                if (GeoWrangler.distanceBetweenPoints(source[pt_Check], source[0]) > 0.01)
                {
                    firstLast_same = false; // stop at the first unmatched point.
                }
            }

            if (keepLast)
            {
                source = pClose(source);
            }

            return source;
        }

        public static List<GeoLibPointF[]> stripColinear(List<GeoLibPointF[]> source)
        {
            return pStripColinear(source);
        }

        static List<GeoLibPointF[]> pStripColinear(List<GeoLibPointF[]> source)
        {
            List<GeoLibPointF[]> ret = new List<GeoLibPointF[]>();
            for (int i = 0; i < source.Count; i++)
            {
                ret.Add(pStripColinear(source[i]));
            }

            return ret;
        }

        public static GeoLibPointF[] stripColinear(GeoLibPointF[] source)
        {
            return pStripColinear(source);
        }

        static GeoLibPointF[] pStripColinear(GeoLibPointF[] source)
        {
            if (source.Length < 3)
            {
                return source;
            }

            List<GeoLibPointF> ret = new List<GeoLibPointF>();

            for (int pt = 0; pt < source.Length; pt++)
            {
                GeoLibPointF interSection_A, interSection_B, interSection_C;
                // Assess angle.
                if (pt == 0)
                {
                    interSection_B = source[source.Length - 1]; // map to last point
                    interSection_C = source[pt];
                    interSection_A = source[pt + 1];
                }
                else if (pt == source.Length - 1) // last point in the list
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

                double theta = pAngleBetweenPoints(interSection_A, interSection_B, interSection_C, false);

                bool addPoint = true;
                if ((pt != 0) && (pt != source.Length - 1))
                {
                    if (theta == 180)
                    {
                        addPoint = false;
                    }
                }

                if (addPoint)
                {
                    ret.Add(new GeoLibPointF(source[pt]));
                }
            }
            return ret.ToArray();
        }

        public static List<GeoLibPointF> stripColinear(List<GeoLibPointF> source)
        {
            return pStripColinear(source);
        }

        static List<GeoLibPointF> pStripColinear(List<GeoLibPointF> source)
        {
            if (source.Count < 3)
            {
                return source;
            }

            List<GeoLibPointF> ret = new List<GeoLibPointF>();

            for (int pt = 0; pt < source.Count; pt++)
            {
                GeoLibPointF interSection_A, interSection_B, interSection_C;
                // Assess angle.
                if (pt == 0)
                {
                    interSection_B = source[source.Count - 1]; // map to last point
                    interSection_C = source[pt];
                    interSection_A = source[pt + 1];
                }
                else if (pt == source.Count - 1) // last point in the list
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

                double theta = pAngleBetweenPoints(interSection_A, interSection_B, interSection_C, false);

                bool addPoint = true;
                if ((pt != 0) && (pt != source.Count - 1))
                {
                    if (theta == 180)
                    {
                        addPoint = false;
                    }
                }

                if (addPoint)
                {
                    ret.Add(new GeoLibPointF(source[pt]));
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

        static GeoLibPointF[] pStripTerminators(GeoLibPointF[] source, bool keepLast)
        {
            return pStripTerminators(source.ToList(), keepLast).ToArray();
        }
        static List<GeoLibPointF> pStripTerminators(List<GeoLibPointF> source, bool keepLast)
        {
            bool firstLast_same = false;
            int pt_Check = source.Count - 1;
            if (GeoWrangler.distanceBetweenPoints(source[pt_Check], source[0]) < 0.01)
            {
                firstLast_same = true; // remove duplicated points. The shape will be closed later.
            }
            while (firstLast_same)
            {
                source.RemoveAt(pt_Check); // remove duplicated points. The shape will be closed later
                pt_Check--;
                if (GeoWrangler.distanceBetweenPoints(source[pt_Check], source[0]) > 0.01)
                {
                    firstLast_same = false; // stop at the first unmatched point.
                }
            }

            if (keepLast)
            {
                source = pClose(source);
            }

            return source;
        }

        public static Path stripTerminators(Path source, bool keepLast)
        {
            return pStripTerminators(source, keepLast);
        }

        static Path pStripTerminators(Path source, bool keepLast)
        {
            if (source.Count <= 1)
            {
                return source;
            }

            bool firstLast_same = false;
            int pt_Check = source.Count - 1;
            if (GeoWrangler.distanceBetweenPoints(source[pt_Check], source[0]) < 10)
            {
                firstLast_same = true; // remove duplicated points. The shape will be closed later.
            }
            while (firstLast_same)
            {
                source.RemoveAt(pt_Check); // remove duplicated points. The shape will be closed later
                pt_Check--;
                if (source.Count < 1)
                {
                    return source;
                }
                if (GeoWrangler.distanceBetweenPoints(source[pt_Check], source[0]) > 10)
                {
                    firstLast_same = false; // stop at the first unmatched point.
                }
            }

            if (keepLast)
            {
                source = pClose(source);
            }

            return source;
        }

        public static Paths close(Paths source)
        {
            return pClose(source);
        }

        static Paths pClose(Paths source)
        {
            Paths ret = new Paths();
            for (int i = 0; i < source.Count; i++)
            {
                ret.Add(pClose(source[i]));
            }
            return ret;
        }
        public static Path close(Path source)
        {
            return pClose(source);
        }

        static Path pClose(Path source)
        {
            if (source.Count < 1)
            {
                return source;
            }
            if ((source[0].X != source[source.Count - 1].X) || (source[0].Y != source[source.Count - 1].Y))
            {
                source.Add(new IntPoint(source[0]));
            }
            return source;
        }

        public static List<GeoLibPoint> close(List<GeoLibPoint> source)
        {
            return pClose(source);
        }

        static List<GeoLibPoint> pClose(List<GeoLibPoint> source)
        {
            if (source.Count < 1)
            {
                return source;
            }
            if ((source[0].X != source[source.Count - 1].X) || (source[0].Y != source[source.Count - 1].Y))
            {
                source.Add(new GeoLibPoint(source[0]));
            }
            return source;
        }

        public static GeoLibPoint[] close(GeoLibPoint[] source)
        {
            return pClose(source);
        }

        static GeoLibPoint[] pClose(GeoLibPoint[] source)
        {
            if (source.Length < 1)
            {
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

        static List<GeoLibPointF> pClose(List<GeoLibPointF> source)
        {
            if (source.Count < 1)
            {
                return source;
            }
            if ((source[0].X != source[source.Count - 1].X) || (source[0].Y != source[source.Count - 1].Y))
            {
                source.Add(new GeoLibPointF(source[0]));
            }
            return source;
        }

        public static GeoLibPointF[] close(GeoLibPointF[] source)
        {
            return pClose(source);
        }

        static GeoLibPointF[] pClose(GeoLibPointF[] source)
        {
            if (source.Length < 1)
            {
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

        static Paths pFromSoup(Paths source, WindingRule wr)
        {
            Paths outers = new Paths();
            Paths cutters = new Paths();

            for (int i = 0; i < source.Count; i++)
            {
                if (Clipper.Orientation(source[i]) == Clipper.Orientation(source[0]))
                {
                    outers.Add(new Path(source[i]));
                }
                else
                {
                    cutters.Add(new Path(source[i]));
                }
            }

            // Set up our contours. Since Clipper sets up the subject as the first item, we'll make that clockwise and force the rest to counterclockwise.	
            var tess = new Tess();

            for (int poly = 0; poly < outers.Count; poly++)
            {
                // Skip tiny fragments. The tessellator has trouble with them.	
                /*
                if (Math.Abs(Clipper.Area(outers[poly])) < 1E-16)
                {
                    continue;
                }
                */
                ContourVertex[] contour = new ContourVertex[outers[poly].Count];
                for (int i = 0; i < contour.Length; i++)
                {
                    contour[i].Position = new Vec3 { X = outers[poly][i].X, Y = outers[poly][i].Y, Z = 0 };
                }

                tess.AddContour(contour);
            }

            for (int poly = 0; poly < cutters.Count; poly++)
            {
                // Skip tiny fragments. The tessellator has trouble with them.
                /*
                if (Math.Abs(Clipper.Area(cutters[poly])) < 1E-16)
                {
                    continue;
                }
                */
                ContourVertex[] contour = new ContourVertex[cutters[poly].Count];
                for (int i = 0; i < contour.Length; i++)
                {
                    contour[i].Position = new Vec3 { X = cutters[poly][i].X, Y = cutters[poly][i].Y, Z = 0 };
                }

                tess.AddContour(contour, ContourOrientation.CounterClockwise);
            }

            try
            {

                Paths retPaths = new Paths();

                // Triangulate. This gives us a triangle soup, which may not be entirely helpful for the case where we're punching multiple holes into a mesh. For now, this works, but the limitation will need to be reviewd.	
                // It might be that a PolyTree is needed as the input to the keyhole for these complex cases.... or clockwise paths may need to be evaluated piecewise using all counterclockwise paths for the tessellation.	
                int polysize = 3;
                tess.Tessellate(wr, ElementType.Polygons, polysize);

                // Iterate triangles and create output geometry. We'll use clipper to simplify the output geometry.	
                Clipper c = new Clipper();
                c.PreserveCollinear = true;
                retPaths = new Paths();

                Paths cPaths = new Paths();
                Paths aPaths = new Paths();

                for (int i = 0; i < tess.ElementCount; i++)
                {
                    Path trianglePath = new Path();
                    for (int p = 0; p < polysize; p++)
                    {
                        IntPoint tmpPt = new IntPoint((Int64)tess.Vertices[tess.Elements[(i * polysize) + p]].Position.X, (Int64)tess.Vertices[tess.Elements[(i * polysize) + p]].Position.Y);
                        trianglePath.Add(tmpPt);
                    }

                    if (Clipper.Orientation(trianglePath))
                    {
                        cPaths.Add(trianglePath.ToList());
                    }
                    else
                    {
                        aPaths.Add(trianglePath.ToList());
                    }
                }

                // Add paths to the clipper.	
                c.AddPaths(cPaths, PolyType.ptSubject, true);
                c.AddPaths(aPaths, PolyType.ptClip, true);

                c.Execute(ClipType.ctUnion, retPaths, PolyFillType.pftNonZero, PolyFillType.pftNonZero);

                retPaths = pClose(retPaths);

                return retPaths;
            }
            catch (Exception)
            {
                return source;
            }
        }

        public static List<GeoLibPointF[]> clean_and_flatten(List<GeoLibPointF[]> source, long scaling)
        {
            return pClean_and_flatten(source, scaling);
        }
        static List<GeoLibPointF[]> pClean_and_flatten(List<GeoLibPointF[]> source, long scaling)
        {
            Paths sourcePaths = pPathsFromPointFs(source, scaling);
            Clipper c = new Clipper();
            c.AddPaths(sourcePaths, PolyType.ptSubject, true);
            Paths solution = new Paths();
            c.Execute(ClipType.ctUnion, solution);

            Paths keyHoled = pMakeKeyHole(solution);

            return pPointFsFromPaths(pClockwiseAndReorder(keyHoled), scaling);
        }
    }
}
