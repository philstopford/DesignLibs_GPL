using ClipperLib;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using utility;

namespace geoWrangler
{
    using Path = List<IntPoint>;
    using Paths = List<List<IntPoint>>;

    public class RayCast
    {
        Paths clippedLines;
        Paths castLines;

        public enum forceSingleDirection { no, vertical, horizontal } // No and horizontal are treated the same in this code at the moment; leaving these to make the code more readable and intent clearer.

        void prox_ZFillCallback(IntPoint bot1, IntPoint top1, IntPoint bot2, IntPoint top2, ref IntPoint pt)
        {
            pt.Z = bot1.Z;
        }

        public static readonly List<string> fallOffList = new List<string> {"None", "Linear", "Gaussian", "Cosine"};
        public enum falloff { none, linear, gaussian, cosine }

        // Corner projection is default and takes an orthogonal ray out from the corner. Setting to false causes an averaged normal to be generated.

        public RayCast(Path emissionPath, Paths collisionPaths, Int64 max, bool projectCorners = true, bool invert = false, int multisampleRayCount = 0, bool runOuterLoopThreaded = false, bool runInnerLoopThreaded = false, IntPoint startOffset = new IntPoint(), IntPoint endOffset = new IntPoint(), falloff sideRayFallOff = falloff.none, double sideRayFallOffMultiplier = 1.0f, forceSingleDirection dirOverride = forceSingleDirection.no)
        {
            rayCast(emissionPath, collisionPaths, max, projectCorners, invert, multisampleRayCount, runOuterLoopThreaded, runInnerLoopThreaded, startOffset, endOffset, sideRayFallOff, sideRayFallOffMultiplier, dirOverride);
        }

        public RayCast(Path emissionPath, Path collisionPath, Int64 max, bool projectCorners = true, bool invert = false, int multisampleRayCount = 0, bool runOuterLoopThreaded = false, bool runInnerLoopThreaded = false, IntPoint startOffset = new IntPoint(), IntPoint endOffset = new IntPoint(), falloff sideRayFallOff = falloff.none, double sideRayFallOffMultiplier = 1.0f, forceSingleDirection dirOverride = forceSingleDirection.no)
        {
            rayCast(emissionPath, new Paths { collisionPath }, max, projectCorners, invert, multisampleRayCount, runOuterLoopThreaded, runInnerLoopThreaded, startOffset, endOffset, sideRayFallOff, sideRayFallOffMultiplier, dirOverride);
        }

        public Paths getRays()
        {
            return pGetRays();
        }

        Paths pGetRays()
        {
            return castLines;
        }

        public Paths getClippedRays()
        {
            return pGetClippedRays();
        }

        Paths pGetClippedRays()
        {
            return clippedLines;
        }

        public double getRayLength(int ray)
        {
            if ((ray < 0) || (ray > clippedLines.Count))
            {
                return -1.0f;
            }
            return pGetRayLength(ray);
        }

        double pGetRayLength(int ray)
        {
            return GeoWrangler.distanceBetweenPoints(clippedLines[ray][0], clippedLines[ray][clippedLines[ray].Count - 1]);
        }

        void rayCast(Path emissionPath, Paths collisionPaths, Int64 maxRayLength, bool projectCorners, bool invert, int multisampleRayCount, bool runOuterLoopThreaded, bool runInnerLoopThreaded, IntPoint startOffset, IntPoint endOffset, falloff sideRayFallOff, double sideRayFallOffMultiplier, forceSingleDirection dirOverride)
        {
            // Setting this to true, we shorten rays with the falloff. False means we reduce the contribution to the average instead.
            bool truncateRaysByWeight = false;

            int ptCount = emissionPath.Count;

            // Due to threading and need to tie to polygon point order, we have to use these local storage options and will do the conversion at the end.
            object castLinesLock = new object();
            Paths[] castLines_ = new Paths[ptCount];
            object clippedLinesLock = new object();
            Path[] clippedLines_ = new Path[ptCount];

            // We need to think about the end point case.
            bool closedPathEmitter = (ptCount > 3) && ((emissionPath[0].X == emissionPath[ptCount - 1].X) && (emissionPath[0].Y == emissionPath[ptCount - 1].Y));

            // Get average angle for this vertex based on angles from line segments.
            // http://stackoverflow.com/questions/1243614/how-do-i-calculate-the-normal-vector-of-a-line-segment

            // Pre-calculate these for the threading to be an option.
            // This is a serial evaluation as we need both the previous and the current normal for each point.
            IntPoint[] normals = new IntPoint[ptCount];
            IntPoint[] previousNormals = new IntPoint[ptCount];
            for (Int32 pt = 0; pt < ptCount; pt++)
            {
                // Start point
                Int64 dx;
                Int64 dy;
                if (pt == emissionPath.Count - 1)
                {
                    // Last matches the first. Since we flip the dx and dy tone later, we need to compensate here.
                    if (closedPathEmitter)
                    {
                        dx = -normals[0].X;
                        dy = -normals[0].Y;
                    }
                    else
                    {
                        dx = emissionPath[ptCount - 1].X - endOffset.X;
                        dy = emissionPath[ptCount - 1].Y - endOffset.Y;
                    }
                }
                else if ((!closedPathEmitter) && (pt == 0))
                {
                    dx = emissionPath[0].X - startOffset.X;
                    dy = emissionPath[0].Y - startOffset.Y;
                }
                else
                {
                    dx = emissionPath[pt + 1].X - emissionPath[pt].X;
                    dy = emissionPath[pt + 1].Y - emissionPath[pt].Y;
                }
                normals[pt] = new IntPoint(-dx, -dy);

                // Previous normal
                if (pt == 0)
                {
                    if (closedPathEmitter)
                    {
                        // n-1 identical to the 0-th point, so we need to dig a little deeper.
                        dx = emissionPath[0].X - emissionPath[ptCount - 2].X;
                        dy = emissionPath[0].Y - emissionPath[ptCount - 2].Y;
                    }
                    else
                    {
                        dx = emissionPath[0].X - startOffset.X;
                        dy = emissionPath[0].Y - startOffset.Y;
                    }

                    previousNormals[pt] = new IntPoint(-dx, -dy);
                }
                else
                {
                    previousNormals[pt] = (normals[pt - 1]);
                }
            }

            ParallelOptions po_outer = new ParallelOptions();
            if (!runOuterLoopThreaded)
            {
                po_outer.MaxDegreeOfParallelism = 1;
            }

            ParallelOptions po_inner = new ParallelOptions();
            if (!runInnerLoopThreaded)
            {
                po_inner.MaxDegreeOfParallelism = 1;
            }

            Parallel.For(0, ptCount, po_outer, pt =>
            {
                IntPoint currentEdgeNormal = normals[pt];
                IntPoint previousEdgeNormal = previousNormals[pt];

                IntPoint averagedEdgeNormal;

                IntPoint startPoint = new IntPoint(emissionPath[pt]);

                // Get average angle for this vertex based on angles from line segments.
                // http://stackoverflow.com/questions/1243614/how-do-i-calculate-the-normal-vector-of-a-line-segment

                if (projectCorners && (((currentEdgeNormal.X == 0) && (previousEdgeNormal.Y == 0)) ||
                                       ((currentEdgeNormal.Y == 0) && (previousEdgeNormal.X == 0))))
                {
                    Int64 tX = currentEdgeNormal.X;
                    Int64 tY = currentEdgeNormal.Y;
                    // If we're traversing a 90 degree corner, let's not project a diagonal, but fix on our current edge normal.
                    if (!invert || (dirOverride == forceSingleDirection.vertical))
                    {
                        tX = -tX;
                        tY = -tY;
                    }
                    averagedEdgeNormal = new IntPoint(tX, tY);
                }
                else
                {
                    if (invert)
                    {
                        currentEdgeNormal = new IntPoint(currentEdgeNormal.X, -currentEdgeNormal.Y);
                        previousEdgeNormal = new IntPoint(previousEdgeNormal.X, -previousEdgeNormal.Y);
                    }
                    // Average out our normals
                    averagedEdgeNormal = new IntPoint((previousEdgeNormal.X + currentEdgeNormal.X) / 2, (previousEdgeNormal.Y + currentEdgeNormal.Y) / 2);
                }

                // Normalization. We don't change the original vectors to avoid having to normalize everywhere.
                double length = Math.Sqrt(Utils.myPow(averagedEdgeNormal.X, 2) + Utils.myPow(averagedEdgeNormal.Y, 2));

                double endPointDeltaX = 0;
                double endPointDeltaY = 0;

                // Avoid div-by-zero; 0-length is unimportant. Note that setting this cut-off too high produces artifacts.
                if (length > 0.0001)
                {
                    endPointDeltaX = Convert.ToDouble(averagedEdgeNormal.X) / length;
                    endPointDeltaY = Convert.ToDouble(averagedEdgeNormal.Y) / length;
                }

                // Set to max ray length from callsite.
                endPointDeltaX *= maxRayLength;
                endPointDeltaY *= maxRayLength;

                if (invert)
                {
                    endPointDeltaY *= -1;
                }
                endPointDeltaX *= -1;

                IntPoint endPoint = new IntPoint(endPointDeltaY + startPoint.X, endPointDeltaX + startPoint.Y);

                if (sideRayFallOff != falloff.none)
                {
                    endPoint.Z = (Int64)1E4;
                }

                Paths rays = new Paths();
                Path line = new Path {new IntPoint(startPoint), new IntPoint(endPoint)};
                rays.Add(line/*.ToList()*/);

                double angleStep = 90.0f / (1 + multisampleRayCount);

                for (int sample = 0; sample < multisampleRayCount; sample++)
                {
                    // Add more samples, each n-degrees rotated from the nominal ray
                    double rayAngle = (sample + 1) * angleStep;

                    IntPoint endPoint_f = endPoint;

                    double weight_val = 1.0f;
                    switch (sideRayFallOff)
                    {
                        // Gaussian fall-off
                        case falloff.gaussian:
                            weight_val = Math.Exp(-(Math.Pow(sideRayFallOffMultiplier * (rayAngle / 90.0f), 2)));
                            break;
                        // Linear fall-off
                        case falloff.linear:
                            weight_val = 1.0f - Math.Min(rayAngle / 90.0f, 1.0f);
                            break;
                        // Cosine fall-off
                        case falloff.cosine:
                            double angle = sideRayFallOffMultiplier * rayAngle;
                            if (angle > 90.0)
                            {
                                angle = 90.0;
                            }
                            if (angle < 0)
                            {
                                angle = 0;
                            }
                            weight_val = Math.Cos(Utils.toRadians(angle));
                            // Shift up and flatten to 0-1 range.
                            weight_val += 1.0f;
                            weight_val *= 0.5;
                            break;
                        // No falloff
                    }

                    if (truncateRaysByWeight)
                    {
                        endPoint_f = new IntPoint(startPoint.X + (weight_val * endPointDeltaY), startPoint.Y + (weight_val * endPointDeltaX));
                    }

                    IntPoint sPoint = new IntPoint(startPoint.X, startPoint.Y);

                    if (sideRayFallOff != falloff.none)
                    {
                        endPoint_f.Z = Convert.ToInt64(weight_val * 1E4);
                        sPoint.Z = endPoint_f.Z;
                    }
                    IntPoint endPoint1 = GeoWrangler.Rotate(startPoint, endPoint_f, rayAngle);
                    IntPoint endPoint2 = GeoWrangler.Rotate(startPoint, endPoint_f, -rayAngle);

                    // The order of line1 below is important, but I'm not yet sure why. If you change it, the expansion becomes asymmetrical on a square (lower section gets squashed).
                    Path line1 = new Path {new IntPoint(endPoint1), new IntPoint(sPoint)};
                    rays.Add(line1);
                    Path line2 = new Path {new IntPoint(sPoint), new IntPoint(endPoint2)};
                    rays.Add(line2);
                }

                Monitor.Enter(castLinesLock);
                try
                {
                    castLines_[pt] = new Paths();
                    castLines_[pt].AddRange(rays);
                }
                finally
                {
                    Monitor.Exit(castLinesLock);
                }

                Int64[] resultX = new Int64[rays.Count];
                Int64[] resultY = new Int64[rays.Count];
                double[] weight = new double[rays.Count];

                previousEdgeNormal = new IntPoint(currentEdgeNormal.X, currentEdgeNormal.Y);

                object resultLock = new object();
                Parallel.For(0, rays.Count, po_inner, ray =>
                    {
                        Clipper d = new Clipper();
                        if (sideRayFallOff != falloff.none)
                        {
                            d.ZFillFunction = prox_ZFillCallback;
                        }
                        d.AddPath(rays[ray], PolyType.ptSubject, false);
                        d.AddPaths(collisionPaths, PolyType.ptClip, true);
                        PolyTree polyTree = new PolyTree();
                        if (invert)
                        {
                            d.Execute(ClipType.ctIntersection, polyTree);
                        }
                        else
                        {
                            d.Execute(ClipType.ctDifference, polyTree);
                        }

                        // There is no matching order in the return output here, so we have to take this odd approach.
                        Paths tmpLine = Clipper.OpenPathsFromPolyTree(polyTree);

                        int tmpLineCount = tmpLine.Count;

                        // If we got no result, let's get that sorted out. End and start points are the same.
                        if (tmpLineCount == 0)
                        {
                            Monitor.Enter(resultLock);
                            try
                            {
                                resultX[ray] = startPoint.X;
                                resultY[ray] = startPoint.Y;
                                weight[ray] = 1.0f;
                            }
                            finally
                            {
                                Monitor.Exit(resultLock);
                            }
                        }

                        if (tmpLineCount > 1)
                        {
                            // We got two lines back. Need to check this carefully.
                            // We need to find the line that has a start/end point matched to our origin point for the ray.
                            int index = -1;
                            for (int tL = 0; tL < tmpLine.Count; tL++)
                            {
                                Int64 tL0X = tmpLine[tL][0].X;
                                Int64 tL0Y = tmpLine[tL][0].Y;
                                Int64 tL1X = tmpLine[tL][1].X;
                                Int64 tL1Y = tmpLine[tL][1].Y;
                                if (((tL0X == startPoint.X) && (tL0Y == startPoint.Y)) || ((tL1X == startPoint.X) && (tL1Y == startPoint.Y)))
                                {
                                    index = tL;
                                    break;
                                }
                            }
                            Path tPath = new Path();
                            if (index >= 0)
                            {
                                tPath = tmpLine[index];
                            }
                            else
                            {
                                tPath.Add(startPoint);
                                tPath.Add(startPoint);
                            }
                            tmpLine.Clear();
                            tmpLine.Add(tPath);
                        }

                        tmpLineCount = tmpLine.Count;

                        for (int tL = 0; tL < tmpLineCount; tL++)
                        {
                            // Figure out which end of the result line matches our origin point.
                            if ((tmpLine[tL][0].X == startPoint.X) && (tmpLine[tL][0].Y == startPoint.Y))
                            {
                                Monitor.Enter(resultLock);
                                try
                                {
                                    resultX[ray] = tmpLine[tL][1].X;
                                    resultY[ray] = tmpLine[tL][1].Y;
                                    weight[ray] = Convert.ToDouble(tmpLine[tL][1].Z) / 1E4;
                                }
                                finally
                                {
                                    Monitor.Exit(resultLock);
                                }
                            }
                            else if ((tmpLine[tL][1].X == startPoint.X) && (tmpLine[tL][1].Y == startPoint.Y))
                            {
                                Monitor.Enter(resultLock);
                                try
                                {
                                    // Clipper reversed the line direction, so we need to deal with this.
                                    resultX[ray] = tmpLine[tL][0].X;
                                    resultY[ray] = tmpLine[tL][0].Y;
                                    weight[ray] = Convert.ToDouble(tmpLine[tL][0].Z) / 1E4;
                                }
                                finally
                                {
                                    Monitor.Exit(resultLock);
                                }
                            }
                        }
                    }
                );

                Path resultPath = new Path {startPoint};

                int xCount = 0;
                Int64 xAv = 0;
                int yCount = 0;
                Int64 yAv = 0;

                // If we are not truncating by weight, we do not need to average here - it was done with the normalization above.
                if (sideRayFallOff == falloff.none)
                {
                    // Average the result to give a weighted spacing across the rays.
                    for (int result = 0; result < resultX.Length; result++)
                    {
                        if (Math.Abs(resultX[result]) > 1000)
                        {
                            xCount++;
                            xAv += resultX[result];
                        }
                        if (Math.Abs(resultY[result]) > 1000)
                        {
                            yCount++;
                            yAv += resultY[result];
                        }
                    }

                    if (xCount != 0)
                    {
                        xAv /= xCount;
                    }
                    if (yCount != 0)
                    {
                        yAv /= yCount;
                    }
                }
                else
                {
                    double totalWeight = 0.0f;
                    foreach (double t in weight)
                    {
                        totalWeight += t;
                    }

                    // Average the result to give a weighted spacing across the rays.
                    for (int w = 0; w < weight.Length; w++)
                    {
                        double weight_ = 1.0f;
                        if (sideRayFallOff != falloff.none && (!truncateRaysByWeight) && (totalWeight > 0))
                        {
                            weight_ = weight[w] / totalWeight;
                        }

                        if (Math.Abs(resultX[w]) > 1000)
                        {
                            xCount++;
                            xAv += Convert.ToInt64(weight_ * resultX[w]);
                        }
                        if (Math.Abs(resultY[w]) > 1000)
                        {
                            yCount++;
                            yAv += Convert.ToInt64(weight_ * resultY[w]);
                        }
                    }
                }

                resultPath.Add(new IntPoint(xAv, yAv));
                Monitor.Enter(clippedLinesLock);
                try
                {
                    clippedLines_[pt] = resultPath;
                }
                finally
                {
                    Monitor.Exit(clippedLinesLock);
                }
            });

            // Convert the array back to a list.
            clippedLines = clippedLines_.ToList();
            castLines = new Paths();
            foreach (Paths t in castLines_)
            {
                castLines.AddRange(t);
            }
        }
    }
}
