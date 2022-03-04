using ClipperLib;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using utility;

namespace geoWrangler;

using Path = List<IntPoint>;
using Paths = List<List<IntPoint>>;

public class RayCast
{
    private Paths clippedLines;
    private Paths castLines;

    public enum forceSingleDirection { no, vertical, horizontal } // No and horizontal are treated the same in this code at the moment; leaving these to make the code more readable and intent clearer.

    private void prox_ZFillCallback(IntPoint bot1, IntPoint top1, IntPoint bot2, IntPoint top2, ref IntPoint pt)
    {
        pt.Z = bot1.Z;
    }

    public static readonly List<string> fallOffList = new() {"None", "Linear", "Gaussian", "Cosine"};
    public enum falloff { none, linear, gaussian, cosine }

    // Corner projection is default and takes an orthogonal ray out from the corner. Setting to false causes an averaged normal to be generated.

    public RayCast(Path emissionPath, Paths collisionPaths, long max, bool projectCorners = true, bool invert = false, int multisampleRayCount = 0, bool runOuterLoopThreaded = false, bool runInnerLoopThreaded = false, IntPoint startOffset = new(), IntPoint endOffset = new(), falloff sideRayFallOff = falloff.none, double sideRayFallOffMultiplier = 1.0f, forceSingleDirection dirOverride = forceSingleDirection.no)
    {
        rayCast(emissionPath, collisionPaths, max, projectCorners, invert, multisampleRayCount, runOuterLoopThreaded, runInnerLoopThreaded, startOffset, endOffset, sideRayFallOff, sideRayFallOffMultiplier, dirOverride);
    }

    public RayCast(Path emissionPath, Path collisionPath, long max, bool projectCorners = true, bool invert = false, int multisampleRayCount = 0, bool runOuterLoopThreaded = false, bool runInnerLoopThreaded = false, IntPoint startOffset = new(), IntPoint endOffset = new(), falloff sideRayFallOff = falloff.none, double sideRayFallOffMultiplier = 1.0f, forceSingleDirection dirOverride = forceSingleDirection.no)
    {
        rayCast(emissionPath, new Paths { collisionPath }, max, projectCorners, invert, multisampleRayCount, runOuterLoopThreaded, runInnerLoopThreaded, startOffset, endOffset, sideRayFallOff, sideRayFallOffMultiplier, dirOverride);
    }

    public Paths getRays()
    {
        return pGetRays();
    }

    private Paths pGetRays()
    {
        return castLines;
    }

    public Paths getClippedRays()
    {
        return pGetClippedRays();
    }

    private Paths pGetClippedRays()
    {
        return clippedLines;
    }

    public double getRayLength(int ray)
    {
        if (ray < 0 || ray > clippedLines.Count)
        {
            return -1.0f;
        }
        return pGetRayLength(ray);
    }

    private double pGetRayLength(int ray)
    {
        return GeoWrangler.distanceBetweenPoints(clippedLines[ray][0], clippedLines[ray][clippedLines[ray].Count - 1]);
    }

    private void rayCast(Path emissionPath, Paths collisionPaths, long maxRayLength, bool projectCorners, bool invert, int multisampleRayCount, bool runOuterLoopThreaded, bool runInnerLoopThreaded, IntPoint startOffset, IntPoint endOffset, falloff sideRayFallOff, double sideRayFallOffMultiplier, forceSingleDirection dirOverride)
    {
        // Setting this to true, we shorten rays with the falloff. False means we reduce the contribution to the average instead.
        const bool truncateRaysByWeight = false;

        int ptCount = emissionPath.Count;

        // Due to threading and need to tie to polygon point order, we have to use these local storage options and will do the conversion at the end.
        object castLinesLock = new();
        Paths[] castLines_ = new Paths[ptCount];
        object clippedLinesLock = new();
        Path[] clippedLines_ = new Path[ptCount];

        // We need to think about the end point case.
        bool closedPathEmitter = ptCount > 3 && emissionPath[0].X == emissionPath[ptCount - 1].X && emissionPath[0].Y == emissionPath[ptCount - 1].Y;

        // Get average angle for this vertex based on angles from line segments.
        // http://stackoverflow.com/questions/1243614/how-do-i-calculate-the-normal-vector-of-a-line-segment

        // Pre-calculate these for the threading to be an option.
        // This is a serial evaluation as we need both the previous and the current normal for each point.
        IntPoint[] normals = new IntPoint[ptCount];
        IntPoint[] previousNormals = new IntPoint[ptCount];
        for (int pt = 0; pt < ptCount; pt++)
        {
            // Start point
            long dx;
            long dy;
            if (pt == emissionPath.Count - 1)
            {
                switch (closedPathEmitter)
                {
                    // Last matches the first. Since we flip the dx and dy tone later, we need to compensate here.
                    case true:
                        dx = -normals[0].X;
                        dy = -normals[0].Y;
                        break;
                    default:
                        dx = emissionPath[ptCount - 1].X - endOffset.X;
                        dy = emissionPath[ptCount - 1].Y - endOffset.Y;
                        break;
                }
            }
            else
            {
                switch (closedPathEmitter)
                {
                    case false when pt == 0:
                        dx = emissionPath[0].X - startOffset.X;
                        dy = emissionPath[0].Y - startOffset.Y;
                        break;
                    default:
                        dx = emissionPath[pt + 1].X - emissionPath[pt].X;
                        dy = emissionPath[pt + 1].Y - emissionPath[pt].Y;
                        break;
                }
            }

            normals[pt] = new IntPoint(-dx, -dy);

            switch (pt)
            {
                // Previous normal
                case 0:
                {
                    switch (closedPathEmitter)
                    {
                        case true:
                            // n-1 identical to the 0-th point, so we need to dig a little deeper.
                            dx = emissionPath[0].X - emissionPath[ptCount - 2].X;
                            dy = emissionPath[0].Y - emissionPath[ptCount - 2].Y;
                            break;
                        default:
                            dx = emissionPath[0].X - startOffset.X;
                            dy = emissionPath[0].Y - startOffset.Y;
                            break;
                    }

                    previousNormals[pt] = new IntPoint(-dx, -dy);
                    break;
                }
                default:
                    previousNormals[pt] = normals[pt - 1];
                    break;
            }
        }

        ParallelOptions po_outer = new();
        po_outer.MaxDegreeOfParallelism = runOuterLoopThreaded switch
        {
            false => 1,
            _ => po_outer.MaxDegreeOfParallelism
        };

        ParallelOptions po_inner = new();
        po_inner.MaxDegreeOfParallelism = runInnerLoopThreaded switch
        {
            false => 1,
            _ => po_inner.MaxDegreeOfParallelism
        };

        Parallel.For(0, ptCount, po_outer, pt =>
        {
            IntPoint currentEdgeNormal = normals[pt];
            IntPoint previousEdgeNormal = previousNormals[pt];

            IntPoint averagedEdgeNormal;

            IntPoint startPoint = new(emissionPath[pt]);

            // Get average angle for this vertex based on angles from line segments.
            // http://stackoverflow.com/questions/1243614/how-do-i-calculate-the-normal-vector-of-a-line-segment

            switch (projectCorners)
            {
                case true when currentEdgeNormal.X == 0 && previousEdgeNormal.Y == 0 ||
                               currentEdgeNormal.Y == 0 && previousEdgeNormal.X == 0:
                {
                    long tX = currentEdgeNormal.X;
                    long tY = currentEdgeNormal.Y;
                    // If we're traversing a 90 degree corner, let's not project a diagonal, but fix on our current edge normal.
                    if (!invert || dirOverride == forceSingleDirection.vertical)
                    {
                        tX = -tX;
                        tY = -tY;
                    }
                    averagedEdgeNormal = new IntPoint(tX, tY);
                    break;
                }
                default:
                {
                    switch (invert)
                    {
                        case true:
                            currentEdgeNormal = new IntPoint(currentEdgeNormal.X, -currentEdgeNormal.Y);
                            previousEdgeNormal = new IntPoint(previousEdgeNormal.X, -previousEdgeNormal.Y);
                            break;
                    }
                    // Average out our normals
                    averagedEdgeNormal = new IntPoint((previousEdgeNormal.X + currentEdgeNormal.X) / 2, (previousEdgeNormal.Y + currentEdgeNormal.Y) / 2);
                    break;
                }
            }

            // Normalization. We don't change the original vectors to avoid having to normalize everywhere.
            double length = Math.Sqrt(Utils.myPow(averagedEdgeNormal.X, 2) + Utils.myPow(averagedEdgeNormal.Y, 2));

            double endPointDeltaX = 0;
            double endPointDeltaY = 0;

            switch (length)
            {
                // Avoid div-by-zero; 0-length is unimportant. Note that setting this cut-off too high produces artifacts.
                case > 0.0001:
                    endPointDeltaX = Convert.ToDouble(averagedEdgeNormal.X) / length;
                    endPointDeltaY = Convert.ToDouble(averagedEdgeNormal.Y) / length;
                    break;
            }

            // Set to max ray length from callsite.
            endPointDeltaX *= maxRayLength;
            endPointDeltaY *= maxRayLength;

            switch (invert)
            {
                case true:
                    endPointDeltaY *= -1;
                    break;
            }
            endPointDeltaX *= -1;

            IntPoint endPoint = new(endPointDeltaY + startPoint.X, endPointDeltaX + startPoint.Y);

            if (sideRayFallOff != falloff.none)
            {
                endPoint.Z = (long)1E4;
            }

            Paths rays = new();
            Path line = new() {new IntPoint(startPoint), new IntPoint(endPoint)};
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
                        weight_val = Math.Exp(-Math.Pow(sideRayFallOffMultiplier * (rayAngle / 90.0f), 2));
                        break;
                    // Linear fall-off
                    case falloff.linear:
                        weight_val = 1.0f - Math.Min(rayAngle / 90.0f, 1.0f);
                        break;
                    // Cosine fall-off
                    case falloff.cosine:
                        double angle = sideRayFallOffMultiplier * rayAngle;
                        angle = angle switch
                        {
                            < 0 => 0,
                            _ => angle switch
                            {
                                > 90.0 => 90.0,
                                _ => angle
                            }
                        };
                        weight_val = Math.Cos(Utils.toRadians(angle));
                        // Shift up and flatten to 0-1 range.
                        weight_val += 1.0f;
                        weight_val *= 0.5;
                        break;
                    // No falloff
                }

                endPoint_f = truncateRaysByWeight switch
                {
                    true => new IntPoint(startPoint.X + weight_val * endPointDeltaY,
                        startPoint.Y + weight_val * endPointDeltaX),
                    _ => endPoint_f
                };

                IntPoint sPoint = new(startPoint.X, startPoint.Y);

                if (sideRayFallOff != falloff.none)
                {
                    endPoint_f.Z = Convert.ToInt64(weight_val * 1E4);
                    sPoint.Z = endPoint_f.Z;
                }
                IntPoint endPoint1 = GeoWrangler.Rotate(startPoint, endPoint_f, rayAngle);
                IntPoint endPoint2 = GeoWrangler.Rotate(startPoint, endPoint_f, -rayAngle);

                // The order of line1 below is important, but I'm not yet sure why. If you change it, the expansion becomes asymmetrical on a square (lower section gets squashed).
                Path line1 = new() {new IntPoint(endPoint1), new IntPoint(sPoint)};
                rays.Add(line1);
                Path line2 = new() {new IntPoint(sPoint), new IntPoint(endPoint2)};
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

            long[] resultX = new long[rays.Count];
            long[] resultY = new long[rays.Count];
            double[] weight = new double[rays.Count];

            previousEdgeNormal = new IntPoint(currentEdgeNormal.X, currentEdgeNormal.Y);

            object resultLock = new();
            Parallel.For(0, rays.Count, po_inner, ray =>
                {
                    Clipper d = new();
                    if (sideRayFallOff != falloff.none)
                    {
                        d.ZFillFunction = prox_ZFillCallback;
                    }
                    d.AddPath(rays[ray], PolyType.ptSubject, false);
                    d.AddPaths(collisionPaths, PolyType.ptClip, true);
                    PolyTree polyTree = new();
                    switch (invert)
                    {
                        case true:
                            d.Execute(ClipType.ctIntersection, polyTree);
                            break;
                        default:
                            d.Execute(ClipType.ctDifference, polyTree);
                            break;
                    }

                    // There is no matching order in the return output here, so we have to take this odd approach.
                    Paths tmpLine = Clipper.OpenPathsFromPolyTree(polyTree);

                    int tmpLineCount = tmpLine.Count;

                    switch (tmpLineCount)
                    {
                        // If we got no result, let's get that sorted out. End and start points are the same.
                        case 0:
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

                            break;
                        case > 1:
                        {
                            // We got two lines back. Need to check this carefully.
                            // We need to find the line that has a start/end point matched to our origin point for the ray.
                            int index = -1;
                            for (int tL = 0; tL < tmpLine.Count; tL++)
                            {
                                long tL0X = tmpLine[tL][0].X;
                                long tL0Y = tmpLine[tL][0].Y;
                                long tL1X = tmpLine[tL][1].X;
                                long tL1Y = tmpLine[tL][1].Y;
                                if ((tL0X != startPoint.X || tL0Y != startPoint.Y) &&
                                    (tL1X != startPoint.X || tL1Y != startPoint.Y))
                                {
                                    continue;
                                }

                                index = tL;
                                break;
                            }
                            Path tPath = new();
                            switch (index)
                            {
                                case >= 0:
                                    tPath = tmpLine[index];
                                    break;
                                default:
                                    tPath.Add(startPoint);
                                    tPath.Add(startPoint);
                                    break;
                            }
                            tmpLine.Clear();
                            tmpLine.Add(tPath);
                            break;
                        }
                    }

                    tmpLineCount = tmpLine.Count;

                    for (int tL = 0; tL < tmpLineCount; tL++)
                    {
                        // Figure out which end of the result line matches our origin point.
                        if (tmpLine[tL][0].X == startPoint.X && tmpLine[tL][0].Y == startPoint.Y)
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
                        else if (tmpLine[tL][1].X == startPoint.X && tmpLine[tL][1].Y == startPoint.Y)
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

            Path resultPath = new() {startPoint};

            int xCount = 0;
            long xAv = 0;
            int yCount = 0;
            long yAv = 0;

            switch (sideRayFallOff)
            {
                // If we are not truncating by weight, we do not need to average here - it was done with the normalization above.
                case falloff.none:
                {
                    // Average the result to give a weighted spacing across the rays.
                    for (int result = 0; result < resultX.Length; result++)
                    {
                        switch (Math.Abs(resultX[result]))
                        {
                            case > 1000:
                                xCount++;
                                xAv += resultX[result];
                                break;
                        }

                        switch (Math.Abs(resultY[result]))
                        {
                            case > 1000:
                                yCount++;
                                yAv += resultY[result];
                                break;
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

                    break;
                }
                default:
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
                        if (sideRayFallOff != falloff.none && !truncateRaysByWeight && totalWeight > 0)
                        {
                            weight_ = weight[w] / totalWeight;
                        }

                        switch (Math.Abs(resultX[w]))
                        {
                            case > 1000:
                                xCount++;
                                xAv += Convert.ToInt64(weight_ * resultX[w]);
                                break;
                        }
                        switch (Math.Abs(resultY[w]))
                        {
                            case > 1000:
                                yCount++;
                                yAv += Convert.ToInt64(weight_ * resultY[w]);
                                break;
                        }
                    }

                    break;
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