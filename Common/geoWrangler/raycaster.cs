using Clipper2Lib;
using System;
using System.Collections.Generic;
using System.Threading;
using System.Threading.Tasks;
using utility;

namespace geoWrangler;

public class RayCast
{
    public enum inversionMode { none, x, y }
    private PathsD clippedLines;
    private PathsD castLines;

    public enum forceSingleDirection { no, vertical, horizontal } // No and horizontal are treated the same in this code at the moment; leaving these to make the code more readable and intent clearer.

    private void prox_ZFillCallback(PointD bot1, PointD top1, PointD bot2, PointD top2, ref PointD pt)
    {
        pt.z = bot1.z;
    }

    public static readonly List<string> fallOffList = new() {"None", "Linear", "Gaussian", "Cosine"};
    public enum falloff { none, linear, gaussian, cosine }

    // Corner projection is default and takes an orthogonal ray out from the corner. Setting to false causes an averaged normal to be generated.

    public RayCast(PathD emissionPath, PathsD collisionPaths, long max, bool projectCorners = true, inversionMode invert = inversionMode.none, int multisampleRayCount = 0, bool runOuterLoopThreaded = false, bool runInnerLoopThreaded = false, PointD startOffset = new(), PointD endOffset = new(), falloff sideRayFallOff = falloff.none, double sideRayFallOffMultiplier = 1.0f, forceSingleDirection dirOverride = forceSingleDirection.no)
    {
        pRayCast(emissionPath, collisionPaths, max, projectCorners, invert, multisampleRayCount, runOuterLoopThreaded, runInnerLoopThreaded, startOffset, endOffset, sideRayFallOff, sideRayFallOffMultiplier, dirOverride);
    }

    public RayCast(PathD emissionPath, PathD collisionPath, long max, bool projectCorners = true, inversionMode invert = inversionMode.none, int multisampleRayCount = 0, bool runOuterLoopThreaded = false, bool runInnerLoopThreaded = false, PointD startOffset = new(), PointD endOffset = new(), falloff sideRayFallOff = falloff.none, double sideRayFallOffMultiplier = 1.0f, forceSingleDirection dirOverride = forceSingleDirection.no)
    {
        pRayCast(emissionPath, new () { collisionPath }, max, projectCorners, invert, multisampleRayCount, runOuterLoopThreaded, runInnerLoopThreaded, startOffset, endOffset, sideRayFallOff, sideRayFallOffMultiplier, dirOverride);
    }

    public PathsD getRays()
    {
        return pGetRays();
    }

    private PathsD pGetRays()
    {
        return new(castLines);
    }

    public PathsD getClippedRays()
    {
        return pGetClippedRays();
    }

    private PathsD pGetClippedRays()
    {
        return new(clippedLines);
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

    class NormalsData
    {
        public PointD[] normals;// = new Point64[ptCount];
        public PointD[] previousNormals;// = new Point64[ptCount];
    }

    private static NormalsData pCalculateNormalsData(PathD path, bool closedPathEmitter, PointD startOffset, PointD endOffset)
    {
        NormalsData ret = new();
        int ptCount = path.Count;

        // This is a serial evaluation as we need both the previous and the current normal for each point.
        ret.normals = new PointD[ptCount];
        ret.previousNormals = new PointD[ptCount];
        for (int pt = 0; pt < ptCount; pt++)
        {
            // Start point
            double dx;
            double dy;
            if (pt == path.Count - 1)
            {
                if (closedPathEmitter)
                {
                    // Last matches the first. Since we flip the dx and dy tone later, we need to compensate here.
                    dx = -ret.normals[0].x;
                    dy = -ret.normals[0].y;
                }
                else
                {
                    dx = path[ptCount - 1].x - endOffset.x;
                    dy = path[ptCount - 1].y - endOffset.y;
                }
            }
            else
            {
                if (!closedPathEmitter && pt == 0)
                {
                    dx = path[0].x - startOffset.x;
                    dy = path[0].y - startOffset.y;
                }
                else
                {
                    dx = path[pt + 1].x - path[pt].x;
                    dy = path[pt + 1].y - path[pt].y;
                }
            }

            ret.normals[pt] = new (-dx, -dy);

            // Previous normal
            if (pt == 0)
            {
                if (closedPathEmitter)
                {
                    // n-1 identical to the 0-th point, so we need to dig a little deeper.
                    dx = path[0].x - path[ptCount - 2].x;
                    dy = path[0].y - path[ptCount - 2].y;
                }
                else
                {
                    dx = path[0].x - startOffset.x;
                    dy = path[0].y - startOffset.y;
                }

                ret.previousNormals[pt] = new (-dx, -dy);
            }
            else
            {
                ret.previousNormals[pt] = new(ret.normals[pt - 1]);
            }
        }

        return ret;
    }

    // Setting this to true, we shorten rays with the falloff. False means we reduce the contribution to the average instead.
    const bool truncateRaysByWeight = false;

    private PathsD pGenerateRays(PathD sourcePath, int index, long maxRayLength, bool projectCorners, inversionMode invert, int multisampleRayCount, falloff sideRayFallOff, double sideRayFallOffMultiplier, NormalsData nData, forceSingleDirection dirOverride)
    {
        PointD startPoint = new(sourcePath[index]);
        
        PathsD rays = new();

        PointD averagedEdgeNormal = pGetAveragedNormal(nData, index, projectCorners, invert, dirOverride);

        // Normalization. We don't change the original vectors to avoid having to normalize everywhere.
        double length = Math.Sqrt(Utils.myPow(averagedEdgeNormal.x, 2) + Utils.myPow(averagedEdgeNormal.y, 2));

        double endPointDeltaX = 0;
        double endPointDeltaY = 0;

        if (length > 0.001)
        {
            // Avoid div-by-zero; 0-length is unimportant. Note that setting this cut-off too high produces artifacts.
            endPointDeltaX = Convert.ToDouble(averagedEdgeNormal.x) / length;
            endPointDeltaY = Convert.ToDouble(averagedEdgeNormal.y) / length;
        }

        // Set to max ray length from callsite.
        endPointDeltaX *= maxRayLength;
        endPointDeltaY *= maxRayLength;

        switch (invert)
        {
            case inversionMode.x:
            case inversionMode.y:
                endPointDeltaY *= -1;
                break;
            case inversionMode.none:
            default:
                break;
        }
        endPointDeltaX *= -1;

        PointD endPoint = new(endPointDeltaY + startPoint.x, endPointDeltaX + startPoint.y);

        if (sideRayFallOff != falloff.none)
        {
            endPoint.z = (long)1E4;
        }
        
        rays.Add(new() {new (startPoint), new (endPoint)});

        double angleStep = 90.0f / (1 + multisampleRayCount);

        for (int sample = 0; sample < multisampleRayCount; sample++)
        {
            // Add more samples, each n-degrees rotated from the nominal ray
            double rayAngle = (sample + 1) * angleStep;

            PointD endPoint_f = new(endPoint);

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
                case falloff.none:
                default:
                    break;
            }

            if (truncateRaysByWeight)
            {
                endPoint_f = new (startPoint.x + weight_val * endPointDeltaY,
                                         startPoint.y + weight_val * endPointDeltaX);
            }

            PointD sPoint = new(startPoint.x, startPoint.y);

            if (sideRayFallOff != falloff.none)
            {
                endPoint_f.z = Convert.ToInt64(weight_val * 1E4);
                sPoint.z = endPoint_f.z;
            }
            PointD endPoint1 = GeoWrangler.Rotate(startPoint, endPoint_f, rayAngle);
            PointD endPoint2 = GeoWrangler.Rotate(startPoint, endPoint_f, -rayAngle);

            // The order of line1 below is important, but I'm not yet sure why. If you change it, the expansion becomes asymmetrical on a square (lower section gets squashed).
            rays.Add(new() {new (endPoint1), new (sPoint)});
            rays.Add(new() {new (sPoint), new (endPoint2)});
        }

        return rays;
    }

    private PointD pGetAveragedNormal(NormalsData nData, int index, bool projectCorners, inversionMode invert, forceSingleDirection dirOverride)
    {
        PointD averagedEdgeNormal;
        PointD currentEdgeNormal = new(nData.normals[index]);
        PointD previousEdgeNormal = new(nData.previousNormals[index]);

        // Get average angle for this vertex based on angles from line segments.
        // http://stackoverflow.com/questions/1243614/how-do-i-calculate-the-normal-vector-of-a-line-segment

        if (projectCorners && (Math.Abs(currentEdgeNormal.x) < constants.tolerance && Math.Abs(previousEdgeNormal.y) < constants.tolerance ||
                       Math.Abs(currentEdgeNormal.y) < constants.tolerance && Math.Abs(previousEdgeNormal.x) < constants.tolerance))
        {
            double tX = currentEdgeNormal.x;
            double tY = currentEdgeNormal.y;
            // If we're traversing a 90 degree corner, let's not project a diagonal, but fix on our current edge normal.
            if (invert == 0 || dirOverride == forceSingleDirection.vertical)
            {
                tX = -tX;
                tY = -tY;
            }
            averagedEdgeNormal = new (tX, tY);
        }
        else
        {
            switch (invert)
            {
                case inversionMode.x:
                    currentEdgeNormal = new (-currentEdgeNormal.x, currentEdgeNormal.y);
                    previousEdgeNormal = new (-previousEdgeNormal.x, previousEdgeNormal.y);
                    break;
                case inversionMode.y:
                    currentEdgeNormal = new (currentEdgeNormal.x, -currentEdgeNormal.y);
                    previousEdgeNormal = new (previousEdgeNormal.x, -previousEdgeNormal.y);
                    break;
                case inversionMode.none:
                default:
                    break;
            }
            // Average out our normals
            averagedEdgeNormal = new ((previousEdgeNormal.x + currentEdgeNormal.x) / 2, (previousEdgeNormal.y + currentEdgeNormal.y) / 2);
        }

        return averagedEdgeNormal;
    }

    private PathsD pCutRay(PathD ray, PathsD collisionPaths, inversionMode invert, falloff sideRayFallOff)
    {
        ClipperD d = new();
        if (sideRayFallOff != falloff.none)
        {
            d.ZCallback = prox_ZFillCallback;
        }
        d.AddOpenSubject(ray);
        d.AddClip(collisionPaths);
        PathsD unused = new();
        PathsD tmpLine = new();
        switch (invert)
        {
            default:
                try
                {
                    d.Execute(ClipType.Intersection, FillRule.EvenOdd, unused, tmpLine);
                }
                catch (Exception e)
                {
                    int yyy = 2;
                }
                break;
            case 0:
                d.Execute(ClipType.Difference, FillRule.EvenOdd, unused, tmpLine);
                break;
        }

        return tmpLine;
    }

    readonly object resultLock = new();
    
    private void pEvaluateCutRay(PathsD ray, int outputIndex, PathD emissionPath, int pt, ref double[] resultX, ref double[] resultY, ref double[] weight)
    {
        PointD startPoint = new(emissionPath[pt]);
        int rayPtCount = ray.Count;

        // There is no matching order in the ray here, so we have to take this odd approach.
        switch (rayPtCount)
        {
            // If we got no result, let's get that sorted out. End and start points are the same.
            case 0:
                Monitor.Enter(resultLock);
                try
                {
                    resultX[outputIndex] = startPoint.x;
                    resultY[outputIndex] = startPoint.y;
                    weight[outputIndex] = 1.0f;
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
                for (int tL = 0; tL < ray.Count; tL++)
                {
                    double tL0X = ray[tL][0].x;
                    double tL0Y = ray[tL][0].y;
                    double tL1X = ray[tL][1].x;
                    double tL1Y = ray[tL][1].y;
                    if ((Math.Abs(tL0X - startPoint.x) > constants.tolerance || Math.Abs(tL0Y - startPoint.y) > constants.tolerance) &&
                        (Math.Abs(tL1X - startPoint.x) > constants.tolerance || Math.Abs(tL1Y - startPoint.y) > constants.tolerance))
                    {
                        continue;
                    }

                    index = tL;
                    break;
                }
                PathD tPath = new();
                switch (index)
                {
                    case >= 0:
                        tPath = new(ray[index]);
                        break;
                    default:
                        tPath.Add(startPoint);
                        tPath.Add(startPoint);
                        break;
                }
                ray.Clear();
                ray.Add(tPath);
                break;
            }
        }
        
        rayPtCount = ray.Count;

        for (int tL = 0; tL < rayPtCount; tL++)
        {
            // Figure out which end of the result line matches our origin point.
            if (ray[tL][0].z == startPoint.z)
            {
                Monitor.Enter(resultLock);
                try
                {
                    resultX[outputIndex] = ray[tL][1].x;
                    resultY[outputIndex] = ray[tL][1].y;
                    weight[outputIndex] = 1E4;
                }
                finally
                {
                    Monitor.Exit(resultLock);
                }
            }
            else if (ray[tL][1].z == startPoint.z)
            {
                Monitor.Enter(resultLock);
                try
                {
                    // Clipper reversed the line direction, so we need to deal with this.
                    resultX[outputIndex] = ray[tL][0].x;
                    resultY[outputIndex] = ray[tL][0].y;
                    weight[outputIndex] = 1E4;
                }
                finally
                {
                    Monitor.Exit(resultLock);
                }
            }
        }
    }

    class ResultData
    {
        public double xAv;
        public double yAv;
    }
    private ResultData pComputeWeightedResult(falloff sideRayFallOff, ref double[] resultX, ref double[] resultY, ref double[] weight)
    {
        int xCount = 0;
        int yCount = 0;

        ResultData res = new();
        
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
                        // Ignore values that are below our floating point tolerance.
                        case > constants.tolerance:
                            xCount++;
                            res.xAv += resultX[result];
                            break;
                    }

                    switch (Math.Abs(resultY[result]))
                    {
                        // Ignore values that are below our floating point tolerance.
                        case > constants.tolerance:
                            yCount++;
                            res.yAv += resultY[result];
                            break;
                    }
                }

                if (xCount != 0)
                {
                    res.xAv /= xCount;
                }

                if (yCount != 0)
                {
                    res.yAv /= yCount;
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
                        // Ignore values that are below our floating point tolerance.
                        case > constants.tolerance:
                            xCount++;
                            res.xAv += weight_ * resultX[w];
                            break;
                    }

                    switch (Math.Abs(resultY[w]))
                    {
                        // Ignore values that are below our floating point tolerance.
                        case > constants.tolerance:
                            yCount++;
                            res.yAv += weight_ * resultY[w];
                            break;
                    }
                }

                break;
            }
        }

        return res;
    }

    // invert used to be a bool, but we need to handle X and Y normal inversions separately, so this had to move to an enum for clarity.
    private void pRayCast(PathD emissionPath, PathsD collisionPaths, long maxRayLength, bool projectCorners, inversionMode invert, int multisampleRayCount, bool runOuterLoopThreaded, bool runInnerLoopThreaded, PointD startOffset, PointD endOffset, falloff sideRayFallOff, double sideRayFallOffMultiplier, forceSingleDirection dirOverride)
    {
        int ptCount = emissionPath.Count;

        // Due to threading and need to tie to polygon point order, we have to use these local storage options and will do the conversion at the end.
        object castLinesLock = new();
        PathsD[] castLines_ = new PathsD[ptCount];
        object clippedLinesLock = new();
        PathD[] clippedLines_ = new PathD[ptCount];

        // We need to think about the end point case.
        bool closedPathEmitter = ptCount > 3 && Math.Abs(emissionPath[0].x - emissionPath[ptCount - 1].x) < constants.tolerance && Math.Abs(emissionPath[0].y - emissionPath[ptCount - 1].y) < constants.tolerance;

        // Get average angle for this vertex based on angles from line segments.
        // http://stackoverflow.com/questions/1243614/how-do-i-calculate-the-normal-vector-of-a-line-segment

        // Pre-calculate these for the threading to be an option.
        NormalsData nData = pCalculateNormalsData(emissionPath, closedPathEmitter, startOffset, endOffset);

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
            PointD startPoint = new(emissionPath[pt]);
            
            PathsD rays = pGenerateRays(emissionPath, pt, maxRayLength, projectCorners, invert, multisampleRayCount, sideRayFallOff, sideRayFallOffMultiplier, nData, dirOverride);
            
            Monitor.Enter(castLinesLock);
            try
            {
                castLines_[pt] = new (rays);
            }
            finally
            {
                Monitor.Exit(castLinesLock);
            }

            double[] resultX = new double[rays.Count];
            double[] resultY = new double[rays.Count];
            double[] weight = new double[rays.Count];
            
            Parallel.For(0, rays.Count, po_inner, ray =>
                {

                    PathsD tmpLine = pCutRay(rays[ray], collisionPaths, invert, sideRayFallOff);
                    pEvaluateCutRay(tmpLine, ray, emissionPath, pt, ref resultX, ref resultY, ref weight);
                }
            );

            PathD resultPath = new() {new(startPoint)};

            // Note that results below constants.tolerance will be ignored in the below call.
            ResultData rData = pComputeWeightedResult(sideRayFallOff, ref resultX, ref resultY, ref weight);

            resultPath.Add(new (rData.xAv, rData.yAv));
            Monitor.Enter(clippedLinesLock);
            try
            {
                clippedLines_[pt] = new(resultPath);
            }
            finally
            {
                Monitor.Exit(clippedLinesLock);
            }
        });

        // Convert the array back to a list.
        clippedLines = new(clippedLines_);
        castLines = new ();
        foreach (PathsD t in castLines_)
        {
            castLines.AddRange(t);
        }
    }
}