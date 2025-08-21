﻿using System;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Collections.Generic;

using PathD = System.Collections.Generic.List<PointD>;
using PathsD = System.Collections.Generic.List<System.Collections.Generic.List<PointD>>;

static class QuadraticBezierSamplingSwitcher_Polygon
{
    /*
     * Added: a QuinticC2 easing strategy that builds a sampled quintic-Hermite
     * (C2) replacement for the start/end blended segments around short-edge diagonals.
     *
     * Implementation notes:
     * - The new EasingStrategy.QuinticC2 option attempts to produce a C2-like
     *   local replacement by estimating second-derivative vectors from the
     *   Hermite first-derivatives of the segment endpoints.
     * - The second-derivative estimates are conservative (scaled difference of the
     *   first derivatives) to avoid overshoot. This is primarily for interactive
     *   testing — you can tweak scaling if you want stronger or weaker curvature.
     *
     * How to test:
     * - In Main() set easingMode = EasingStrategy.QuinticC2 and run to emit assembled.svg.
     * - Tune the smoothing parameters (insetFraction, neighborhood/sample resolution).
     */

    // -------------------------------------------------------------------------
    // Enums and small configuration-like types
    // -------------------------------------------------------------------------
    enum SamplingMode
    {
        ByMaxSegmentLength,
        ByMaxAngle
    }

    private enum CornerType
    {
        Concave,
        Convex,
        ShortEdge
    }

    private enum SmoothStepMode
    {
        None,
        Cubic,
        Quintic
    }

    enum EasingStrategy
    {
        None,
        CubicBezierHermite, // cubic bezier constructed from Hermite derivatives (no overshoot)
        QuinticHermite,      // legacy quintic hermite (C2) placeholder (kept)
        SmoothBlend,         // quintic Hermite sampled (soft)
        CircularArc,         // tangent arc with cubic fallback
        QuinticC2            // NEW: use quintic Hermite with estimated second derivatives (C2-like)
    }

    // -------------------------------------------------------------------------
    // Example/demo entry point
    // -------------------------------------------------------------------------
    static void Main()
    {
        // Example polygon (closed - first == last)
        List<PointD> original_path = new List<PointD>
        {
            new PointD(0, 0),
            new PointD(0, 120),
            new PointD(20, 120),
            new PointD(20, 100),
            new PointD(40, 100),
            new PointD(40, 80),
            new PointD(59, 80),
            new PointD(59, 59),
            new PointD(110, 59),
            new PointD(110, 0),
            new PointD(0, 0),
        };

        // Parameters that control curvature, resolution and handling of short edges
        double concaveRadius = 20.0;
        double convexRadius = 50.0;
        double edgeResolution = 0.05;   // smaller => more points on curves
        double angularResolution = 1.0; // degrees
        double shortEdgeLength = 20.0;
        double maxShortEdgeLength = 30.0;

        // smoothing to apply to the blend factor across the short-edge run
        SmoothStepMode smoothMode = SmoothStepMode.Quintic;

        // 1) classify corners
        int[] corner_types = CategorizeCorners(original_path, shortEdgeLength);

        // 2) build sampled corner polylines and record midpoints used for diagonals
        PathsD processed = new PathsD();
        List<PointD> cornerMidpoints = new List<PointD>(); // midpoint on outgoing edge for each corner
        List<PointD> cornerVertices = new List<PointD>();  // original polygon vertex positions

        for (int i = 0; i < original_path.Count - 1; i++)
        {
            cornerVertices.Add(original_path[i]);

            // start-line: vertex -> midpoint of incoming edge
            PointD prevMid = (i == 0)
                ? Mid(original_path[i], original_path[^2])
                : Mid(original_path[i], original_path[i - 1]);

            var startLine = new PathD { original_path[i], prevMid };

            // end-line: vertex -> midpoint of outgoing edge
            PointD nextMid = (i == original_path.Count - 2)
                ? Mid(original_path[i], original_path[0])
                : Mid(original_path[i], original_path[i + 1]);

            var endLine = new PathD { original_path[i], nextMid };

            // store midpoint on the outgoing edge for index i
            cornerMidpoints.Add(nextMid);

            // choose radius depending on concavity/convexity
            double radius = convexRadius;
            if (corner_types[i] == (int)CornerType.Concave)
                radius = concaveRadius;

            // produce sampled corner polyline
            PathD current_corner = ProcessCorner(startLine, endLine, radius, angularResolution, edgeResolution,
                SamplingMode.ByMaxSegmentLength);
            processed.Add(current_corner);
        }

        // Easing configuration for assembling diagonals between short-edge runs and neighbours
        // Try the new option here:
        EasingStrategy easingMode = EasingStrategy.QuinticC2; // <-- switch to QuinticC2 to test
        double insetFraction = 0.12; // relative to diagonal length (0..0.5)
        double minInset = 2.0;
        double maxInset = 40.0;
        int diagStraightSample = 10; // samples for straight center section (0 = just endpoints)

        // 3) assemble final polyline with easing across short-edge runs
        PathD assembled = AssembleWithEasing(processed, corner_types, easingMode, insetFraction, minInset, maxInset,
            diagStraightSample, cornerMidpoints, cornerVertices, shortEdgeLength, maxShortEdgeLength, smoothMode,
            edgeResolution);

        // 4) optional final smoothing pass to ease any remaining seams between diagonals and curves
        PathD smoothed = SmoothPolylineSeams(assembled, angleThresholdDeg: 8.0, neighborhoodRadius: 3, sampleLenHint: edgeResolution);

        // 5) write an SVG for inspection
        string svg = BuildDetailedSvg(original_path, smoothed);
        File.WriteAllText("assembled.svg", svg, Encoding.UTF8);
        Console.WriteLine("Wrote assembled.svg (easingMode = " + easingMode + ", smoothMode = " + smoothMode + ")");
    }

    // -------------------------------------------------------------------------
    // Corner classification
    // -------------------------------------------------------------------------
    static int[] CategorizeCorners(PathD path_, double shortEdgeLength)
    {
        PathD path = new PathD(path_);
        bool trimmed_path = false;

        if (path.Count > 1 && path[0].x == path[^1].x && path[0].y == path[^1].y)
        {
            trimmed_path = true;
            path.RemoveAt(path.Count - 1);
        }

        double area2 = 0;
        for (int i = 0; i < path.Count; i++)
        {
            PointD p1 = path[i];
            PointD p2 = path[(i + 1) % path.Count];
            area2 += p1.x * p2.y - p2.x * p1.y;
        }

        bool isCCW = area2 > 0;
        int[] status = new int[path.Count];

        for (int i = 0; i < path.Count; i++)
        {
            PointD prev = path[(i - 1 + path.Count) % path.Count];
            PointD curr = path[i];
            PointD next = path[(i + 1) % path.Count];

            PointD v1 = Minus(curr, prev);
            PointD v2 = Minus(next, curr);

            double len1 = Length(v1);
            double len2 = Length(v2);

            if (len1 <= shortEdgeLength && len2 <= shortEdgeLength)
            {
                status[i] = (int)CornerType.ShortEdge;
                continue;
            }

            double crossZ = v1.x * v2.y - v1.y * v2.x;
            bool isVertexConvex = isCCW ? (crossZ > 0) : (crossZ < 0);
            status[i] = (int)(isVertexConvex ? CornerType.Convex : CornerType.Concave);
        }

        if (trimmed_path)
        {
            Array.Resize(ref status, status.Length + 1);
            status[^1] = status[0];
        }

        return status;
    }

    // -------------------------------------------------------------------------
    // Corner processing
    // -------------------------------------------------------------------------
    static PathD ProcessCorner(PathD startLine, PathD endLine, double radius, double angular_resolution,
        double edge_resolution, SamplingMode mode = SamplingMode.ByMaxAngle)
    {
        PointD startVertex = startLine[0];
        PointD startMid = startLine[1];
        PointD endVertex = endLine[0];
        PointD endMid = endLine[1];

        PointD startDir = Normalized(Minus(startMid, startVertex));
        PointD endDir = Normalized(Minus(endMid, endVertex));

        double startEdgeLen = Length(Minus(startMid, startVertex));
        double startRadius = Math.Min(radius, Math.Abs(startEdgeLen));

        double endEdgeLen = Length(Minus(endMid, endVertex));
        double endRadius = Math.Min(radius, Math.Abs(endEdgeLen));

        PointD curveStartPoint = Add(startVertex, Mul(startDir, startRadius));
        PointD curveEndPoint = Add(endVertex, Mul(endDir, endRadius));

        double det = startDir.x * endDir.y - startDir.y * endDir.x;
        if (Math.Abs(det) < 1e-9)
        {
            var straight = new PathD { curveStartPoint, curveEndPoint };
            return straight;
        }

        double dx = curveEndPoint.x - curveStartPoint.x;
        double dy = curveEndPoint.y - curveStartPoint.y;
        double tParam = (dx * endDir.y - dy * endDir.x) / det;
        PointD controlPoint = Add(curveStartPoint, Mul(startDir, tParam));

        PathD samples = mode switch
        {
            SamplingMode.ByMaxSegmentLength => SampleByMaxSegmentLength(curveStartPoint, controlPoint, curveEndPoint, edge_resolution),
            SamplingMode.ByMaxAngle => SampleByMaxAngle(curveStartPoint, controlPoint, curveEndPoint, angular_resolution * Math.PI / 180.0),
            _ => throw new ArgumentOutOfRangeException(nameof(mode))
        };

        return samples;
    }

    // -------------------------------------------------------------------------
    // Assembly with easing (uses the new QuinticC2 option)
    // -------------------------------------------------------------------------
    static PathD AssembleWithEasing(
        PathsD processedCorners,
        int[] corner_types,
        EasingStrategy strategy,
        double insetFraction,
        double minInset,
        double maxInset,
        int diagStraightSample,
        List<PointD> cornerMidpoints,
        List<PointD> cornerVertices,
        double short_edge_length,
        double max_short_edge_length,
        SmoothStepMode smoothMode,
        double edgeResolution)
    {
        int n = processedCorners.Count;
        if (cornerMidpoints == null || cornerMidpoints.Count != n)
            throw new ArgumentException("cornerMidpoints must be same length as processedCorners");
        if (cornerVertices == null || cornerVertices.Count != n)
            throw new ArgumentException("cornerVertices must be same length as processedCorners");

        PathD outPts = new PathD();

        static double CubicSmoothstep(double t) => t <= 0 ? 0 : t >= 1 ? 1 : (3.0 * t * t - 2.0 * t * t * t);
        static double QuinticSmoothstep(double t) =>
            t <= 0 ? 0 : t >= 1 ? 1 : (6.0 * Math.Pow(t, 5) - 15.0 * Math.Pow(t, 4) + 10.0 * Math.Pow(t, 3));
        static double ApplySmooth(SmoothStepMode mode, double t) =>
            mode switch
            {
                SmoothStepMode.Cubic => CubicSmoothstep(t),
                SmoothStepMode.Quintic => QuinticSmoothstep(t),
                _ => Math.Clamp(t, 0.0, 1.0)
            };

        bool[] isShort = new bool[n];
        for (int i = 0; i < n; i++)
            isShort[i] = (i < corner_types.Length && corner_types[i] == (int)CornerType.ShortEdge);

        int idx = 0;
        while (idx < n)
        {
            if (!isShort[idx])
            {
                PathD poly = processedCorners[idx];
                if (poly.Count > 0)
                {
                    if (outPts.Count == 0) outPts.AddRange(poly);
                    else
                    {
                        if (PointsEqual(outPts[^1], poly[0]))
                            outPts.AddRange(poly.Skip(1));
                        else
                            outPts.AddRange(poly);
                    }
                }
                idx++;
                continue;
            }

            int runStart = idx;
            while (idx < n && isShort[idx]) idx++;
            int runEnd = idx - 1;

            int prevIdx = (runStart - 1 + n) % n;
            while (isShort[prevIdx]) prevIdx = (prevIdx - 1 + n) % n;

            int nextIdx = (runEnd + 1) % n;
            while (isShort[nextIdx]) nextIdx = (nextIdx + 1) % n;

            PointD processedStartPt = processedCorners[prevIdx].Last();
            PointD processedEndPt = processedCorners[nextIdx].First();

            PointD diagStartMid = cornerMidpoints[prevIdx];
            PointD diagEndMid = cornerMidpoints[nextIdx];
            PointD vertexPrev = cornerVertices[prevIdx];
            PointD vertexNext = cornerVertices[nextIdx];

            double dPrev = Length(Minus(diagStartMid, vertexPrev));
            double dNext = Length(Minus(diagEndMid, vertexNext));

            double denom = (max_short_edge_length - short_edge_length);
            double tPrev = denom == 0 ? 0.0 : Math.Clamp((dPrev - short_edge_length) / denom, 0.0, 1.0);
            double tNext = denom == 0 ? 0.0 : Math.Clamp((dNext - short_edge_length) / denom, 0.0, 1.0);

            double blendPrev = ApplySmooth(smoothMode, tPrev);
            double blendNext = ApplySmooth(smoothMode, tNext);

            if (blendPrev <= 1e-9 && blendNext <= 1e-9)
            {
                if (outPts.Count == 0) outPts.Add(processedStartPt);
                else if (!PointsEqual(outPts[^1], processedStartPt)) outPts.Add(processedStartPt);

                if (diagStraightSample <= 0)
                {
                    outPts.Add(processedEndPt);
                }
                else
                {
                    PointD diagFull = Minus(processedEndPt, processedStartPt);
                    for (int s = 1; s <= diagStraightSample; s++)
                    {
                        double t = (double)s / (diagStraightSample + 1);
                        PointD p = new PointD(processedStartPt.x + diagFull.x * t, processedStartPt.y + diagFull.y * t);
                        outPts.Add(p);
                    }
                    outPts.Add(processedEndPt);
                }

                continue;
            }

            PointD diag = Minus(processedEndPt, processedStartPt);
            double diagLen = Length(diag);
            if (diagLen <= 1e-9)
            {
                if (outPts.Count == 0 || !PointsEqual(outPts[^1], processedEndPt)) outPts.Add(processedEndPt);
                continue;
            }

            PointD diagDir = Normalized(diag);

            double inset = Math.Max(minInset, Math.Min(maxInset, diagLen * insetFraction));
            if (inset * 2.0 > diagLen) inset = diagLen * 0.5 * 0.999;

            PointD S = Add(processedStartPt, Mul(diagDir, inset));
            PointD E = Minus(processedEndPt, Mul(diagDir, inset));

            PointD prevTangent = EstimateOutgoingTangent(processedCorners[prevIdx]);
            PointD nextTangent = EstimateIncomingTangent(processedCorners[nextIdx]);

            if (Length(prevTangent) < 1e-12) prevTangent = diagDir;
            if (Length(nextTangent) < 1e-12) nextTangent = diagDir;

            PointD prevDirUnit = Normalized(prevTangent);
            PointD nextDirUnit = Normalized(nextTangent);

            double tPrevBlend = ApplySmooth(smoothMode, blendPrev);
            double tNextBlend = ApplySmooth(smoothMode, blendNext);

            PointD blendedPrevDir = Normalized(Add(Mul(prevDirUnit, 1 - tPrevBlend), Mul(diagDir, tPrevBlend)));
            PointD blendedNextDir = Normalized(Add(Mul(nextDirUnit, 1 - tNextBlend), Mul(diagDir, tNextBlend)));

            double Lstart = Math.Max(1e-9, Length(Minus(S, processedStartPt)));
            double Lend = Math.Max(1e-9, Length(Minus(processedEndPt, E)));

            PointD H_blended_start = Mul(blendedPrevDir, Lstart);
            PointD H_blended_end = Mul(blendedNextDir, Lend);

            double magStart = Length(H_blended_start);
            double magEnd = Length(H_blended_end);

            PointD finalDirAtS = Normalized(Add(Mul(blendedPrevDir, 1 - tPrevBlend), Mul(diagDir, tPrevBlend)));
            PointD finalDirAtE = Normalized(Add(Mul(blendedNextDir, 1 - tNextBlend), Mul(Neg(diagDir), tNextBlend)));

            PointD finalDerivAtS = Mul(finalDirAtS, magStart);
            PointD finalDerivAtE = Mul(finalDirAtE, magEnd);

            PointD T0_for_start = H_blended_start;
            PointD T1_for_start = finalDerivAtS;

            PointD T0_for_end = finalDerivAtE;
            PointD T1_for_end = H_blended_end;

            if (strategy == EasingStrategy.None)
            {
                if (outPts.Count == 0) outPts.Add(processedStartPt);
                else if (!PointsEqual(outPts[^1], processedStartPt)) outPts.Add(processedStartPt);

                PointD diagControl = new PointD((processedStartPt.x + processedEndPt.x) / 2.0,
                    (processedStartPt.y + processedEndPt.y) / 2.0);
                PathD blendedCurve = SampleByMaxSegmentLength(processedStartPt, diagControl, processedEndPt, 0.5);
                outPts.AddRange(blendedCurve);
                continue;
            }

            PathD startSeg = null;
            PathD endSeg = null;

            switch (strategy)
            {
                case EasingStrategy.CubicBezierHermite:
                    startSeg = BuildCubicHermiteSampled(processedStartPt, S, T0_for_start, T1_for_start, edgeResolution);
                    endSeg = BuildCubicHermiteSampled(E, processedEndPt, T0_for_end, T1_for_end, edgeResolution);
                    break;

                case EasingStrategy.QuinticHermite:
                case EasingStrategy.SmoothBlend:
                    // fallback to cubic hermite sampling for now
                    startSeg = BuildCubicHermiteSampled(processedStartPt, S, T0_for_start, T1_for_start, edgeResolution);
                    endSeg = BuildCubicHermiteSampled(E, processedEndPt, T0_for_end, T1_for_end, edgeResolution);
                    break;

                case EasingStrategy.CircularArc:
                    startSeg = BuildCircularArcOrNull(processedStartPt, S, T0_for_start, diagDir, inset)
                               ?? BuildCubicHermiteSampled(processedStartPt, S, T0_for_start, T1_for_start, edgeResolution);

                    endSeg = BuildCircularArcOrNull(E, processedEndPt, Neg(diagDir), T1_for_end, inset)
                             ?? BuildCubicHermiteSampled(E, processedEndPt, T0_for_end, T1_for_end, edgeResolution);
                    break;

                case EasingStrategy.QuinticC2:
                    // Build quintic C2-like segments using estimated second derivatives.
                    // We estimate second-derivative vectors from the difference of Hermite first-derivatives,
                    // scaled conservatively to avoid overshoot.
                    startSeg = BuildQuinticC2Segment(processedStartPt, S, T0_for_start, T1_for_start, inset);
                    endSeg = BuildQuinticC2Segment(E, processedEndPt, T0_for_end, T1_for_end, inset);
                    break;

                default:
                    startSeg = BuildCubicHermiteSampled(processedStartPt, S, T0_for_start, T1_for_start, edgeResolution);
                    endSeg = BuildCubicHermiteSampled(E, processedEndPt, T0_for_end, T1_for_end, edgeResolution);
                    break;
            }

            if (startSeg != null && startSeg.Count > 0)
            {
                if (outPts.Count > 0 && PointsEqual(outPts[^1], startSeg[0]))
                    outPts.AddRange(startSeg.Skip(1));
                else
                    outPts.AddRange(startSeg);
            }
            else
            {
                if (outPts.Count == 0 || !PointsEqual(outPts[^1], S)) outPts.Add(S);
            }

            if (diagStraightSample <= 0)
            {
                if (!PointsEqual(outPts[^1], E)) outPts.Add(E);
            }
            else
            {
                for (int s = 1; s <= diagStraightSample; s++)
                {
                    double t = (double)s / (diagStraightSample + 1);
                    PointD p = Add(S, Mul(diagDir, (diagLen - 2 * inset) * t));
                    outPts.Add(p);
                }
                outPts.Add(E);
            }

            if (endSeg != null && endSeg.Count > 0)
            {
                if (PointsEqual(outPts[^1], endSeg[0])) outPts.AddRange(endSeg.Skip(1));
                else outPts.AddRange(endSeg);
            }
            else
            {
                if (!PointsEqual(outPts[^1], processedEndPt)) outPts.Add(processedEndPt);
            }
        }

        if (outPts.Count > 0)
        {
            if (Length(Minus(outPts[0], outPts[^1])) < 1e-9)
                outPts[^1] = outPts[0];
            else
                outPts.Add(outPts[0]);
        }

        return outPts;
    }

    // -------------------------------------------------------------------------
    // Quintic C2 helper: estimate second derivatives and sample quintic hermite
    // -------------------------------------------------------------------------
    static PathD BuildQuinticC2Segment(PointD a, PointD b, PointD m0, PointD m1, double inset)
    {
        // Conservative scale to avoid overshoot for second derivatives
        double secondDerivScale = 0.5;

        // Estimate second derivative vectors as scaled difference between endpoint first derivatives
        PointD diff = Minus(m1, m0);
        // normalize by segment length to keep scale reasonable
        double segLen = Math.Max(1e-9, Length(Minus(b, a)));
        PointD s0 = Mul(diff, (secondDerivScale / segLen));
        PointD s1 = Mul(Neg(diff), (secondDerivScale / segLen));

        // Ensure magnitudes are not huge by clamping
        double maxMag = inset * 0.5;
        if (Length(s0) > maxMag) s0 = Mul(Normalized(s0), maxMag);
        if (Length(s1) > maxMag) s1 = Mul(Normalized(s1), maxMag);

        int samples = SampleCountForInset(inset);
        return SampleQuinticHermite(a, b, m0, m1, s0, s1, samples);
    }

    // -------------------------------------------------------------------------
    // Final smoothing pass
    // -------------------------------------------------------------------------
    static PathD SmoothPolylineSeams(PathD pts, double angleThresholdDeg = 8.0, int neighborhoodRadius = 3, double sampleLenHint = 0.5)
    {
        if (pts == null || pts.Count < 5) return pts;

        bool closed = PointsEqual(pts[0], pts[^1]);
        PathD work = new PathD(pts);
        if (closed) work.RemoveAt(work.Count - 1);

        int n = work.Count;
        if (n < 4)
        {
            if (closed) work.Add(work[0]);
            return work;
        }

        double threshRad = Math.Abs(angleThresholdDeg) * Math.PI / 180.0;
        var ranges = new List<(int s, int e)>();

        for (int i = 0; i < n; i++)
        {
            int im1 = (i - 1 + n) % n;
            int ip1 = (i + 1) % n;

            PointD vprev = Minus(work[i], work[im1]);
            PointD vnext = Minus(work[ip1], work[i]);
            double lprev = Length(vprev);
            double lnext = Length(vnext);
            if (lprev < 1e-9 || lnext < 1e-9) continue;

            double dot = Math.Max(-1.0, Math.Min(1.0, Dot(Normalized(vprev), Normalized(vnext))));
            double ang = Math.Acos(dot);

            if (Math.Abs(ang) > threshRad)
            {
                int s = i - neighborhoodRadius;
                int e = i + neighborhoodRadius;
                s = Math.Max(0, s);
                e = Math.Min(n - 1, e);

                if (ranges.Count > 0 && s <= ranges[^1].e)
                {
                    var last = ranges[^1];
                    ranges[^1] = (last.s, Math.Max(last.e, e));
                }
                else
                {
                    ranges.Add((s, e));
                }
            }
        }

        if (ranges.Count == 0)
        {
            if (closed) work.Add(work[0]);
            return work;
        }

        for (int ri = ranges.Count - 1; ri >= 0; ri--)
        {
            int s = ranges[ri].s;
            int e = ranges[ri].e;
            if (e - s < 2) continue;

            PointD A = work[s];
            PointD B = work[e];

            PointD neighborA = work[Math.Min(s + 1, n - 1)];
            PointD neighborB = work[Math.Max(e - 1, 0)];

            PointD dirA = Normalized(Minus(neighborA, A));
            PointD dirB = Normalized(Minus(B, neighborB));
            if (Length(dirA) < 1e-12) dirA = Normalized(Minus(B, A));
            if (Length(dirB) < 1e-12) dirB = Normalized(Minus(B, A));

            double magA = Length(Minus(neighborA, A));
            double magB = Length(Minus(B, neighborB));
            magA = Math.Max(magA, Length(Minus(B, A)) * 0.1);
            magB = Math.Max(magB, Length(Minus(B, A)) * 0.1);

            PointD T0 = Mul(dirA, magA);
            PointD T1 = Mul(dirB, magB);

            PathD replacement = BuildCubicHermiteSampled(A, B, T0, T1, sampleLenHint);

            var newWork = new PathD();
            for (int k = 0; k <= s; k++) newWork.Add(work[k]);
            if (replacement.Count >= 2)
            {
                for (int k = 1; k < replacement.Count - 1; k++) newWork.Add(replacement[k]);
            }
            for (int k = e; k < n; k++) newWork.Add(work[k]);

            work = newWork;
            n = work.Count;
        }

        if (closed && work.Count > 0) work.Add(work[0]);
        return work;
    }

    // -------------------------------------------------------------------------
    // Sampling & Builder helper functions
    // -------------------------------------------------------------------------
    static PathD BuildCubicHermiteSampled(PointD A, PointD B, PointD T0, PointD T1, double maxSegLen)
    {
        HermiteToBezier(A, B, T0, T1, out PointD b0, out PointD b1, out PointD b2, out PointD b3);
        return SampleBezierByMaxSegmentLength(b0, b1, b2, b3, maxSegLen);
    }

    static void HermiteToBezier(PointD P0, PointD P1, PointD T0, PointD T1,
        out PointD B0, out PointD B1, out PointD B2, out PointD B3)
    {
        B0 = P0;
        B1 = Add(P0, Mul(T0, 1.0 / 3.0));
        B2 = Minus(P1, Mul(T1, 1.0 / 3.0));
        B3 = P1;
    }

    static PathD SampleBezierByMaxSegmentLength(PointD b0, PointD b1, PointD b2, PointD b3, double maxLen)
    {
        PathD outp = new PathD { b0 };

        void Recurse(PointD p0, PointD p1, PointD p2, PointD p3)
        {
            double chord = Length(Minus(p3, p0));
            double contLen = Length(Minus(p1, p0)) + Length(Minus(p2, p1)) + Length(Minus(p3, p2));

            if (contLen - chord <= maxLen || chord <= maxLen)
            {
                outp.Add(p3);
                return;
            }

            PointD p01 = Mul(Add(p0, p1), 0.5);
            PointD p12 = Mul(Add(p1, p2), 0.5);
            PointD p23 = Mul(Add(p2, p3), 0.5);
            PointD p012 = Mul(Add(p01, p12), 0.5);
            PointD p123 = Mul(Add(p12, p23), 0.5);
            PointD p0123 = Mul(Add(p012, p123), 0.5);

            Recurse(p0, p01, p012, p0123);
            Recurse(p0123, p123, p23, p3);
        }

        Recurse(b0, b1, b2, b3);
        return outp;
    }

    static PathD SampleByMaxSegmentLength(PointD P0, PointD P1, PointD P2, double maxSegLen)
    {
        PathD pts = new PathD { P0 };
        SubdivideByLength(P0, P1, P2, maxSegLen, pts);
        pts.Add(P2);
        return pts;
    }

    static void SubdivideByLength(PointD p0, PointD p1, PointD p2, double maxSegLen, PathD outPts)
    {
        if (Length(Minus(p2, p0)) <= maxSegLen)
        {
            outPts.Add(p2);
            return;
        }

        PointD p01 = Mid(p0, p1), p12 = Mid(p1, p2), p012 = Mid(p01, p12);
        SubdivideByLength(p0, p01, p012, maxSegLen, outPts);
        SubdivideByLength(p012, p12, p2, maxSegLen, outPts);
    }

    static PathD SampleByMaxAngle(PointD P0, PointD P1, PointD P2, double maxAngle)
    {
        PathD pts = new PathD { P0 };
        SubdivideByAngle(P0, P1, P2, maxAngle, pts);
        pts.Add(P2);
        return pts;
    }

    static void SubdivideByAngle(PointD p0, PointD p1, PointD p2, double maxAngle, PathD outPts)
    {
        PointD tan0 = Normalized(Minus(p1, p0));
        PointD tan1 = Normalized(Minus(p2, p1));
        double dot = Math.Max(-1.0, Math.Min(1.0, Dot(tan0, tan1)));
        double angle = Math.Acos(dot);
        if (angle <= maxAngle)
        {
            outPts.Add(p2);
            return;
        }

        PointD p01 = Mid(p0, p1), p12 = Mid(p1, p2), p012 = Mid(p01, p12);
        SubdivideByAngle(p0, p01, p012, maxAngle, outPts);
        SubdivideByAngle(p012, p12, p2, maxAngle, outPts);
    }

    static PointD EstimateOutgoingTangent(PathD poly)
    {
        if (poly == null || poly.Count < 2) return new PointD(1, 0);
        PointD a = poly[^2];
        PointD b = poly[^1];
        PointD t = Minus(b, a);
        double l = Length(t);
        return l > 0 ? Mult(t, 1.0 / l) : new PointD(1, 0);
    }

    static PointD EstimateIncomingTangent(PathD poly)
    {
        if (poly == null || poly.Count < 2) return new PointD(1, 0);
        PointD a = poly[0];
        PointD b = poly[1];
        PointD t = Minus(b, a);
        double l = Length(t);
        return l > 0 ? Mult(t, 1.0 / l) : new PointD(1, 0);
    }

    // -------------------------------------------------------------------------
    // Additional builders (Cubic/Quintic/CircularArc)
    // -------------------------------------------------------------------------
    static PathD SampleCubicBezier(PointD p0, PointD b1, PointD b2, PointD p3, int samples)
    {
        PathD outp = new PathD();
        for (int i = 0; i <= samples; i++)
        {
            double t = (double)i / samples;
            double mt = 1 - t;
            double w0 = mt * mt * mt;
            double w1 = 3 * mt * mt * t;
            double w2 = 3 * mt * t * t;
            double w3 = t * t * t;
            PointD p = new PointD(
                p0.x * w0 + b1.x * w1 + b2.x * w2 + p3.x * w3,
                p0.y * w0 + b1.y * w1 + b2.y * w2 + p3.y * w3
            );
            outp.Add(p);
        }
        return outp;
    }

    static PathD SampleQuinticHermite(PointD p0, PointD p1, PointD m0, PointD m1, PointD s0, PointD s1, int samples)
    {
        PathD outp = new PathD();
        for (int i = 0; i <= samples; i++)
        {
            double t = (double)i / samples;
            double t2 = t * t;
            double t3 = t2 * t;
            double t4 = t3 * t;
            double t5 = t4 * t;
            double h0 = 1 - 10 * t3 + 15 * t4 - 6 * t5;
            double h1 = t - 6 * t3 + 8 * t4 - 3 * t5;
            double h2 = 0.5 * (t2 - 3 * t3 + 3 * t4 - t5);
            double h3 = 10 * t3 - 15 * t4 + 6 * t5;
            double h4 = -4 * t3 + 7 * t4 - 3 * t5;
            double h5 = 0.5 * (t3 - t4);
            PointD p = new PointD(
                p0.x * h0 + m0.x * h1 + s0.x * h2 + p1.x * h3 + m1.x * h4 + s1.x * h5,
                p0.y * h0 + m0.y * h1 + s0.y * h2 + p1.y * h3 + m1.y * h4 + s1.y * h5
            );
            outp.Add(p);
        }
        return outp;
    }

    static int SampleCountForInset(double ins)
    {
        int sc = Math.Max(8, (int)Math.Ceiling(ins * 40));
        return sc;
    }

    static PathD BuildCubicFromHermite(PointD p0, PointD p1, PointD m0, PointD m1, int samples)
    {
        PointD b1 = new PointD(p0.x + m0.x / 3.0, p0.y + m0.y / 3.0);
        PointD b2 = new PointD(p1.x - m1.x / 3.0, p1.y - m1.y / 3.0);
        return SampleCubicBezier(p0, b1, b2, p1, samples);
    }

    static PathD BuildC1Cubic(PointD a, PointD b, PointD tanA, PointD tanB, double inset)
    {
        double HandleScale = 0.3;
        double maxHandle = inset * HandleScale;
        PointD dirA = Normalized(tanA);
        PointD dirB = Normalized(tanB);
        if (Length(dirA) < 1e-12) dirA = Normalized(Minus(b, a));
        if (Length(dirB) < 1e-12) dirB = Normalized(Minus(b, a));
        PointD m0 = Mul(dirA, maxHandle);
        PointD m1 = Mul(dirB, maxHandle);
        int samples = SampleCountForInset(inset);
        return BuildCubicFromHermite(a, b, m0, m1, samples);
    }

    static PathD BuildQuinticC2(PointD a, PointD b, PointD tanA, PointD tanB, double inset)
    {
        double HandleScale = 0.3;
        double maxHandle = inset * HandleScale;
        PointD dirA = Normalized(tanA);
        PointD dirB = Normalized(tanB);
        if (Length(dirA) < 1e-12) dirA = Normalized(Minus(b, a));
        if (Length(dirB) < 1e-12) dirB = Normalized(Minus(b, a));
        PointD m0 = Mul(dirA, maxHandle);
        PointD m1 = Mul(dirB, maxHandle);
        PointD s0 = new PointD(0, 0);
        PointD s1 = new PointD(0, 0);
        int samples = SampleCountForInset(inset);
        return SampleQuinticHermite(a, b, m0, m1, s0, s1, samples);
    }

    static PathD BuildSmoothBlendC1(PointD a, PointD b, PointD tanA, PointD tanB, double inset)
    {
        double HandleScale = 0.35;
        double maxHandle = inset * HandleScale;
        PointD dirA = Normalized(tanA);
        PointD dirB = Normalized(tanB);
        if (Length(dirA) < 1e-12) dirA = Normalized(Minus(b, a));
        if (Length(dirB) < 1e-12) dirB = Normalized(Minus(b, a));
        PointD m0 = Mul(dirA, maxHandle);
        PointD m1 = Mul(dirB, maxHandle);
        PointD s0 = new PointD(0, 0);
        PointD s1 = new PointD(0, 0);
        int samples = SampleCountForInset(inset);
        return SampleQuinticHermite(a, b, m0, m1, s0, s1, samples);
    }

    static PathD BuildCircularArcOrNull(PointD a, PointD b, PointD tanA, PointD tanB, double inset)
    {
        if (!TryBuildTangentArc(a, tanA, b, tanB, out PointD center, out double radius, out double startAng,
                out double sweep))
        {
            return null;
        }

        int samples = SampleCountForInset(inset);
        PathD outp = new PathD();
        for (int i = 0; i <= samples; i++)
        {
            double t = (double)i / samples;
            double ang = startAng + sweep * t;
            PointD p = new PointD(center.x + radius * Math.Cos(ang), center.y + radius * Math.Sin(ang));
            outp.Add(p);
        }
        return outp;
    }

    static bool TryBuildTangentArc(PointD p0, PointD dir0, PointD p1, PointD dir1, out PointD center, out double radius,
        out double startAng, out double sweep)
    {
        center = new PointD(0, 0);
        radius = 0;
        startAng = 0;
        sweep = 0;

        double a11 = -dir0.x, a12 = -dir0.y;
        double b1 = -(-dir0.x * p0.x - dir0.y * p0.y);
        double a21 = -dir1.x, a22 = -dir1.y;
        double b2 = -(-dir1.x * p1.x - dir1.y * p1.y);
        double det = a11 * a22 - a12 * a21;
        if (Math.Abs(det) < 1e-12) return false;

        double cx = (b1 * a22 - a12 * b2) / det;
        double cy = (a11 * b2 - b1 * a21) / det;
        center = new PointD(cx, cy);

        PointD v0 = Minus(p0, center);
        PointD v1 = Minus(p1, center);
        double r0 = Length(v0);
        double r1 = Length(v1);
        if (r0 < 1e-9 || r1 < 1e-9) return false;
        if (Math.Abs(r0 - r1) > Math.Max(r0, r1) * 1e-3) return false;

        radius = 0.5 * (r0 + r1);
        startAng = Math.Atan2(v0.y, v0.x);
        double endAng = Math.Atan2(v1.y, v1.x);
        double rawSweep = endAng - startAng;
        while (rawSweep <= -Math.PI) rawSweep += 2 * Math.PI;
        while (rawSweep > Math.PI) rawSweep -= 2 * Math.PI;

        PointD tangentAtStart = new PointD(-v0.y, v0.x);
        double ssign = Math.Sign(Dot(tangentAtStart, dir0));
        if (ssign == 0) ssign = 1;
        if (ssign < 0 && rawSweep > 0) rawSweep -= 2 * Math.PI;
        if (ssign > 0 && rawSweep < 0) rawSweep += 2 * Math.PI;

        sweep = rawSweep;
        if (Math.Abs(sweep) > Math.PI * 1.5) return false;
        return true;
    }

    // -------------------------------------------------------------------------
    // CSV & SVG Utilities
    // -------------------------------------------------------------------------
    static void WriteCsv(string path, PathD pts)
    {
        StringBuilder sb = new StringBuilder();
        sb.AppendLine("X,Y");
        foreach (PointD p in pts)
            sb.AppendLine($"{p.x.ToString(CultureInfo.InvariantCulture)},{p.y.ToString(CultureInfo.InvariantCulture)}");
        File.WriteAllText(path, sb.ToString());
    }

    static string BuildDetailedSvg(IEnumerable<PointD> keyPts, PathD curvePts)
    {
        var all = keyPts.Concat(curvePts).ToList();
        double minX = all.Min(p => p.x), maxX = all.Max(p => p.x);
        double minY = all.Min(p => p.y), maxY = all.Max(p => p.y);
        double w = Math.Max(1e-6, maxX - minX), h = Math.Max(1e-6, maxY - minY);
        double mX = w * 0.1, mY = h * 0.1;
        minX -= mX;
        maxX += mX;
        minY -= mY;
        maxY += mY;

        string viewBox = $"{minX} {(-maxY)} {maxX - minX} {maxY - minY}";
        StringBuilder sb = new StringBuilder();
        sb.AppendLine($"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"600\" height=\"600\" viewBox=\"{viewBox}\">");

        double xs = NiceStep(w), ys = NiceStep(h);
        sb.AppendLine("  <g stroke=\"#ddd\" stroke-width=\"0.5\">");
        for (double x = Math.Ceiling(minX); x <= maxX; x += xs)
            sb.AppendLine($"    <line x1=\"{x}\" y1=\"{-maxY}\" x2=\"{x}\" y2=\"{-minY}\"/>");
        for (double y = Math.Floor(minY); y <= maxY; y += ys)
            sb.AppendLine($"    <line x1=\"{minX}\" y1=\"{-y}\" x2=\"{maxX}\" y2=\"{-y}\"/>");
        sb.AppendLine("  </g>");

        sb.AppendLine("  <g stroke=\"#000\" stroke-width=\"1\">");
        sb.AppendLine($"    <line x1=\"{minX}\" y1=\"0\" x2=\"{maxX}\" y2=\"0\"/>");
        sb.AppendLine($"    <line x1=\"0\" y1=\"{-maxY}\" x2=\"0\" y2=\"{-minY}\"/>");
        sb.AppendLine("  </g>");

        var pts = keyPts.ToList();
        if (pts.Count >= 4)
        {
            sb.AppendLine(DrawLine(pts[0], pts[1], "#888", 1.5));
            sb.AppendLine(DrawLine(pts[2], pts[3], "#888", 1.5));
        }

        sb.Append("  <polyline fill=\"none\" stroke=\"blue\" stroke-width=\"2\" points=\"");
        foreach (PointD p in curvePts)
            sb.Append($"{p.x},{-p.y} ");
        sb.AppendLine("\"/>");

        sb.AppendLine("</svg>");
        return sb.ToString();
    }

    static double NiceStep(double r)
    {
        double e = Math.Floor(Math.Log10(r));
        double f = r / Math.Pow(10, e);
        double bs = f <= 2 ? 0.5 : f <= 5 ? 1 : 2;
        return bs * Math.Pow(10, e);
    }

    // -------------------------------------------------------------------------
    // Small geometric helpers
    // -------------------------------------------------------------------------
    static bool PointsEqual(PointD a, PointD b) =>
        Math.Abs(a.x - b.x) < 1e-9 && Math.Abs(a.y - b.y) < 1e-9;

    static PointD Add(PointD a, PointD b) => new PointD(a.x + b.x, a.y + b.y);
    static PointD Minus(PointD a, PointD b) => new PointD(a.x - b.x, a.y - b.y);
    static PointD Mult(PointD a, double s) => new PointD(a.x * s, a.y * s);
    static double Length(PointD p) => Math.Sqrt(p.x * p.x + p.y * p.y);
    static PointD Mul(PointD a, double s) => new PointD(a.x * s, a.y * s);

    static PointD Neg(PointD a) => new PointD(-a.x, -a.y);
    static double Dot(PointD a, PointD b) => a.x * b.x + a.y * b.y;
    static PointD Mid(PointD a, PointD b) => new PointD((a.x + b.x) / 2.0, (a.y + b.y) / 2.0);

    static PointD Normalized(PointD p)
    {
        double len = Length(p);
        return len > 0 ? new PointD(p.x / len, p.y / len) : new PointD(0, 0);
    }

    static string DrawLine(PointD a, PointD b, string c, double w) =>
        $"  <line x1=\"{a.x}\" y1=\"{-a.y}\" x2=\"{b.x}\" y2=\"{-b.y}\" stroke=\"{c}\" stroke-width=\"{w}\"/>";

    static string DrawCircle(PointD p, int r, string fill) =>
        $"  <circle cx=\"{p.x}\" cy=\"{-p.y}\" r=\"{r}\" fill=\"{fill}\"/>";

    static string DrawText(string txt, PointD p, int dx, int dy) =>
        $"  <text x=\"{p.x + dx}\" y=\"{-p.y + dy}\" font-size=\"10\" fill=\"#000\">{txt}</text>";
}

// -----------------------------------------------------------------------------
// PointD: small 2D point struct with convenience operators
// -----------------------------------------------------------------------------
struct PointD
{
    public double x, y;

    public PointD(double x_, double y_)
    {
        x = x_;
        y = y_;
    }

    public static PointD operator +(PointD a, PointD b) => new PointD(a.x + b.x, a.y + b.y);
    public static PointD operator -(PointD a, PointD b) => new PointD(a.x - b.x, a.y - b.y);
    public static PointD operator *(PointD a, double s) => new PointD(a.x * s, a.y * s);
    public static PointD operator *(double s, PointD a) => new PointD(a.x * s, a.y * s);
    public double Length() => Math.Sqrt(x * x + y * y);

    public PointD Normalized()
    {
        double len = Length();
        return len > 0 ? new PointD(x / len, y / len) : new PointD(0, 0);
    }
}