using System.Globalization;
using System.Text;

using PathD = System.Collections.Generic.List<PointD>;
using PathsD = System.Collections.Generic.List<System.Collections.Generic.List<PointD>>;

static class QuadraticBezierSamplingSwitcher_Polygon
{
    enum SamplingMode
    {
        ByMaxSegmentLength,
        ByMaxAngle
    }

private enum types { concave, convex, shortedge }
private enum SmoothStepMode { None, Cubic, Quintic }
static void Main()
{
    // Input polygon (closed: first == last). Clockwise or CCW works.
    List<PointD> original_path = new List<PointD>();
    original_path.Add(new PointD(0, 0));
    original_path.Add(new PointD(0, 120));
    original_path.Add(new PointD(20, 120));
    original_path.Add(new PointD(20, 100));
    original_path.Add(new PointD(40, 100));
    original_path.Add(new PointD(40, 80));
    original_path.Add(new PointD(59, 80));
    original_path.Add(new PointD(59, 59));
    original_path.Add(new PointD(110, 59));
    original_path.Add(new PointD(110, 0));
    original_path.Add(new PointD(0, 0));

    double concave_radius = 20.0;
    double convex_radius = 50.0;
    double edge_resolution = 0.5; // smaller -> more points on curves
    double angular_resolution = 1.0;
    double short_edge_length = 20.0;
    double max_short_edge_length = 30.0;

    // Choose smoothing mode here: None (linear), Cubic, or Quintic
    SmoothStepMode smoothMode = SmoothStepMode.Quintic;

    int[] corner_types = categorizeCorners(original_path, short_edge_length);

    // Build per-corner sampled polylines AND collect midpoints used for each corner
    PathsD processed = new PathsD();
    List<PointD> cornerMidpoints = new List<PointD>(); // midpoint on the outgoing edge for each corner
    List<PointD> cornerVertices = new List<PointD>();  // store vertex positions for each corner (original_path[i])
    for (int i = 0; i < original_path.Count - 1; i++)
    {
        cornerVertices.Add(original_path[i]);

        PathD startLine = new PathD();
        startLine.Add(original_path[i]);

        PointD prevMid;
        if (i == 0)
            prevMid = new PointD((original_path[i].x + original_path[^2].x) * 0.5, (original_path[i].y + original_path[^2].y) * 0.5);
        else
            prevMid = new PointD((original_path[i].x + original_path[i - 1].x) * 0.5, (original_path[i].y + original_path[i - 1].y) * 0.5);
        startLine.Add(prevMid);

        PathD endLine = new PathD();
        endLine.Add(original_path[i]);

        PointD nextMid;
        if (i == original_path.Count - 2)
            nextMid = new PointD((original_path[i].x + original_path[0].x) * 0.5, (original_path[i].y + original_path[0].y) * 0.5);
        else
            nextMid = new PointD((original_path[i].x + original_path[i + 1].x) * 0.5, (original_path[i].y + original_path[i + 1].y) * 0.5);
        endLine.Add(nextMid);

        // store the midpoint on the outgoing edge for index i (we'll use these for diagonal length)
        cornerMidpoints.Add(nextMid);

        double radius = convex_radius;
        if (corner_types[i] == (int)types.concave)
            radius = concave_radius;

        PathD current_corner = processCorner(startLine, endLine, radius, angular_resolution, edge_resolution, SamplingMode.ByMaxSegmentLength);
        processed.Add(current_corner);
    }

    // Easing configuration:
    EasingStrategy easingMode = EasingStrategy.CubicBezierHermite;
    double insetFraction = 0.12; // fraction of diagonal length to inset at each end (0..0.5)
    double minInset = 2.0; // minimum inset in world units
    double maxInset = 40.0; // maximum inset
    int diagStraightSample = 5; // samples for straight center section (0 = just endpoints)

    // Assemble, replacing runs of short-edge corners with diagonals
    PathD assembled = AssembleWithDiagonals(processed, corner_types, cornerMidpoints, cornerVertices,
        short_edge_length, max_short_edge_length, smoothMode, easingMode, insetFraction, minInset, maxInset,
        diagStraightSample);

    // Write SVG for inspection
    string svg = BuildDetailedSvg(original_path, assembled);
    File.WriteAllText("assembled.svg", svg, Encoding.UTF8);
    Console.WriteLine("Wrote assembled.svg (smoothMode = " + smoothMode + ")");
}

static int[] categorizeCorners(PathD path_, double short_edge_length)
{
    int[] status = new int[path_.Count];

    PathD path = new PathD(path_);
    bool trimmed_path = false;
    if (path.Count > 1 && path[0].x == path[^1].x && path[0].y == path[^1].y)
    {
        trimmed_path = true;
        path.RemoveAt(path.Count - 1);
    }

    // signed area (positive = CCW)
    double area2 = 0;
    for (int i = 0; i < path.Count; i++)
    {
        PointD p1 = path[i];
        PointD p2 = path[(i + 1) % path.Count];
        area2 += p1.x * p2.y - p2.x * p1.y;
    }
    bool isCCW = area2 > 0;

    for (int i = 0; i < path.Count; i++)
    {
        PointD prev = path[(i - 1 + path.Count) % path.Count];
        PointD curr = path[i];
        PointD next = path[(i + 1) % path.Count];

        double vx1 = curr.x - prev.x;
        double vy1 = curr.y - prev.y;
        double vx2 = next.x - curr.x;
        double vy2 = next.y - curr.y;

        double len1 = Math.Sqrt(vx1 * vx1 + vy1 * vy1);
        double len2 = Math.Sqrt(vx2 * vx2 + vy2 * vy2);

        if (len1 <= short_edge_length && len2 <= short_edge_length)
        {
            status[i] = (int)types.shortedge;
            continue;
        }

        double crossZ = vx1 * vy2 - vy1 * vx2;
        bool isVertexConvex = isCCW ? (crossZ > 0) : (crossZ < 0);
        status[i] = (int)(isVertexConvex ? types.convex : types.concave);
    }

    if (trimmed_path)
    {
        Array.Resize(ref status, status.Length + 1);
        status[^1] = status[0];
    }

    return status;
}

static PathD processCorner(PathD startLine, PathD endLine, double radius, double angular_resolution, double edge_resolution, SamplingMode mode = SamplingMode.ByMaxAngle)
{
    PointD startLineStart = startLine[0];
    PointD startLineEnd = startLine[1];
    PointD endLineStart = endLine[0];
    PointD endLineEnd = endLine[1];

    PointD startLength = Minus(startLineEnd, startLineStart);
    PointD startDir = Normalized(startLength);
    PointD endLength = Minus(endLineEnd, endLineStart);
    PointD endDir = Normalized(endLength);

    double start_radius = radius;
    double half_edge_length = Math.Abs(Math.Sqrt(startLength.x * startLength.x + startLength.y * startLength.y));
    if (start_radius > half_edge_length) start_radius = half_edge_length;
    PointD curveStartPoint = Add(startLineStart, Mult(startDir, start_radius));

    double end_radius = radius;
    half_edge_length = Math.Abs(Math.Sqrt(endLength.x * endLength.x + endLength.y * endLength.y));
    if (end_radius > half_edge_length) end_radius = half_edge_length;
    PointD curveEndPoint = Add(endLineStart, Mult(endDir, end_radius));

    double dx = curveEndPoint.x - curveStartPoint.x;
    double dy = curveEndPoint.y - curveStartPoint.y;
    double det = startDir.x * endDir.y - startDir.y * endDir.x;
    if (Math.Abs(det) < 1e-9)
    {
        // fallback: straight line between curve endpoints
        PathD straight = new PathD();
        straight.Add(curveStartPoint);
        straight.Add(curveEndPoint);
        return straight;
    }

    double tParam = (dx * endDir.y - dy * endDir.x) / det;
    PointD controlPoint = Add(curveStartPoint, Mult(startDir, tParam));

    PathD samples;
    switch (mode)
    {
        case SamplingMode.ByMaxSegmentLength:
            samples = SampleByMaxSegmentLength(curveStartPoint, controlPoint, curveEndPoint, edge_resolution);
            break;
        case SamplingMode.ByMaxAngle:
            double maxAngleRad = angular_resolution * Math.PI / 180.0;
            samples = SampleByMaxAngle(curveStartPoint, controlPoint, curveEndPoint, maxAngleRad);
            break;
        default:
            throw new ArgumentOutOfRangeException();
    }

    return samples;
}

enum EasingStrategy
{
    None,
    CubicBezierHermite, // cubic bezier constructed from Hermite derivatives (no overshoot)
    QuinticHermite, // quintic Hermite (C2)
    SmoothBlend, // quintic Hermite sampled (soft)
    CircularArc // tangent arc with cubic fallback
}

static PathD AssembleWithDiagonals(
    PathsD processedCorners,
    int[] corner_types,
    List<PointD> cornerMidpoints,
    List<PointD> cornerVertices,
    double short_edge_length,
    double max_short_edge_length,
    SmoothStepMode smoothMode,
    EasingStrategy easingStrategy,
    double insetFraction,
    double minInset,
    double maxInset,
    int diagStraightSample)
{
    int n = processedCorners.Count;
    if (cornerMidpoints == null || cornerMidpoints.Count != n)
        throw new ArgumentException("cornerMidpoints must be same length as processedCorners");
    if (cornerVertices == null || cornerVertices.Count != n)
        throw new ArgumentException("cornerVertices must be same length as processedCorners");
    PathD outPts = new PathD();

    // helper smoothstep functions (kept from original)
    static double CubicSmoothstep(double t) => t <= 0 ? 0 : t >= 1 ? 1 : (3.0 * t * t - 2.0 * t * t * t);

    static double QuinticSmoothstep(double t) =>
        t <= 0 ? 0 : t >= 1 ? 1 : (6.0 * Math.Pow(t, 5) - 15.0 * Math.Pow(t, 4) + 10.0 * Math.Pow(t, 3));

    static double ApplySmooth(SmoothStepMode mode, double t)
    {
        return mode switch
        {
            SmoothStepMode.Cubic => CubicSmoothstep(t),
            SmoothStepMode.Quintic => QuinticSmoothstep(t),
            _ => Math.Clamp(t, 0.0, 1.0)
        };
    }

    // Mark short corners (only consider first n entries of corner_types)
    bool[] isShort = new bool[n];
    for (int k = 0; k < n; k++)
        isShort[k] = (k < corner_types.Length && corner_types[k] == (int)types.shortedge);

    int i = 0;
    while (i < n)
    {
        if (!isShort[i])
        {
            // Append corner polyline, avoiding duplicating a point
            PathD poly = processedCorners[i];
            if (poly.Count > 0)
            {
                if (outPts.Count == 0)
                    outPts.AddRange(poly);
                else
                {
                    PointD last = outPts[^1];
                    PointD first = poly[0];
                    if (Math.Abs(last.x - first.x) < 1e-9 && Math.Abs(last.y - first.y) < 1e-9)
                        outPts.AddRange(poly.Skip(1));
                    else
                        outPts.AddRange(poly);
                }
            }

            i++;
            continue;
        }

        // Start of short run
        int runStart = i;
        while (i < n && isShort[i]) i++;
        int runEnd = i - 1;

        // Find previous non-short index (wrap)
        int prevIdx = (runStart - 1 + n) % n;
        while (isShort[prevIdx]) prevIdx = (prevIdx - 1 + n) % n;

        // Find next non-short index (wrap)
        int nextIdx = (runEnd + 1) % n;
        while (isShort[nextIdx]) nextIdx = (nextIdx + 1) % n;

        // Use processed endpoints for actual geometry connection
        PointD processedStartPt = processedCorners[prevIdx].Last();
        PointD processedEndPt = processedCorners[nextIdx].First();

        // Diagonal endpoints for length measurement from midpoints (these are the midpoints on edges)
        PointD diagStartMid = cornerMidpoints[prevIdx]; // midpoint on edge leaving prevIdx
        PointD diagEndMid = cornerMidpoints[nextIdx]; // midpoint on edge leaving nextIdx

        // Corner vertex positions adjacent to the short run:
        PointD vertexPrev = cornerVertices[prevIdx];
        PointD vertexNext = cornerVertices[nextIdx];

        // Distances from vertex to its adjacent midpoint (per-side lengths)
        double dPrev = Length(Minus(diagStartMid, vertexPrev)); // distance from prev vertex to its outgoing midpoint
        double dNext = Length(Minus(diagEndMid, vertexNext)); // distance from next vertex to its outgoing midpoint

        // Normalized t values in [0,1]
        double tPrev = Math.Clamp((dPrev - short_edge_length) / (max_short_edge_length - short_edge_length), 0.0, 1.0);
        double tNext = Math.Clamp((dNext - short_edge_length) / (max_short_edge_length - short_edge_length), 0.0, 1.0);

        // Apply chosen smoothing
        double blendPrev = ApplySmooth(smoothMode, tPrev);
        double blendNext = ApplySmooth(smoothMode, tNext);

        // Ensure processedStartPt is present as last point to connect diagonal from
        if (outPts.Count == 0)
            outPts.Add(processedStartPt);
        else
        {
            PointD last = outPts[^1];
            if (!(Math.Abs(last.x - processedStartPt.x) < 1e-9 && Math.Abs(last.y - processedStartPt.y) < 1e-9))
                outPts.Add(processedStartPt);
        }

        // If both blends are effectively zero, use straight diagonal from processedStartPt to processedEndPt
        if (blendPrev <= 1e-9 && blendNext <= 1e-9)
        {
            // preserve original optional sampling behavior
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

        // Compute blended endpoints by lerping between diagonal endpoint and bezier endpoint using each side blend.
        PointD bezierStartEndpoint = processedCorners[prevIdx].Last(); // actual end of prev processed corner
        PointD bezierEndEndpoint = processedCorners[nextIdx].First(); // actual start of next processed corner

        PointD blendedStart = new PointD(
            processedStartPt.x * (1 - blendPrev) + bezierStartEndpoint.x * blendPrev,
            processedStartPt.y * (1 - blendPrev) + bezierStartEndpoint.y * blendPrev
        );
        PointD blendedEnd = new PointD(
            processedEndPt.x * (1 - blendNext) + bezierEndEndpoint.x * blendNext,
            processedEndPt.y * (1 - blendNext) + bezierEndEndpoint.y * blendNext
        );

        // Compute control point using average blend to move control from diagonal midpoint toward bezier midpoint
        PointD diagControl = new PointD((processedStartPt.x + processedEndPt.x) / 2.0,
            (processedStartPt.y + processedEndPt.y) / 2.0);
        PointD bezierControl = new PointD((bezierStartEndpoint.x + bezierEndEndpoint.x) / 2.0,
            (bezierStartEndpoint.y + bezierEndEndpoint.y) / 2.0);
        double avgBlend = (blendPrev + blendNext) * 0.5;
        PointD blendedControl = new PointD(
            diagControl.x * (1 - avgBlend) + bezierControl.x * avgBlend,
            diagControl.y * (1 - avgBlend) + bezierControl.y * avgBlend
        );

        // --- New easing integration starts here ---
        // We'll compute an inset and create easing segments near start and end, then keep center straight region.
        // Compute diagonal vector and inset values (respecting provided min/max and insetFraction)
        PointD fullDiag = Minus(processedEndPt, processedStartPt);
        double fullDiagLen = Length(fullDiag);
        PointD fullDiagDir = fullDiag;
        if (fullDiagLen > 1e-12) fullDiagDir = Normalized(fullDiag);
        double inset = Math.Max(minInset, Math.Min(maxInset, fullDiagLen * insetFraction));
        if (inset * 2.0 > fullDiagLen) inset = fullDiagLen * 0.5 * 0.999;

        PointD S = Add(processedStartPt, Mul(fullDiagDir, inset)); // inward from start along diagonal
        PointD E = Minus(processedEndPt, Mul(fullDiagDir, inset)); // inward from end along diagonal

        // If easing disabled explicitly, fall back to previous blended cubic sampling
        if (easingStrategy == EasingStrategy.None)
        {
            // sample cubic using blendedStart/blendedControl/blendedEnd like original
            PathD blendedCurve = SampleByMaxSegmentLength(blendedStart, blendedControl, blendedEnd, 0.5);
            outPts.AddRange(blendedCurve);
            continue;
        }

        // Estimate tangents from processed polyline endpoints to match continuity
        PointD prevTangent = EstimateOutgoingTangent(processedCorners[prevIdx]);
        PointD nextTangent = EstimateIncomingTangent(processedCorners[nextIdx]);

        if (Length(prevTangent) < 1e-12) prevTangent = fullDiagDir;
        if (Length(nextTangent) < 1e-12) nextTangent = fullDiagDir;

        // Build easing segments for start and end based on chosen strategy.
        // Each Build* helper is expected to return a PathD starting at the first argument and ending at the second.
        PathD startSeg = null;
        PathD endSeg = null;

        switch (easingStrategy)
        {
            case EasingStrategy.CubicBezierHermite:
                // create cubic Hermite segments from processedStartPt -> S and E -> processedEndPt
                // control lengths scaled by inset to keep locality
                startSeg = BuildC1Cubic(processedStartPt, S, prevTangent, fullDiagDir, inset);
                endSeg = BuildC1Cubic(E, processedEndPt, Neg(fullDiagDir), nextTangent, inset);
                break;

            case EasingStrategy.QuinticHermite:
                // produce C2 quintic Hermite transitions
                startSeg = BuildQuinticC2(processedStartPt, S, prevTangent, fullDiagDir, inset);
                endSeg = BuildQuinticC2(E, processedEndPt, Neg(fullDiagDir), nextTangent, inset);
                break;

            case EasingStrategy.SmoothBlend:
                // softer quintic-like sampled blend
                startSeg = BuildSmoothBlendC1(processedStartPt, S, prevTangent, fullDiagDir, inset);
                endSeg = BuildSmoothBlendC1(E, processedEndPt, Neg(fullDiagDir), nextTangent, inset);
                break;

            case EasingStrategy.CircularArc:
                // Prefer circular tangent arc; fallback to cubic hermite when invalid
                startSeg = BuildCircularArcOrNull(processedStartPt, S, prevTangent, fullDiagDir, inset);
                if (startSeg == null) startSeg = BuildC1Cubic(processedStartPt, S, prevTangent, fullDiagDir, inset);
                endSeg = BuildCircularArcOrNull(E, processedEndPt, Neg(fullDiagDir), nextTangent, inset);
                if (endSeg == null) endSeg = BuildC1Cubic(E, processedEndPt, Neg(fullDiagDir), nextTangent, inset);
                break;

            default:
                startSeg = BuildC1Cubic(processedStartPt, S, prevTangent, fullDiagDir, inset);
                endSeg = BuildC1Cubic(E, processedEndPt, Neg(fullDiagDir), nextTangent, inset);
                break;
        }

        // Append start segment avoiding duplicate first point
        if (startSeg != null)
        {
            if (outPts.Count > 0 && PointsEqual(outPts[^1], startSeg[0])) outPts.AddRange(startSeg.Skip(1));
            else outPts.AddRange(startSeg);
        }
        else
        {
            if (!PointsEqual(outPts[^1], S)) outPts.Add(S);
        }

        // Center straight: sample from S to E (preserve diagStraightSample param)
        double centerLen = Length(Minus(E, S));
        if (centerLen <= 1e-12)
        {
            // Degenerate: skip to E
            if (!PointsEqual(outPts[^1], E)) outPts.Add(E);
        }
        else
        {
            if (diagStraightSample <= 0)
            {
                if (!PointsEqual(outPts[^1], E)) outPts.Add(E);
            }
            else
            {
                for (int s = 1; s <= diagStraightSample; s++)
                {
                    double t = (double)s / (diagStraightSample + 1);
                    PointD p = Add(S, Mul(fullDiagDir, centerLen * t));
                    outPts.Add(p);
                }

                outPts.Add(E);
            }
        }

        // Append end segment avoiding duplicate E
        if (endSeg != null)
        {
            if (PointsEqual(outPts[^1], endSeg[0])) outPts.AddRange(endSeg.Skip(1));
            else outPts.AddRange(endSeg);
        }
        else
        {
            if (!PointsEqual(outPts[^1], processedEndPt)) outPts.Add(processedEndPt);
        }
    }

    return outPts;
}

static PathD SampleByMaxSegmentLength(PointD P0, PointD P1, PointD P2, double maxSegLen)
    {
        PathD pts = new PathD();
        pts.Add(P0);
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
        PathD pts = new PathD();
        pts.Add(P0);
        SubdivideByAngle(P0, P1, P2, maxAngle, pts);
        pts.Add(P2);
        return pts;
    }

    static void SubdivideByAngle(PointD p0, PointD p1, PointD p2, double maxAngle, PathD outPts)
    {
        PointD tan0 = Normalized((Minus(p1, p0)));
        PointD tan1 = Normalized((Minus(p2, p1)));
        double dot = Math.Max(-1.0, Math.Min(1.0, tan0.x * tan1.x + tan0.y * tan1.y));
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

    // Estimate outgoing tangent at end of a corner polyline (last - secondLast normalized)
    static PointD EstimateOutgoingTangent(PathD poly)
    {
        if (poly == null || poly.Count < 2) return new PointD(1, 0);
        PointD a = poly[^2];
        PointD b = poly[^1];
        PointD t = Minus(b, a);
        double l = Length(t);
        return l > 0 ? Mult(t, 1.0 / l) : new PointD(1, 0);
    }

    // Estimate incoming tangent at start of a corner polyline (second - first normalized)
    static PointD EstimateIncomingTangent(PathD poly)
    {
        if (poly == null || poly.Count < 2) return new PointD(1, 0);
        PointD a = poly[0];
        PointD b = poly[1];
        PointD t = Minus(b, a);
        double l = Length(t);
        return l > 0 ? Mult(t, 1.0 / l) : new PointD(1, 0);
    }

    // Create quadratic samples between p0->p2 with control p1. Use max segment length as "sampleLenHint".
    static PathD SampleQuadratic(PointD p0, PointD p1, PointD p2, double sampleLenHint)
    {
        // For these small easings, interpret sampleLenHint as desired maximum chord length
        double maxSeg = Math.Max(0.5, sampleLenHint); // guard smallness
        return SampleByMaxSegmentLength(p0, p1, p2, maxSeg);
    }

    // Convert Hermite (p0,p1,m0,m1) to cubic bezier and sample it
    static PathD HermiteToCubicAndSample(PointD p0, PointD p1, PointD m0, PointD m1, double sampleLenHint)
    {
        // Cubic Bezier control points:
        // b0 = p0
        // b1 = p0 + m0/3
        // b2 = p1 - m1/3
        // b3 = p1
        PointD b0 = p0;
        PointD b1 = p0 + (1.0 / 3.0) * m0;
        PointD b2 = p1 - (1.0 / 3.0) * m1;
        PointD b3 = p1;
        // Convert cubic to a polyline by simple recursive subdivision using chord length threshold
        return SampleCubicByLength(b0, b1, b2, b3, Math.Max(0.5, sampleLenHint));
    }

    // Simple cubic subdivision sampling by chord length
    static PathD SampleCubicByLength(PointD b0, PointD b1, PointD b2, PointD b3, double maxSegLen)
    {
        PathD pts = new PathD();
        pts.Add(b0);
        SubdivideCubicByLength(b0, b1, b2, b3, maxSegLen, pts);
        pts.Add(b3);
        return pts;
    }

    static void SubdivideCubicByLength(PointD p0, PointD p1, PointD p2, PointD p3, double maxSegLen, PathD outPts)
    {
        // If chord length small enough, emit p3
        if (Length(Minus(p3, p0)) <= maxSegLen)
        {
            outPts.Add(p3);
            return;
        }

        // de Casteljau split cubic
        PointD p01 = Mid(p0, p1), p12 = Mid(p1, p2), p23 = Mid(p2, p3);
        PointD p012 = Mid(p01, p12), p123 = Mid(p12, p23);
        PointD p0123 = Mid(p012, p123);
        SubdivideCubicByLength(p0, p01, p012, p0123, maxSegLen, outPts);
        SubdivideCubicByLength(p0123, p123, p23, p3, maxSegLen, outPts);
    }

        // Solve circle center for tangency at p0 with dir0 and p1 with dir1 (used by CircularArc)
    static bool TryBuildTangentArc(PointD p0, PointD dir0, PointD p1, PointD dir1, out PointD center, out double radius,
        out double startAng, out double sweep)
    {
        center = new PointD(0, 0);
        radius = 0;
        startAng = 0;
        sweep = 0;
        // Equations: (p0 - C) dot dir0 = 0 and (p1 - C) dot dir1 = 0
        double a11 = -dir0.x;
        double a12 = -dir0.y;
        double b1 = -(-dir0.x * p0.x - dir0.y * p0.y);
        double a21 = -dir1.x;
        double a22 = -dir1.y;
        double b2 = -(-dir1.x * p1.x - dir1.y * p1.y);
        double det = a11 * a22 - a12 * a21;
        if (Math.Abs(det) < 1e-12) return false;
        double cx = (b1 * a22 - a12 * b2) / det;
        double cy = (a11 * b2 - b1 * a21) / det;
        center = new PointD(cx, cy);
        PointD v0 = new PointD(p0.x - cx, p0.y - cy);
        PointD v1 = new PointD(p1.x - cx, p1.y - cy);
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
        // Ensure orientation matches tangent directions
        PointD tangentAtStart = new PointD(-v0.y, v0.x);
        double ssign = Math.Sign(Dot(tangentAtStart, dir0));
        if (ssign == 0) ssign = 1;
        if (ssign < 0 && rawSweep > 0) rawSweep -= 2 * Math.PI;
        if (ssign > 0 && rawSweep < 0) rawSweep += 2 * Math.PI;
        sweep = rawSweep;
        if (Math.Abs(sweep) > Math.PI * 1.5) return false;
        return true;
    }

    
    // uniform cubic Bezier sampler (inclusive endpoints)
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

    // quintic hermite sampler (C2)
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

    // quintic smoothstep (C2) — used only when needed
    static double SmoothStep5(double t) =>
        t <= 0 ? 0 : t >= 1 ? 1 : (6 * Math.Pow(t, 5) - 15 * Math.Pow(t, 4) + 10 * Math.Pow(t, 3));
    // sample count heuristic
    static int SampleCountForInset(double ins)
    {
        int sc = Math.Max(8, (int)Math.Ceiling(ins * 40));
        return sc;
    }
        // ---- Builders that guarantee C1 and avoid overshoot ----
    // Convert Hermite derivatives m0,m1 to cubic Bezier controls (exact)
    static PathD BuildCubicFromHermite(PointD p0, PointD p1, PointD m0, PointD m1, int samples)
    {
        // conversion: b1 = p0 + m0/3, b2 = p1 - m1/3
        PointD b1 = new PointD(p0.x + m0.x / 3.0, p0.y + m0.y / 3.0);
        PointD b2 = new PointD(p1.x - m1.x / 3.0, p1.y - m1.y / 3.0);
        return SampleCubicBezier(p0, b1, b2, p1, samples);
    }

    // C1 cubic bezier builder that uses conservative handle magnitudes to avoid overshoot
    static PathD BuildC1Cubic(PointD a, PointD b, PointD tanA, PointD tanB, double inset)
    {
        double HandleScale = 0.3; // conservative default — reduce to avoid overshoot, increase to soften
        double maxHandle = inset * HandleScale;
        // compute derivative vectors m0 (at a) and m1 (at b). We want these to be actual derivative vectors, not unit directions.
        PointD dirA = Normalized(tanA);
        PointD dirB = Normalized(tanB);
        if (Length(dirA) < 1e-12) dirA = Normalized(Minus(b, a)); // fallback
        if (Length(dirB) < 1e-12) dirB = Normalized(Minus(b, a));
        PointD m0 = Mul(dirA, maxHandle);
        PointD m1 = Mul(dirB, maxHandle);
        int samples = SampleCountForInset(inset);
        return BuildCubicFromHermite(a, b, m0, m1, samples);
    }

    // Quintic Hermite builder (C2) with conservative derivative magnitudes
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

    // Smooth blend implemented by sampling a quintic hermite (no remap with zero endpoint slope)
    static PathD BuildSmoothBlendC1(PointD a, PointD b, PointD tanA, PointD tanB, double inset)
    {
        // Use quintic hermite but slightly bigger handle for a softer shape while keeping it conservative.
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

    // Circular arc builder; returns null if cannot build. Caller should fallback to BuildC1Cubic.
    static PathD BuildCircularArcOrNull(PointD a, PointD b, PointD tanA, PointD tanB, double inset)
    {
        // For tangency to diag direction at b we will pass tanB = -diagDir typically.
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

    // --- CSV & SVG Utilities ---
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
        minX -= mX; maxX += mX; minY -= mY; maxY += mY;

        string viewBox = $"{minX} {(-maxY)} {maxX - minX} {maxY - minY}";
        StringBuilder sb = new StringBuilder();
        sb.AppendLine($"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"600\" height=\"600\" viewBox=\"{viewBox}\">");

        // grid
        double xs = NiceStep(w), ys = NiceStep(h);
        sb.AppendLine("  <g stroke=\"#ddd\" stroke-width=\"0.5\">");
        for (double x = Math.Ceiling(minX); x <= maxX; x += xs)
            sb.AppendLine($"    <line x1=\"{x}\" y1=\"{-maxY}\" x2=\"{x}\" y2=\"{-minY}\"/>");
        for (double y = Math.Floor(minY); y <= maxY; y += ys)
            sb.AppendLine($"    <line x1=\"{minX}\" y1=\"{-y}\" x2=\"{maxX}\" y2=\"{-y}\"/>");
        sb.AppendLine("  </g>");

        // axes
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

    static bool PointsEqual(PointD a, PointD b)
    {
        return Math.Abs(a.x - b.x) < 1e-9 && Math.Abs(a.y - b.y) < 1e-9;
    }

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