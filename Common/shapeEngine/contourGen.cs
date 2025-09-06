using Clipper2Lib;
using System.Linq;
using System.Threading.Tasks;
using utility;
using geoWrangler;

namespace shapeEngine;

public static class contourGen
{
    /*
     * Improved blending from diagonal to bezier:
     * - Replaced the linear vector interpolation used to blend directions
     *   with a spherical linear interpolation (slerp) driven by an S-curve
     *   easing function. This produces a smoother rotational transition
     *   between the corner tangent and the diagonal direction (an "S-curve"
     *   ease-in/ease-out interpolation rather than a linear mix).
     *
     * - Added the missing HermiteToBezier helper (maps Hermite endpoints + tangents
     *   to cubic Bezier control points) so BuildCubicHermiteSampled and related helpers
     *   can operate.
     *
     * - Replace the sharp straight diagonal between S and E with a soft
     *   quintic Hermite segment (C2) built from the eased derivatives. This
     *   removes the abrupt curvature change where short-edge diagonals meet
     *   the corner fillets and produces a visually much softer transition.
     *
     * Everything else is preserved so you can test with existing easingMode,
     * including the previously added QuinticC2 option.
     */

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
        CubicBezierHermite,
        QuinticHermite,
        SmoothBlend,
        CircularArc,
        QuinticC2
    }

    // Fast math approximations for performance-critical operations
    static class FastMath
    {
        /// <summary>
        /// Fast approximation of acos using polynomial approximation
        /// Error is typically less than 0.0001 radians for most inputs
        /// </summary>
        [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
        public static double FastAcos(double x)
        {
            // For very high precision requirements, fall back to Math.Acos
            // This ensures compatibility with existing code that may rely on high precision
            return Math.Acos(x);
            
            // Note: The fast approximation below can be enabled for maximum performance
            // when high precision is not critical:
            /*
            // Clamp input to valid range
            x = Math.Max(-1.0, Math.Min(1.0, x));
            
            // Use polynomial approximation for better performance
            if (x >= 0.5)
            {
                return Math.Sqrt(2.0 * (1.0 - x)) * (1.5707963267948966 - x * (0.2145993 + x * 0.088972));
            }
            else if (x >= -0.5)
            {
                return 1.5707963267948966 - x * (1.0 - x * x * (0.16666667 - x * x * 0.075));
            }
            else
            {
                return 3.141592653589793 - Math.Sqrt(2.0 * (1.0 + x)) * (1.5707963267948966 + x * (0.2145993 - x * 0.088972));
            }
            */
        }

        /// <summary>
        /// Fast sine/cosine calculation using the same angle to avoid redundant calculations
        /// </summary>
        [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
        public static (double sin, double cos) FastSinCos(double angle)
        {
            return (Math.Sin(angle), Math.Cos(angle));
        }
    }

    /// <summary>
    /// Calculate tension-adjusted edge midpoint similar to ShapeLibrary's edgeMidpoints approach
    /// </summary>
    /// <param name="corner">Current corner point</param>
    /// <param name="prevCorner">Previous corner point</param>  
    /// <param name="nextCorner">Next corner point</param>
    /// <param name="afterNextCorner">Corner after next (for edge length calculation)</param>
    /// <param name="edgeTension">Edge tension value</param>
    /// <returns>Tension-adjusted midpoint between corner and nextCorner</returns>
    static PointD CalculateTensionAdjustedMidpoint(PointD corner, PointD prevCorner, PointD nextCorner, PointD afterNextCorner, double edgeTension)
    {
        // If tension is very close to default (1.0), use simple midpoint for performance
        if (Math.Abs(edgeTension - 1.0) < 1e-9)
        {
            return Helper.Mid(corner, nextCorner);
        }

        // Calculate edge lengths exactly like ShapeLibrary
        // previousEdgeLength: from current corner to previous corner  
        double previousEdgeLength = GeoWrangler.distanceBetweenPoints(corner, prevCorner);
        // currentEdgeLength: from current corner to next corner
        double currentEdgeLength = GeoWrangler.distanceBetweenPoints(corner, nextCorner);
        // nextEdgeLength: from next corner to corner after that
        double nextEdgeLength = GeoWrangler.distanceBetweenPoints(nextCorner, afterNextCorner);

        // Default to simple midpoint
        double offset = currentEdgeLength * 0.5;
        bool reverseSlide = true;

        // Apply tension calculation if previous and next edge lengths are valid
        if (previousEdgeLength > 0 && nextEdgeLength > 0 && currentEdgeLength > 0)
        {
            // Calculate ratio exactly like ShapeLibrary
            double ratio = Math.Abs(nextEdgeLength / previousEdgeLength);

            if (ratio < 1)
            {
                reverseSlide = false;
                if (ratio < 1E-2)
                {
                    ratio = 1E-2; // clamp
                }
                ratio = 1 / ratio; // normalize into expected range
            }

            // Apply sigmoid function with tension control (from ShapeLibrary)
            const double center = 1.0;
            double sigmoidResult = 1 / (1 + Math.Exp(-edgeTension * (center - ratio)));
            offset = currentEdgeLength * sigmoidResult;
        }

        // Calculate the direction vector from corner to nextCorner
        PointD direction = Helper.Minus(nextCorner, corner);
        double dirLength = Helper.Length(direction);
        
        if (dirLength < 1e-12)
        {
            return Helper.Mid(corner, nextCorner); // fallback for degenerate case
        }

        // Normalize direction and apply offset
        PointD normalizedDir = Helper.Mul(direction, 1.0 / dirLength);
        double finalOffset = reverseSlide ? (currentEdgeLength - offset) : offset;
        
        return Helper.Add(corner, Helper.Mul(normalizedDir, finalOffset));
    }

    public static PathD makeContour(PathD original_path, double concaveRadius, double convexRadius, double edgeResolution, double angularResolution, double shortEdgeLength, double maxShortEdgeLength, int optimizeCorners, double edgeTension)
    {
        return makeContour(original_path, concaveRadius, convexRadius, edgeResolution, angularResolution, shortEdgeLength, maxShortEdgeLength, optimizeCorners, enableParallel: true, edgeTension);
    }

    /// <summary>
    /// Generate contour with optional parallel processing for corner computation
    /// </summary>
    public static PathD makeContour(PathD original_path, double concaveRadius, double convexRadius, double edgeResolution, double angularResolution, double shortEdgeLength, double maxShortEdgeLength, int optimizeCorners, double edgeTension, bool enableParallel)
    {
        return makeContour(original_path, concaveRadius, convexRadius, edgeResolution, angularResolution, shortEdgeLength, maxShortEdgeLength, optimizeCorners, enableParallel, edgeTension);
    }

    /// <summary>
    /// Generate contour with optional parallel processing and edge tension control
    /// </summary>
    /// <param name="original_path">Original polygon path</param>
    /// <param name="concaveRadius">Radius for concave corners</param>
    /// <param name="convexRadius">Radius for convex corners</param>
    /// <param name="edgeResolution">Edge resolution for sampling</param>
    /// <param name="angularResolution">Angular resolution for sampling</param>
    /// <param name="shortEdgeLength">Threshold for short edges</param>
    /// <param name="maxShortEdgeLength">Maximum short edge length</param>
    /// <param name="optimizeCorners">Corner optimization level</param>
    /// <param name="enableParallel">Enable parallel processing</param>
    /// <param name="edgeTension">Edge tension control (similar to ShapeLibrary eTension)</param>
    /// <returns>Processed contour path</returns>
    public static PathD makeContour(PathD original_path, double concaveRadius, double convexRadius, double edgeResolution, double angularResolution, double shortEdgeLength, double maxShortEdgeLength, int optimizeCorners, bool enableParallel, double edgeTension)
    {
        SmoothStepMode smoothMode = SmoothStepMode.Quintic;

        int[] corner_types = CategorizeCorners(original_path, shortEdgeLength);
        
        PathsD processed = new PathsD();
        List<PointD> cornerMidpoints = new List<PointD>();
        List<PointD> cornerVertices = new List<PointD>();

        // Pre-allocate collections for thread safety
        int cornerCount = original_path.Count - 1;
        processed = new PathsD(new PathD[cornerCount]);
        cornerMidpoints = new List<PointD>(new PointD[cornerCount]);
        cornerVertices = new List<PointD>(new PointD[cornerCount]);

        // Process corners - parallel processing for performance on complex shapes
        if (enableParallel && cornerCount > 4) // Use parallel for more than 4 corners
        {
            var parallelOptions = new ParallelOptions
            {
                MaxDegreeOfParallelism = Math.Min(Environment.ProcessorCount, cornerCount)
            };

            ParallelProcessing.OptimizedParallelFor(0, cornerCount, parallelOptions, i =>
            {
                cornerVertices[i] = original_path[i];

                // Calculate corner indices for tension calculation
                int prevIdx = (i == 0) ? cornerCount - 1 : i - 1;
                int nextIdx = (i == cornerCount - 1) ? 0 : i + 1;
                int afterNextIdx = (nextIdx == cornerCount - 1) ? 0 : nextIdx + 1;

                PointD prevCorner = original_path[prevIdx];
                PointD currentCorner = original_path[i];
                PointD nextCorner = original_path[nextIdx];
                PointD afterNextCorner = original_path[afterNextIdx];

                // Calculate tension-adjusted midpoints
                PointD prevMid = CalculateTensionAdjustedMidpoint(currentCorner, 
                    (i == 0) ? original_path[^2] : original_path[prevIdx], 
                    prevCorner, currentCorner, edgeTension);

                PointD nextMid = CalculateTensionAdjustedMidpoint(currentCorner, prevCorner, nextCorner, afterNextCorner, edgeTension);

                var startLine = new PathD { original_path[i], prevMid };
                var endLine = new PathD { original_path[i], nextMid };

                cornerMidpoints[i] = nextMid;

                double radius = convexRadius;
                if (corner_types[i] == (int)CornerType.Concave)
                    radius = concaveRadius;

                PathD current_corner = ProcessCorner(startLine, endLine, radius, angularResolution, edgeResolution,
                    SamplingMode.ByMaxAngle);
                processed[i] = current_corner;
            });
        }
        else
        {
            // Sequential processing for simple cases or when parallel is disabled
            for (int i = 0; i < cornerCount; i++)
            {
                cornerVertices[i] = original_path[i];

                // Calculate corner indices for tension calculation
                int prevIdx = (i == 0) ? cornerCount - 1 : i - 1;
                int nextIdx = (i == cornerCount - 1) ? 0 : i + 1;
                int afterNextIdx = (nextIdx == cornerCount - 1) ? 0 : nextIdx + 1;

                PointD prevCorner = original_path[prevIdx];
                PointD currentCorner = original_path[i];
                PointD nextCorner = original_path[nextIdx];
                PointD afterNextCorner = original_path[afterNextIdx];

                // Calculate tension-adjusted midpoints
                PointD prevMid = CalculateTensionAdjustedMidpoint(currentCorner, 
                    (i == 0) ? original_path[^2] : original_path[prevIdx], 
                    prevCorner, currentCorner, edgeTension);

                PointD nextMid = CalculateTensionAdjustedMidpoint(currentCorner, prevCorner, nextCorner, afterNextCorner, edgeTension);

                var startLine = new PathD { original_path[i], prevMid };
                var endLine = new PathD { original_path[i], nextMid };

                cornerMidpoints[i] = nextMid;

                double radius = convexRadius;
                if (corner_types[i] == (int)CornerType.Concave)
                    radius = concaveRadius;

                PathD current_corner = ProcessCorner(startLine, endLine, radius, angularResolution, edgeResolution,
                    SamplingMode.ByMaxAngle);
                processed[i] = current_corner;
            }
        }
        
        EasingStrategy easingMode = EasingStrategy.QuinticC2; // try the new C2 option
        double insetFraction = 0.12;
        double minInset = 2.0;
        double maxInset = 40.0;
        int diagStraightSample = 10;

        PathD assembled = AssembleWithEasing(processed, corner_types, easingMode, insetFraction, minInset, maxInset,
            diagStraightSample, cornerMidpoints, cornerVertices, shortEdgeLength, maxShortEdgeLength, smoothMode,
            edgeResolution);

        PathD smoothed = SmoothPolylineSeams(assembled, angleThresholdDeg: 8.0, neighborhoodRadius: 3, sampleLenHint: edgeResolution);

        if (optimizeCorners == 1)
        {
            return DecimateToTargetResolution(smoothed, edgeResolution);
        }
        else
        {
            return smoothed;
            
        }
    }

    /// <summary>
    /// Decimates a high-resolution point list to meet a target segment length resolution.
    /// Preserves curve fidelity while reducing point count to match the user-specified resolution.
    /// </summary>
    /// <param name="finePoints">High-resolution curve points</param>
    /// <param name="targetResolution">Target maximum segment length</param>
    /// <returns>Decimated point list meeting target resolution</returns>
    static PathD DecimateToTargetResolution(PathD finePoints, double targetResolution)
    {
        if (finePoints.Count <= 2)
            return new PathD(finePoints);

        // Pre-allocate with estimated capacity to reduce reallocations
        int estimatedCapacity = Math.Max(4, finePoints.Count / 3);
        PathD decimated = new PathD(estimatedCapacity);
        decimated.Add(finePoints[0]); // Always keep start point

        // Calculate total curve length to determine if we need special handling for small curves
        double totalLength = 0;
        for (int i = 1; i < finePoints.Count; i++)
        {
            totalLength += Helper.Length(Helper.Minus(finePoints[i], finePoints[i-1]));
        }

        // For very small curves relative to target resolution, ensure we keep some intermediate points
        if (totalLength < targetResolution * 2)
        {
            // For small curves, use a more conservative approach
            // Keep every nth point to ensure reasonable representation
            int step = Math.Max(1, finePoints.Count / 4); // Keep roughly 4 points total
            for (int i = step; i < finePoints.Count - 1; i += step)
            {
                decimated.Add(finePoints[i]);
            }
        }
        else
        {
            // Normal decimation for larger curves
            PointD lastKept = finePoints[0];

            for (int i = 1; i < finePoints.Count - 1; i++)
            {
                PointD current = finePoints[i];
                double distanceFromLast = Helper.Length(Helper.Minus(current, lastKept));

                // Keep point if it's far enough from the last kept point
                if (distanceFromLast >= targetResolution)
                {
                    decimated.Add(current);
                    lastKept = current;
                }
            }
        }

        // Always keep end point
        decimated.Add(finePoints[^1]);

        return decimated;
    }
    static int[] CategorizeCorners(PathD path_, double short_edge_length)
    {
        int[] status = new int[path_.Count];

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

        for (int i = 0; i < path.Count; i++)
        {
            PointD prev = path[(i - 1 + path.Count) % path.Count];
            PointD curr = path[i];
            PointD next = path[(i + 1) % path.Count];

            double vx1 = curr.x - prev.x;
            double vy1 = curr.y - prev.y;
            double vx2 = next.x - curr.x;
            double vy2 = next.y - curr.y;

            double len1Sq = vx1 * vx1 + vy1 * vy1;
            double len2Sq = vx2 * vx2 + vy2 * vy2;
            double shortEdgeLengthSq = short_edge_length * short_edge_length;

            if (len1Sq <= shortEdgeLengthSq && len2Sq <= shortEdgeLengthSq)
            {
                status[i] = (int)CornerType.ShortEdge;
                continue;
            }

            double crossZ = vx1 * vy2 - vy1 * vx2;
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

    static PathD ProcessCorner(PathD startLine, PathD endLine, double radius, double angular_resolution,
        double edge_resolution, SamplingMode mode = SamplingMode.ByMaxAngle)
    {
        PointD startLineStart = startLine[0];
        PointD startLineEnd = startLine[1];
        PointD endLineStart = endLine[0];
        PointD endLineEnd = endLine[1];

        PointD startLength = Helper.Minus(startLineEnd, startLineStart);
        PointD startDir = Helper.Normalized(startLength);
        PointD endLength = Helper.Minus(endLineEnd, endLineStart);
        PointD endDir = Helper.Normalized(endLength);

        double start_radius = radius;
        double startLengthSq = startLength.x * startLength.x + startLength.y * startLength.y;
        double half_edge_length = Math.Sqrt(startLengthSq);
        if (start_radius > half_edge_length) start_radius = half_edge_length;
        PointD curveStartPoint = Helper.Add(startLineStart, Helper.Mult(startDir, start_radius));

        double end_radius = radius;
        double endLengthSq = endLength.x * endLength.x + endLength.y * endLength.y;
        half_edge_length = Math.Sqrt(endLengthSq);
        if (end_radius > half_edge_length) end_radius = half_edge_length;
        PointD curveEndPoint = Helper.Add(endLineStart, Helper.Mult(endDir, end_radius));

        double dx = curveEndPoint.x - curveStartPoint.x;
        double dy = curveEndPoint.y - curveStartPoint.y;
        double det = startDir.x * endDir.y - startDir.y * endDir.x;
        if (Math.Abs(det) < 1e-9)
        {
            PathD straight = new PathD();
            straight.Add(curveStartPoint);
            straight.Add(curveEndPoint);
            return straight;
        }

        double tParam = (dx * endDir.y - dy * endDir.x) / det;
        PointD controlPoint = Helper.Add(curveStartPoint, Helper.Mult(startDir, tParam));

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

    /// <summary>
    /// Connects midpoints when all edges are classified as short.
    /// This solves the problem where the algorithm would loop infinitely 
    /// trying to find non-short corners when none exist.
    /// For example, a square with all short edges becomes a diamond.
    /// </summary>
    /// <param name="cornerMidpoints">List of corner midpoints</param>
    /// <returns>Path connecting all midpoints</returns>
    static PathD ConnectCornerMidpoints(List<PointD> cornerMidpoints)
    {
        return new PathD(cornerMidpoints);
    }

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

        static double QuinticSmoothstep(double t) 
        {
            if (t <= 0) return 0;
            if (t >= 1) return 1;
            double t2 = t * t;
            double t3 = t2 * t;
            double t4 = t3 * t;
            double t5 = t4 * t;
            return 6.0 * t5 - 15.0 * t4 + 10.0 * t3;
        }

        static double ApplySmooth(SmoothStepMode mode, double t) =>
            mode switch
            {
                SmoothStepMode.Cubic => CubicSmoothstep(t),
                SmoothStepMode.Quintic => QuinticSmoothstep(t),
                _ => Math.Clamp(t, 0.0, 1.0)
            };

        // Spherical linear interpolation between two unit direction vectors.
        static PointD SlerpDir(PointD a, PointD b, double t)
        {
            // a, b should be unit length; handle small angles gracefully
            double dot = Math.Max(-1.0, Math.Min(1.0, Helper.Dot(a, b)));
            
            // Fast exit for nearly identical vectors
            if (dot > 0.9999) return a;
            
            double theta = FastMath.FastAcos(dot);
            if (Math.Abs(theta) < 1e-6) return a; // nearly identical: no rotation needed
            double sinTheta = Math.Sin(theta);
            if (Math.Abs(sinTheta) < 1e-12) return a;
            double wA = Math.Sin((1 - t) * theta) / sinTheta;
            double wB = Math.Sin(t * theta) / sinTheta;
            return Helper.Normalized(Helper.Add(Helper.Mul(a, wA), Helper.Mul(b, wB)));
        }

        bool[] isShort = new bool[n];
        for (int i = 0; i < n; i++)
            isShort[i] = (i < corner_types.Length && corner_types[i] == (int)CornerType.ShortEdge);

        // Special case: all edges are short - connect midpoints to avoid infinite loops
        bool allShort = isShort.All(x => x);
        if (allShort)
        {
            return ConnectCornerMidpoints(cornerMidpoints);
        }

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
                        PointD last = outPts[^1];
                        PointD first = poly[0];
                        if (Math.Abs(last.x - first.x) < 1e-9 && Math.Abs(last.y - first.y) < 1e-9)
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

            double dPrev = Helper.Length(Helper.Minus(diagStartMid, vertexPrev));
            double dNext = Helper.Length(Helper.Minus(diagEndMid, vertexNext));

            double denom = (max_short_edge_length - short_edge_length);
            double tPrev = denom == 0 ? 0.0 : Math.Clamp((dPrev - short_edge_length) / denom, 0.0, 1.0);
            double tNext = denom == 0 ? 0.0 : Math.Clamp((dNext - short_edge_length) / denom, 0.0, 1.0);

            double blendPrev = ApplySmooth(smoothMode, tPrev);
            double blendNext = ApplySmooth(smoothMode, tNext);

            if (blendPrev <= 1e-9 && blendNext <= 1e-9)
            {
                if (outPts.Count == 0) outPts.Add(processedStartPt);
                else if (!Helper.PointsEqual(outPts[^1], processedStartPt)) outPts.Add(processedStartPt);

                if (diagStraightSample <= 0)
                {
                    outPts.Add(processedEndPt);
                }
                else
                {
                    PointD diagFull = Helper.Minus(processedEndPt, processedStartPt);
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

            PointD diag = Helper.Minus(processedEndPt, processedStartPt);
            double diagLen = Helper.Length(diag);
            if (diagLen <= 1e-9)
            {
                if (outPts.Count == 0 || !Helper.PointsEqual(outPts[^1], processedEndPt)) outPts.Add(processedEndPt);
                continue;
            }

            PointD diagDir = Helper.Normalized(diag);
            double inset = Math.Max(minInset, Math.Min(maxInset, diagLen * insetFraction));
            if (inset * 2.0 > diagLen) inset = diagLen * 0.5 * 0.999;

            PointD S = Helper.Add(processedStartPt, Helper.Mul(diagDir, inset));
            PointD E = Helper.Minus(processedEndPt, Helper.Mul(diagDir, inset));

            PointD prevTangent = EstimateOutgoingTangent(processedCorners[prevIdx]);
            PointD nextTangent = EstimateIncomingTangent(processedCorners[nextIdx]);

            if (Helper.Length(prevTangent) < 1e-12) prevTangent = diagDir;
            if (Helper.Length(nextTangent) < 1e-12) nextTangent = diagDir;

            PointD prevDirUnit = Helper.Normalized(prevTangent);
            PointD nextDirUnit = Helper.Normalized(nextTangent);

            // previous code applied ApplySmooth twice — keep that behavior but then use slerp
            double tPrevBlend = ApplySmooth(smoothMode, blendPrev);
            double tNextBlend = ApplySmooth(smoothMode, blendNext);

            // Use slerp with S-curve easing instead of linear mixing for direction blending.
            // blendedPrevDir rotates from prevDirUnit -> diagDir using S-curve tPrevBlend.
            PointD blendedPrevDir = SlerpDir(prevDirUnit, diagDir, tPrevBlend);
            PointD blendedNextDir = SlerpDir(nextDirUnit, diagDir, tNextBlend);

            // Hermite magnitudes
            double Lstart = Math.Max(1e-9, Helper.Length(Helper.Minus(S, processedStartPt)));
            double Lend = Math.Max(1e-9, Helper.Length(Helper.Minus(processedEndPt, E)));

            PointD H_blended_start = Helper.Mul(blendedPrevDir, Lstart);
            PointD H_blended_end = Helper.Mul(blendedNextDir, Lend);

            double magStart = Helper.Length(H_blended_start);
            double magEnd = Helper.Length(H_blended_end);

            // Ensure tangent-matching to the diagonal uses an S-curve rotation (not linear).
            // Rotate blendedPrevDir -> diagDir by eased blend again to get final directions at S/E.
            PointD finalDirAtS = SlerpDir(blendedPrevDir, diagDir, tPrevBlend);
            PointD finalDirAtE = SlerpDir(blendedNextDir, Helper.Neg(diagDir), tNextBlend);

            PointD finalDerivAtS = Helper.Mul(finalDirAtS, magStart);
            PointD finalDerivAtE = Helper.Mul(finalDirAtE, magEnd);

            PointD T0_for_start = H_blended_start;
            PointD T1_for_start = finalDerivAtS;

            PointD T0_for_end = finalDerivAtE;
            PointD T1_for_end = H_blended_end;

            if (strategy == EasingStrategy.None)
            {
                if (outPts.Count == 0) outPts.Add(processedStartPt);
                else if (!Helper.PointsEqual(outPts[^1], processedStartPt)) outPts.Add(processedStartPt);

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
                    startSeg = BuildCubicHermiteSampled(processedStartPt, S, T0_for_start, T1_for_start,
                        edgeResolution);
                    endSeg = BuildCubicHermiteSampled(E, processedEndPt, T0_for_end, T1_for_end, edgeResolution);
                    break;

                case EasingStrategy.QuinticHermite:
                case EasingStrategy.SmoothBlend:
                    startSeg = BuildCubicHermiteSampled(processedStartPt, S, T0_for_start, T1_for_start,
                        edgeResolution);
                    endSeg = BuildCubicHermiteSampled(E, processedEndPt, T0_for_end, T1_for_end, edgeResolution);
                    break;

                case EasingStrategy.CircularArc:
                    startSeg = BuildCircularArcOrNull(processedStartPt, S, T0_for_start, diagDir, inset)
                        ?? BuildCubicHermiteSampled(processedStartPt, S, T0_for_start, T1_for_start, edgeResolution);

                    endSeg = BuildCircularArcOrNull(E, processedEndPt, Helper.Neg(diagDir), T1_for_end, inset)
                        ?? BuildCubicHermiteSampled(E, processedEndPt, T0_for_end, T1_for_end, edgeResolution);
                    break;

                case EasingStrategy.QuinticC2:
                    startSeg = BuildQuinticC2Segment(processedStartPt, S, T0_for_start, T1_for_start, inset);
                    endSeg = BuildQuinticC2Segment(E, processedEndPt, T0_for_end, T1_for_end, inset);
                    break;

                default:
                    startSeg = BuildCubicHermiteSampled(processedStartPt, S, T0_for_start, T1_for_start,
                        edgeResolution);
                    endSeg = BuildCubicHermiteSampled(E, processedEndPt, T0_for_end, T1_for_end, edgeResolution);
                    break;
            }

            if (startSeg != null)
            {
                if (outPts.Count > 0 && Helper.PointsEqual(outPts[^1], startSeg[0])) outPts.AddRange(startSeg.Skip(1));
                else outPts.AddRange(startSeg);
            }
            else
            {
                if (outPts.Count == 0 || !Helper.PointsEqual(outPts[^1], S)) outPts.Add(S);
            }

            // --- Replace the previous straight diagonal with a soft Quintic Hermite middle segment.
            // Build a C2 quintic bridging S -> E using the eased derivatives finalDerivAtS/finalDerivAtE.
            // This yields a much softer curvature transition into/out-of the diagonal.
            PathD midSeg = null;
            try
            {
                midSeg = BuildQuinticC2Segment(S, E, finalDerivAtS, finalDerivAtE, inset);
            }
            catch
            {
                midSeg = null;
            }

            if (midSeg != null && midSeg.Count >= 2)
            {
                if (!Helper.PointsEqual(outPts[^1], midSeg[0])) outPts.Add(midSeg[0]);
                outPts.AddRange(midSeg.Skip(1));
            }
            else
            {
                // fallback to the previous straight sampling (should rarely happen)
                if (diagStraightSample <= 0)
                {
                    if (!Helper.PointsEqual(outPts[^1], E)) outPts.Add(E);
                }
                else
                {
                    for (int s = 1; s <= diagStraightSample; s++)
                    {
                        double t = (double)s / (diagStraightSample + 1);
                        PointD p = Helper.Add(S, Helper.Mul(diagDir, (diagLen - 2 * inset) * t));
                        outPts.Add(p);
                    }

                    outPts.Add(E);
                }
            }

            if (endSeg != null)
            {
                if (Helper.PointsEqual(outPts[^1], endSeg[0])) outPts.AddRange(endSeg.Skip(1));
                else outPts.AddRange(endSeg);
            }
            else
            {
                if (!Helper.PointsEqual(outPts[^1], processedEndPt)) outPts.Add(processedEndPt);
            }
        }

        if (outPts.Count > 0)
        {
            if (Helper.Length(Helper.Minus(outPts[0], outPts[^1])) < 1e-9)
                outPts[^1] = outPts[0];
            else
                outPts.Add(outPts[0]);
        }

        return outPts;
    }

    // Inserted helper: convert Hermite (p0,p1,T0,T1) to cubic Bezier control points.
    // Hermite->Bezier mapping (given derivative dP/dt) is:
    // B0 = P0
    // B1 = P0 + T0 / 3
    // B2 = P1 - T1 / 3
    // B3 = P1
    static void HermiteToBezier(PointD P0, PointD P1, PointD T0, PointD T1,
        out PointD B0, out PointD B1, out PointD B2, out PointD B3)
    {
        B0 = P0;
        B1 = Helper.Add(P0, Helper.Mul(T0, 1.0 / 3.0));
        B2 = Helper.Minus(P1, Helper.Mul(T1, 1.0 / 3.0));
        B3 = P1;
    }

    static PathD BuildQuinticC2Segment(PointD a, PointD b, PointD m0, PointD m1, double inset)
    {
        double secondDerivScale = 0.5;
        PointD diff = Helper.Minus(m1, m0);
        double segLen = Math.Max(1e-9, Helper.Length(Helper.Minus(b, a)));
        PointD s0 = Helper.Mul(diff, (secondDerivScale / segLen));
        PointD s1 = Helper.Mul(Helper.Neg(diff), (secondDerivScale / segLen));

        double maxMag = inset * 0.5;
        if (Helper.Length(s0) > maxMag) s0 = Helper.Mul(Helper.Normalized(s0), maxMag);
        if (Helper.Length(s1) > maxMag) s1 = Helper.Mul(Helper.Normalized(s1), maxMag);

        int samples = SampleCountForInset(inset);
        return SampleQuinticHermite(a, b, m0, m1, s0, s1, samples);
    }

    static PathD SmoothPolylineSeams(PathD pts, double angleThresholdDeg = 8.0, int neighborhoodRadius = 3, double sampleLenHint = 0.5)
    {
        if (pts == null || pts.Count < 5) return pts;

        bool closed = Helper.PointsEqual(pts[0], pts[^1]);
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

            PointD vprev = Helper.Minus(work[i], work[im1]);
            PointD vnext = Helper.Minus(work[ip1], work[i]);
            double lprev = Helper.Length(vprev);
            double lnext = Helper.Length(vnext);
            if (lprev < 1e-9 || lnext < 1e-9) continue;

            double dot = Math.Max(-1.0, Math.Min(1.0, Helper.Dot(Helper.Normalized(vprev), Helper.Normalized(vnext))));
            double ang = FastMath.FastAcos(dot);

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

            PointD dirA = Helper.Normalized(Helper.Minus(neighborA, A));
            PointD dirB = Helper.Normalized(Helper.Minus(B, neighborB));
            if (Helper.Length(dirA) < 1e-12) dirA = Helper.Normalized(Helper.Minus(B, A));
            if (Helper.Length(dirB) < 1e-12) dirB = Helper.Normalized(Helper.Minus(B, A));

            double magA = Helper.Length(Helper.Minus(neighborA, A));
            double magB = Helper.Length(Helper.Minus(B, neighborB));
            magA = Math.Max(magA, Helper.Length(Helper.Minus(B, A)) * 0.1);
            magB = Math.Max(magB, Helper.Length(Helper.Minus(B, A)) * 0.1);

            PointD T0 = Helper.Mul(dirA, magA);
            PointD T1 = Helper.Mul(dirB, magB);

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

    static PathD BuildCubicHermiteSampled(PointD A, PointD B, PointD T0, PointD T1, double maxSegLen)
    {
        HermiteToBezier(A, B, T0, T1, out PointD b0, out PointD b1, out PointD b2, out PointD b3);
        return SampleBezierByMaxSegmentLength(b0, b1, b2, b3, maxSegLen);
    }

    static PathD SampleBezierByMaxSegmentLength(PointD b0, PointD b1, PointD b2, PointD b3, double maxLen)
    {
        PathD outp = new PathD();
        outp.Add(b0);

        void Recurse(PointD p0, PointD p1, PointD p2, PointD p3)
        {
            double chord = Helper.Length(Helper.Minus(p3, p0));
            double contLen = Helper.Length(Helper.Minus(p1, p0)) + Helper.Length(Helper.Minus(p2, p1)) + Helper.Length(Helper.Minus(p3, p2));
            if (contLen - chord <= maxLen || chord <= maxLen)
            {
                outp.Add(p3);
                return;
            }

            PointD p01 = Helper.Mul(Helper.Add(p0, p1), 0.5);
            PointD p12 = Helper.Mul(Helper.Add(p1, p2), 0.5);
            PointD p23 = Helper.Mul(Helper.Add(p2, p3), 0.5);
            PointD p012 = Helper.Mul(Helper.Add(p01, p12), 0.5);
            PointD p123 = Helper.Mul(Helper.Add(p12, p23), 0.5);
            PointD p0123 = Helper.Mul(Helper.Add(p012, p123), 0.5);
            Recurse(p0, p01, p012, p0123);
            Recurse(p0123, p123, p23, p3);
        }

        Recurse(b0, b1, b2, b3);
        return outp;
    }

    static PathD SampleByMaxSegmentLength(PointD P0, PointD P1, PointD P2, double maxSegLen)
    {
        // Use a much finer resolution for initial sampling to better represent the curve
        double fineResolution = maxSegLen / 2.0;

        // Generate high-resolution curve points
        PathD finePts = new PathD();
        finePts.Add(P0);
        SubdivideByLength(P0, P1, P2, fineResolution, finePts);
        finePts.Add(P2);

        // Decimate the fine points to meet the target resolution
        return DecimateToTargetResolution(finePts, maxSegLen);
    }

    static void SubdivideByLength(PointD p0, PointD p1, PointD p2, double maxSegLen, PathD outPts)
    {
        if (Helper.Length(Helper.Minus(p2, p0)) <= maxSegLen)
        {
            outPts.Add(p2);
            return;
        }

        PointD p01 = Helper.Mid(p0, p1), p12 = Helper.Mid(p1, p2), p012 = Helper.Mid(p01, p12);
        SubdivideByLength(p0, p01, p012, maxSegLen, outPts);
        SubdivideByLength(p012, p12, p2, maxSegLen, outPts);
    }

    static PathD SampleByMaxAngle(PointD P0, PointD P1, PointD P2, double maxAngle)
    {
        // Estimate capacity based on angle - more points for sharper curves
        double distance = Helper.Length(Helper.Minus(P2, P0));
        int estimatedCapacity = Math.Max(4, (int)(distance / maxAngle * 2));
        PathD pts = new PathD(estimatedCapacity);
        pts.Add(P0);
        SubdivideByAngle(P0, P1, P2, maxAngle, pts);
        pts.Add(P2);
        return pts;
    }

    static void SubdivideByAngle(PointD p0, PointD p1, PointD p2, double maxAngle, PathD outPts)
    {
        const int maxDepth = 3; // Increased from 1 for better quality while maintaining performance
        SubdivideByAngleRecursive(p0, p1, p2, maxAngle, outPts, 0, maxDepth);
    }

    static void SubdivideByAngleRecursive(PointD p0, PointD p1, PointD p2, double maxAngle, PathD outPts, int depth, int maxDepth)
    { 
        // Prevent stack overflow by limiting recursion depth
        if (depth >= maxDepth)
        {
            outPts.Add(p2);
            return;
        }

        // Additional fallback: if segment is very small, don't subdivide further
        double segmentLength = Helper.Length(Helper.Minus(p2, p0));
        if (segmentLength < 1e-1)
        {
            outPts.Add(p2);
            return;
        }        

        PointD tan0 = Helper.Normalized((Helper.Minus(p1, p0)));
        PointD tan1 = Helper.Normalized((Helper.Minus(p2, p1)));
        double dot = Math.Max(-1.0, Math.Min(1.0, tan0.x * tan1.x + tan0.y * tan1.y));
        double angle = FastMath.FastAcos(dot);
        if (angle <= maxAngle)
        {
            outPts.Add(p2);
            return;
        }

        PointD p01 = Helper.Mid(p0, p1), p12 = Helper.Mid(p1, p2), p012 = Helper.Mid(p01, p12);
        SubdivideByAngle(p0, p01, p012, maxAngle, outPts);
        SubdivideByAngle(p012, p12, p2, maxAngle, outPts);
    }

    static PointD EstimateOutgoingTangent(PathD poly)
    {
        if (poly == null || poly.Count < 2) return new PointD(1, 0);
        PointD a = poly[^2];
        PointD b = poly[^1];
        PointD t = Helper.Minus(b, a);
        double l = Helper.Length(t);
        return l > 0 ? Helper.Mult(t, 1.0 / l) : new PointD(1, 0);
    }

    static PointD EstimateIncomingTangent(PathD poly)
    {
        if (poly == null || poly.Count < 2) return new PointD(1, 0);
        PointD a = poly[0];
        PointD b = poly[1];
        PointD t = Helper.Minus(b, a);
        double l = Helper.Length(t);
        return l > 0 ? Helper.Mult(t, 1.0 / l) : new PointD(1, 0);
    }

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

    static double SmoothStep5(double t) 
    {
        if (t <= 0) return 0;
        if (t >= 1) return 1;
        double t2 = t * t;
        double t3 = t2 * t;
        double t4 = t3 * t;
        double t5 = t4 * t;
        return 6 * t5 - 15 * t4 + 10 * t3;
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
        PointD dirA = Helper.Normalized(tanA);
        PointD dirB = Helper.Normalized(tanB);
        if (Helper.Length(dirA) < 1e-12) dirA = Helper.Normalized(Helper.Minus(b, a));
        if (Helper.Length(dirB) < 1e-12) dirB = Helper.Normalized(Helper.Minus(b, a));
        PointD m0 = Helper.Mul(dirA, maxHandle);
        PointD m1 = Helper.Mul(dirB, maxHandle);
        int samples = SampleCountForInset(inset);
        return BuildCubicFromHermite(a, b, m0, m1, samples);
    }

    static PathD BuildQuinticC2(PointD a, PointD b, PointD tanA, PointD tanB, double inset)
    {
        double HandleScale = 0.3;
        double maxHandle = inset * HandleScale;
        PointD dirA = Helper.Normalized(tanA);
        PointD dirB = Helper.Normalized(tanB);
        if (Helper.Length(dirA) < 1e-12) dirA = Helper.Normalized(Helper.Minus(b, a));
        if (Helper.Length(dirB) < 1e-12) dirB = Helper.Normalized(Helper.Minus(b, a));
        PointD m0 = Helper.Mul(dirA, maxHandle);
        PointD m1 = Helper.Mul(dirB, maxHandle);
        PointD s0 = new PointD(0, 0);
        PointD s1 = new PointD(0, 0);
        int samples = SampleCountForInset(inset);
        return SampleQuinticHermite(a, b, m0, m1, s0, s1, samples);
    }

    static PathD BuildSmoothBlendC1(PointD a, PointD b, PointD tanA, PointD tanB, double inset)
    {
        double HandleScale = 0.35;
        double maxHandle = inset * HandleScale;
        PointD dirA = Helper.Normalized(tanA);
        PointD dirB = Helper.Normalized(tanB);
        if (Helper.Length(dirA) < 1e-12) dirA = Helper.Normalized(Helper.Minus(b, a));
        if (Helper.Length(dirB) < 1e-12) dirB = Helper.Normalized(Helper.Minus(b, a));
        PointD m0 = Helper.Mul(dirA, maxHandle);
        PointD m1 = Helper.Mul(dirB, maxHandle);
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
            var (sinAng, cosAng) = FastMath.FastSinCos(ang);
            PointD p = new PointD(center.x + radius * cosAng, center.y + radius * sinAng);
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
        double r0 = Helper.Length(v0);
        double r1 = Helper.Length(v1);
        if (r0 < 1e-9 || r1 < 1e-9) return false;
        if (Math.Abs(r0 - r1) > Math.Max(r0, r1) * 1e-3) return false;
        radius = 0.5 * (r0 + r1);
        startAng = Math.Atan2(v0.y, v0.x);
        double endAng = Math.Atan2(v1.y, v1.x);
        double rawSweep = endAng - startAng;
        while (rawSweep <= -Math.PI) rawSweep += 2 * Math.PI;
        while (rawSweep > Math.PI) rawSweep -= 2 * Math.PI;
        PointD tangentAtStart = new PointD(-v0.y, v0.x);
        double ssign = Math.Sign(Helper.Dot(tangentAtStart, dir0));
        if (ssign == 0) ssign = 1;
        if (ssign < 0 && rawSweep > 0) rawSweep -= 2 * Math.PI;
        if (ssign > 0 && rawSweep < 0) rawSweep += 2 * Math.PI;
        sweep = rawSweep;
        if (Math.Abs(sweep) > Math.PI * 1.5) return false;
        return true;
    }
}