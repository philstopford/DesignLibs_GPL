using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;

using PathD = System.Collections.Generic.List<PointD>;
using PathsD = System.Collections.Generic.List<System.Collections.Generic.List<PointD>>;

/// <summary>
/// Represents a 2D point with double precision coordinates for polygon processing.
/// This version has the corrected multiplication operator.
/// </summary>
public struct PointD
{
    public double x, y;
    
    public PointD(double x_, double y_) { x = x_; y = y_; }
    
    public static PointD operator +(PointD a, PointD b) => new PointD(a.x + b.x, a.y + b.y);
    public static PointD operator -(PointD a, PointD b) => new PointD(a.x - b.x, a.y - b.y);
    public static PointD operator *(PointD a, double s) => new PointD(a.x * s, a.y * s);  // Corrected
    
    public double Length() => Math.Sqrt(x * x + y * y);
    
    public PointD Normalized()
    {
        double len = Length();
        return len > 0 ? new PointD(x / len, y / len) : new PointD(0, 0);
    }
}

/// <summary>
/// Advanced polygon corner processing with short edge detection and adaptive curve fitting.
/// 
/// This prototype extends the basic corner processing algorithm by adding special handling
/// for "short edges" - cases where adjacent polygon edges are significantly shorter than
/// the specified curve radius. Instead of applying the standard curve fitting algorithm
/// (which might create curves larger than the available edge space), the algorithm detects
/// these cases and applies appropriate geometric constraints.
/// 
/// Algorithm Enhancement over Basic Version:
/// 1. **Short Edge Detection**: Identifies corners where adjacent edges are shorter than a threshold
/// 2. **Adaptive Radius Limiting**: Automatically reduces curve radius to fit available edge space
/// 3. **Special Case Handling**: Different processing strategies for short vs. normal edges
/// 4. **Geometric Validation**: Ensures generated curves don't exceed edge boundaries
/// 
/// Key Applications:
/// - PCB trace routing with fine-pitch components (small features mixed with large)
/// - Architectural drawings with mixed detail levels (doors/windows vs. walls)
/// - CNC machining toolpaths where material removal constraints vary by feature size
/// - 3D printing where layer adhesion requirements differ for small vs. large features
/// 
/// Mathematical Approach:
/// - Edge length analysis using Euclidean distance computation
/// - Threshold-based classification of edges as "short" or "normal"
/// - Proportional radius scaling to maintain geometric validity
/// - Preservation of curve continuity across mixed edge types
/// 
/// This approach ensures that the corner processing algorithm remains robust when
/// applied to polygons containing features at multiple scales.
/// </summary>
public static class QuadraticBezierSamplingSwitcher_Polygon
{
    /// <summary>
    /// Curve sampling strategies for discretization.
    /// </summary>
    public enum SamplingMode
    {
        /// <summary>Subdivide based on maximum segment length</summary>
        ByMaxSegmentLength,
        /// <summary>Subdivide based on maximum angular change</summary>
        ByMaxAngle
    }

    /// <summary>
    /// Corner classification types including short edge detection.
    /// </summary>
    private enum CornerType { Concave, Convex, ShortEdge }

    static void Main()
    {
        // Test polygon with mixed edge lengths to demonstrate short edge handling
        List<PointD> original_path = new List<PointD>();
        original_path.Add(new PointD(0, 0));
        original_path.Add(new PointD(0, 100));      // Long vertical edge
        original_path.Add(new PointD(20, 100));     // Medium horizontal edge
        original_path.Add(new PointD(20, 80));      // Short vertical edge (20 units)
        original_path.Add(new PointD(40, 80));      // Short horizontal edge (20 units)
        original_path.Add(new PointD(40, 60));      // Short vertical edge (20 units)
        original_path.Add(new PointD(60, 60));      // Short horizontal edge (20 units)
        original_path.Add(new PointD(60, 40));      // Short vertical edge (20 units)
        original_path.Add(new PointD(80, 40));      // Short horizontal edge (20 units)
        original_path.Add(new PointD(80, 0));       // Medium vertical edge
        original_path.Add(new PointD(0, 0));        // Long horizontal edge

        // Processing parameters
        double concave_radius = 20.0;        // Base radius for inward corners
        double convex_radius = 50.0;         // Base radius for outward corners
        double edge_resolution = 0.5;        // Sampling density
        double angular_resolution = 1.0;     // Angular sampling threshold  
        double short_edge_length = 20.0;     // Threshold for short edge detection

        // Classify corners with short edge awareness
        int[] corner_types = CategorizeCorners(original_path, short_edge_length);

        // Process each corner with adaptive radius limiting
        PathsD processed = ProcessAllCornersWithShortEdgeHandling(original_path, corner_types,
            concave_radius, convex_radius, edge_resolution, angular_resolution, short_edge_length);

        // Assemble and export final result
        PathD assembled = AssembleProcessedCorners(processed);
        
        string svg = BuildDetailedSvg(original_path, assembled);
        File.WriteAllText("assembled_shortedge.svg", svg, Encoding.UTF8);
        
        WriteCsv("shortedge_corners.csv", assembled);
    }

    /// <summary>
    /// Processes all corners with enhanced short edge detection and handling.
    /// </summary>
    /// <param name="original_path">Input polygon vertices</param>
    /// <param name="corner_types">Pre-classified corner types</param>
    /// <param name="concave_radius">Base radius for concave corners</param>
    /// <param name="convex_radius">Base radius for convex corners</param>
    /// <param name="edge_resolution">Sampling density parameter</param>
    /// <param name="angular_resolution">Angular sampling threshold</param>
    /// <param name="short_edge_threshold">Length threshold for short edge detection</param>
    /// <returns>List of processed corner curves</returns>
    public static PathsD ProcessAllCornersWithShortEdgeHandling(List<PointD> original_path, 
        int[] corner_types, double concave_radius, double convex_radius, 
        double edge_resolution, double angular_resolution, double short_edge_threshold)
    {
        PathsD processed = new PathsD();
        
        for (int i = 0; i < original_path.Count - 1; i++)
        {
            // Create tangent line definitions for this corner
            PathD startLine = CreateTangentLine(original_path, i, true);
            PathD endLine = CreateTangentLine(original_path, i, false);

            // Select base radius based on corner type
            double baseRadius = corner_types[i] == (int)CornerType.Concave ? concave_radius : convex_radius;
            
            // Apply short edge limitations
            double adaptiveRadius = ComputeAdaptiveRadius(startLine, endLine, baseRadius, short_edge_threshold);

            // Generate corner curve with adaptive parameters
            PathD corner_curve = ProcessCornerWithAdaptiveRadius(startLine, endLine, 
                adaptiveRadius, angular_resolution, edge_resolution);
            
            processed.Add(corner_curve);
        }
        
        return processed;
    }

    /// <summary>
    /// Processes a corner with an adaptive radius, using the enhanced algorithm.
    /// </summary>
    /// <param name="startLine">Incoming edge definition</param>
    /// <param name="endLine">Outgoing edge definition</param>
    /// <param name="adaptiveRadius">Radius adjusted for edge constraints</param>
    /// <param name="angular_resolution">Angular sampling threshold</param>
    /// <param name="edge_resolution">Edge sampling density</param>
    /// <returns>Processed corner curve points</returns>
    public static PathD ProcessCornerWithAdaptiveRadius(PathD startLine, PathD endLine, 
        double adaptiveRadius, double angular_resolution, double edge_resolution)
    {
        return processCorner(startLine, endLine, adaptiveRadius, angular_resolution, edge_resolution, SamplingMode.ByMaxSegmentLength);
    }

    /// <summary>
    /// Computes an adaptive radius that respects edge length constraints.
    /// 
    /// This is the key innovation of the short edge handling algorithm.
    /// It ensures that the curve radius never exceeds what the available
    /// edge geometry can support, preventing self-intersection and
    /// maintaining geometric validity.
    /// </summary>
    /// <param name="startLine">Incoming edge definition</param>
    /// <param name="endLine">Outgoing edge definition</param>
    /// <param name="desiredRadius">Base radius from corner type</param>
    /// <param name="shortEdgeThreshold">Length threshold for short edge detection</param>
    /// <returns>Geometrically valid radius value</returns>
    public static double ComputeAdaptiveRadius(PathD startLine, PathD endLine, 
        double desiredRadius, double shortEdgeThreshold)
    {
        // Compute actual edge lengths
        double startEdgeLength = ComputeEdgeLength(startLine);
        double endEdgeLength = ComputeEdgeLength(endLine);
        
        // Check if either edge qualifies as "short"
        bool hasShortEdge = (startEdgeLength <= shortEdgeThreshold) || (endEdgeLength <= shortEdgeThreshold);
        
        if (hasShortEdge)
        {
            // For short edges, limit radius to a fraction of the shortest edge
            double minEdgeLength = Math.Min(startEdgeLength, endEdgeLength);
            double maxAllowableRadius = minEdgeLength * 0.4;  // Use 40% of edge length
            return Math.Min(desiredRadius, maxAllowableRadius);
        }
        
        // For normal edges, use desired radius with standard half-edge limiting
        double halfStartEdge = startEdgeLength * 0.5;
        double halfEndEdge = endEdgeLength * 0.5;
        double maxStandardRadius = Math.Min(halfStartEdge, halfEndEdge);
        
        return Math.Min(desiredRadius, maxStandardRadius);
    }

    /// <summary>
    /// Computes the length of an edge defined by a two-point line.
    /// </summary>
    private static double ComputeEdgeLength(PathD line)
    {
        if (line.Count < 2) return 0;
        
        var delta = line[1] - line[0];
        return delta.Length();
    }

    /// <summary>
    /// Creates a tangent line definition for corner processing.
    /// </summary>
    private static PathD CreateTangentLine(List<PointD> path, int cornerIndex, bool isStartLine)
    {
        PathD line = new PathD();
        line.Add(path[cornerIndex]);

        PointD midPoint;
        if (isStartLine)
        {
            // Previous edge midpoint
            if (cornerIndex == 0)
                midPoint = new PointD((path[cornerIndex].x + path[^2].x) * 0.5, 
                                    (path[cornerIndex].y + path[^2].y) * 0.5);
            else
                midPoint = new PointD((path[cornerIndex].x + path[cornerIndex - 1].x) * 0.5, 
                                    (path[cornerIndex].y + path[cornerIndex - 1].y) * 0.5);
        }
        else
        {
            // Next edge midpoint
            if (cornerIndex == path.Count - 2)
                midPoint = new PointD((path[cornerIndex].x + path[0].x) * 0.5, 
                                    (path[cornerIndex].y + path[0].y) * 0.5);
            else
                midPoint = new PointD((path[cornerIndex].x + path[cornerIndex + 1].x) * 0.5, 
                                    (path[cornerIndex].y + path[cornerIndex + 1].y) * 0.5);
        }

        line.Add(midPoint);
        return line;
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

    static PathD SampleByMaxSegmentLength(PointD P0, PointD P1, PointD P2, double maxSegLen)
    {
        // Use a much finer resolution for initial sampling to better represent the curve
        double fineResolution = maxSegLen / 10.0;
        
        // Generate high-resolution curve points
        PathD finePts = new PathD();
        finePts.Add(P0);
        SubdivideByLength(P0, P1, P2, fineResolution, finePts);
        finePts.Add(P2);
        
        // Decimate the fine points to meet the target resolution
        return DecimateToTargetResolution(finePts, maxSegLen);
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

        PathD decimated = new PathD();
        decimated.Add(finePoints[0]); // Always keep start point

        // Calculate total curve length to determine if we need special handling for small curves
        double totalLength = 0;
        for (int i = 1; i < finePoints.Count; i++)
        {
            totalLength += Length(Minus(finePoints[i], finePoints[i-1]));
        }
        
        // For very small curves relative to target resolution, ensure we keep some intermediate points
        int minIntermediatePoints = 1; // At least 1 intermediate point for any curve
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
                double distanceFromLast = Length(Minus(current, lastKept));
                
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

    static PointD Add(PointD a, PointD b) => new PointD(a.x + b.x, a.y + b.y);
    static PointD Minus(PointD a, PointD b) => new PointD(a.x - b.x, a.y - b.y);
    static PointD Mult(PointD a, double s) => new PointD(a.x * s, a.y * s);
    static double Length(PointD p) => Math.Sqrt(p.x * p.x + p.y * p.y);

    static PointD Normalized(PointD p)
    {
        double len = Length(p);
        return len > 0 ? new PointD(p.x / len, p.y / len) : new PointD(0, 0);
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

    static PointD Mid(PointD a, PointD b) => new PointD((a.x + b.x) / 2.0, (a.y + b.y) / 2.0);

    static PathD AssembleWithDiagonals(PathsD processedCorners, int[] corner_types)
    {
        int n = processedCorners.Count;
        PathD outPts = new PathD();

        // Mark short corners (only consider first n entries of corner_types)
        bool[] isShort = new bool[n];
        for (int k = 0; k < n; k++)
            isShort[k] = (k < corner_types.Length && corner_types[k] == (int)CornerType.ShortEdge);

        int i = 0;
        while (i < n)
        {
            if (!isShort[i])
            {
                // append corner polyline, avoiding duplicating a point
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

            PointD startPt = processedCorners[prevIdx].Last();
            PointD endPt = processedCorners[nextIdx].First();

            // Ensure startPt is present as last point to connect diagonal from
            if (outPts.Count == 0)
            {
                outPts.Add(startPt);
            }
            else
            {
                PointD last = outPts[^1];
                if (!(Math.Abs(last.x - startPt.x) < 1e-9 && Math.Abs(last.y - startPt.y) < 1e-9))
                    outPts.Add(startPt);
            }

            // Add diagonal by appending endPt (polyline will connect them as straight segment)
            outPts.Add(endPt);

            // Loop continues; nextIdx may be handled in subsequent iteration if not short
        }

        return outPts;
    }

    /// <summary>
    /// Assembles processed corner curves into a single continuous path.
    /// </summary>
    /// <param name="processedCorners">List of processed corner curves</param>
    /// <returns>Assembled continuous path</returns>
    static PathD AssembleProcessedCorners(PathsD processedCorners)
    {
        PathD assembled = new PathD();
        
        for (int i = 0; i < processedCorners.Count; i++)
        {
            PathD corner = processedCorners[i];
            if (corner.Count > 0)
            {
                if (assembled.Count == 0)
                {
                    assembled.AddRange(corner);
                }
                else
                {
                    // Avoid duplicating points at connections
                    PointD last = assembled[^1];
                    PointD first = corner[0];
                    if (Math.Abs(last.x - first.x) < 1e-9 && Math.Abs(last.y - first.y) < 1e-9)
                        assembled.AddRange(corner.Skip(1));
                    else
                        assembled.AddRange(corner);
                }
            }
        }
        
        return assembled;
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

    static string DrawLine(PointD a, PointD b, string c, double w) =>
        $"  <line x1=\"{a.x}\" y1=\"{-a.y}\" x2=\"{b.x}\" y2=\"{-b.y}\" stroke=\"{c}\" stroke-width=\"{w}\"/>";

    static string DrawCircle(PointD p, int r, string fill) =>
        $"  <circle cx=\"{p.x}\" cy=\"{-p.y}\" r=\"{r}\" fill=\"{fill}\"/>";

    static string DrawText(string txt, PointD p, int dx, int dy) =>
        $"  <text x=\"{p.x + dx}\" y=\"{-p.y + dy}\" font-size=\"10\" fill=\"#000\">{txt}</text>";
}
