using System.Globalization;
using System.Text;

using PathD = System.Collections.Generic.List<PointD>;
using PathsD = System.Collections.Generic.List<System.Collections.Generic.List<PointD>>;

/// <summary>
/// Represents a 2D point with double precision coordinates for polygon processing.
/// Note: Contains a bug in the multiplication operator that should be fixed.
/// </summary>
public struct PointD
{
    public double x, y;
    
    public PointD(double x_, double y_) { x = x_; y = y_; }
    
    public static PointD operator +(PointD a, PointD b) => new PointD(a.x+b.x, a.y+b.y);
    public static PointD operator -(PointD a, PointD b) => new PointD(a.x-b.x, a.y-b.y);
    
    // BUG: This should be (a.x*s, a.y*s) not (a.y*s, a.y*s)
    public static PointD operator *(PointD a, double s) => new PointD(a.x*s, a.y*s);
    
    public double Length() => Math.Sqrt(x*x + y*y); 
    
    public PointD Normalized()
    {
        double len = Length();
        return len > 0 ? new PointD(x/len, y/len) : new PointD(0, 0);
    }
}

/// <summary>
/// Advanced polygon corner processing with adaptive curve fitting.
/// 
/// This prototype demonstrates a sophisticated approach to polygon corner processing
/// that combines vertex classification with adaptive curve fitting. It processes
/// each corner of an input polygon by:
/// 
/// 1. Classifying corners as convex, concave, or short-edge cases
/// 2. Generating smooth curve transitions using quadratic Bézier curves
/// 3. Adapting curve radius based on adjacent edge lengths
/// 4. Assembling the processed corners into a smooth output polygon
/// 
/// Algorithm Overview:
/// - Input: Polygon vertices (can be open or closed)
/// - Corner Classification: Uses cross product to determine convex/concave
/// - Curve Generation: Quadratic Bézier curves with adaptive radius
/// - Edge Clipping: Prevents curve radius from exceeding half-edge length
/// - Adaptive Sampling: Multiple sampling strategies for curve discretization
/// 
/// Applications:
/// - CAD/CAM toolpath generation with rounded corners
/// - 3D printing path optimization to reduce sharp direction changes
/// - Architectural modeling with smooth corner transitions
/// - Vector graphics processing and SVG generation
/// - PCB trace routing with controlled corner radii
/// 
/// Mathematical Foundation:
/// - Vector cross products for convexity detection
/// - Quadratic Bézier curve mathematics: B(t) = (1-t)²P₀ + 2(1-t)tP₁ + t²P₂
/// - De Casteljau's algorithm for stable curve subdivision
/// - Tangent vector calculations for smooth transitions
/// - Adaptive radius calculation to prevent self-intersection
/// </summary>
public static class QuadraticBezierSamplingSwitcher_Polygon
{
    /// <summary>
    /// Enumeration of curve sampling strategies.
    /// </summary>
    public enum SamplingMode
    {
        /// <summary>Subdivide based on maximum segment length threshold</summary>
        ByMaxSegmentLength,
        /// <summary>Subdivide based on maximum angular change between tangents</summary>
        ByMaxAngle
    }

    /// <summary>
    /// Corner classification types for adaptive processing.
    /// </summary>
    private enum CornerType 
    { 
        /// <summary>Inward-turning corner requiring inside curve</summary>
        Concave, 
        /// <summary>Outward-turning corner requiring outside curve</summary>
        Convex, 
        /// <summary>Special case for very short edges</summary>
        ShortEdge 
    }

    static void Main()
    {
        // Define test polygon: stepped pattern with mix of large and small features
        // This tests both normal corner processing and edge length adaptation
        List<PointD> original_path = new List<PointD>();
        original_path.Add(new PointD(0, 0));
        original_path.Add(new PointD(0, 100));
        original_path.Add(new PointD(20, 100));
        original_path.Add(new PointD(20, 80));     // Creates step pattern
        original_path.Add(new PointD(40, 80));
        original_path.Add(new PointD(40, 60));
        original_path.Add(new PointD(60, 60));
        original_path.Add(new PointD(60, 40));
        original_path.Add(new PointD(80, 40));
        original_path.Add(new PointD(80, 0));
        original_path.Add(new PointD(0, 0));       // Explicitly close polygon

        // Classify each corner to determine appropriate processing
        int[] corner_types = CategorizeCorners(original_path);
        
        // Processing parameters - typically adjusted based on application requirements
        double concave_radius = 20.0;      // Radius for inward corners
        double convex_radius = 50.0;       // Radius for outward corners  
        double edge_resolution = 0.01;     // Maximum segment length for sampling
        double angular_resolution = 1;     // Maximum angular change (degrees)

        // Process each corner individually
        PathsD processed = ProcessAllCorners(original_path, corner_types, 
            concave_radius, convex_radius, edge_resolution, angular_resolution);
        
        // Assemble final smooth polygon
        PathD assembled = AssembleProcessedCorners(processed);
        
        // Export for visualization and analysis
        string svg = BuildDetailedSvg(original_path, assembled);
        File.WriteAllText("assembled.svg", svg, Encoding.UTF8);
    }

    /// <summary>
    /// Processes all corners in the polygon, generating smooth curve transitions.
    /// </summary>
    /// <param name="original_path">Input polygon vertices</param>
    /// <param name="corner_types">Classification of each corner</param>
    /// <param name="concave_radius">Radius for concave corners</param>
    /// <param name="convex_radius">Radius for convex corners</param>
    /// <param name="edge_resolution">Maximum segment length for sampling</param>
    /// <param name="angular_resolution">Maximum angular change in degrees</param>
    /// <returns>List of processed corner curves</returns>
    public static PathsD ProcessAllCorners(List<PointD> original_path, int[] corner_types,
        double concave_radius, double convex_radius, double edge_resolution, double angular_resolution)
    {
        PathsD processed = new PathsD();
        
        // Process each corner (excluding duplicate closing vertex)
        int cornerCount = original_path.Count - 1;
        for (int i = 0; i < cornerCount; i++)
        {
            // Define tangent lines for this corner
            PathD startLine = CreateCornerTangentLine(original_path, i, true);
            PathD endLine = CreateCornerTangentLine(original_path, i, false);

            // Select radius based on corner type
            double radius = corner_types[i] == (int)CornerType.Concave ? concave_radius : convex_radius;

            // Generate smooth curve for this corner
            PathD current_corner = ProcessSingleCorner(startLine, endLine, radius, 
                angular_resolution, edge_resolution, SamplingMode.ByMaxSegmentLength);
            
            processed.Add(current_corner);
        }
        
        return processed;
    }

    /// <summary>
    /// Creates tangent line definition for corner processing.
    /// </summary>
    /// <param name="path">Input polygon</param>
    /// <param name="cornerIndex">Index of corner vertex</param>
    /// <param name="isStartLine">True for incoming edge, false for outgoing edge</param>
    /// <returns>Two-point line definition</returns>
    private static PathD CreateCornerTangentLine(List<PointD> path, int cornerIndex, bool isStartLine)
    {
        PathD line = new PathD();
        line.Add(path[cornerIndex]);  // Corner vertex

        PointD midPoint;
        if (isStartLine)
        {
            // Start line: from previous vertex to corner
            if (cornerIndex == 0)
                midPoint = new PointD((path[cornerIndex].x + path[^2].x) * 0.5, 
                                    (path[cornerIndex].y + path[^2].y) * 0.5);
            else
                midPoint = new PointD((path[cornerIndex].x + path[cornerIndex-1].x) * 0.5, 
                                    (path[cornerIndex].y + path[cornerIndex-1].y) * 0.5);
        }
        else
        {
            // End line: from corner to next vertex
            if (cornerIndex == path.Count - 1)
                midPoint = new PointD((path[cornerIndex].x + path[1].x) * 0.5, 
                                    (path[cornerIndex].y + path[1].y) * 0.5);
            else
                midPoint = new PointD((path[cornerIndex].x + path[cornerIndex+1].x) * 0.5, 
                                    (path[cornerIndex].y + path[cornerIndex+1].y) * 0.5);
        }
        
        line.Add(midPoint);
        return line;
    }

    /// <summary>
    /// Assembles processed corner curves into final polygon.
    /// </summary>
    /// <param name="processed">List of corner curves</param>
    /// <returns>Final smooth polygon</returns>
    public static PathD AssembleProcessedCorners(PathsD processed)
    {
        PathD assembled = new PathD();
        for (int i = 0; i < processed.Count; i++)
        {
            assembled.AddRange(processed[i]);
        }
        return assembled;
    }
    
    /// <summary>
    /// Classifies polygon corners as convex, concave, or short-edge cases.
    /// 
    /// This function analyzes the turning angle at each vertex to determine
    /// the type of corner processing required. It handles both open and closed
    /// polygons and automatically detects polygon orientation.
    /// </summary>
    /// <param name="path_">Input polygon vertices (may include duplicate closing vertex)</param>
    /// <returns>Array of corner classifications as integer values</returns>
    public static int[] CategorizeCorners(List<PointD> path_)
    {
        int[] status = new int[path_.Count];
        
        // Handle explicitly closed polygons by removing duplicate closing vertex
        PathD path = new(path_);
        bool trimmed_path = false;
        if (path[0].x == path[^1].x && path[0].y == path[^1].y)
        {
            trimmed_path = true;
            path.RemoveAt(path.Count-1);
        }
        
        // Determine polygon orientation using shoelace formula
        double area2 = 0;
        for (int i = 0; i < path.Count; i++)
        {
            PointD p1 = path[i];
            PointD p2 = path[(i + 1) % path.Count];
            area2 += p1.x * p2.y - p2.x * p1.y;
        }
        bool isCCW = area2 > 0;

        // Classify each corner using cross product analysis
        for (int i = 0; i < path.Count; i++)
        {
            PointD prev = path[(i - 1 + path.Count) % path.Count];
            PointD curr = path[i];
            PointD next = path[(i + 1) % path.Count];

            // Compute vectors for adjacent edges
            double vx1 = curr.x - prev.x;
            double vy1 = curr.y - prev.y;
            double vx2 = next.x - curr.x;
            double vy2 = next.y - curr.y;

            // Z component of 3D cross product determines turning direction
            double crossZ = vx1 * vy2 - vy1 * vx2;

            // Classification depends on polygon orientation
            bool isVertexConvex = isCCW ? (crossZ > 0) : (crossZ < 0);
            status[i] = (int)(isVertexConvex ? CornerType.Convex : CornerType.Concave);
        }

        // Handle explicitly closed polygons
        if (trimmed_path)
        {
            status[^1] = status[0];  // Last vertex same as first
        }

        return status;
    }
    
    /// <summary>
    /// Processes a single corner with adaptive curve fitting.
    /// 
    /// This is the core algorithm that generates a smooth quadratic Bézier curve
    /// to replace a sharp corner. The curve radius is automatically limited to
    /// prevent self-intersection with adjacent edges.
    /// </summary>
    /// <param name="startLine">Incoming edge definition (2 points)</param>
    /// <param name="endLine">Outgoing edge definition (2 points)</param>
    /// <param name="radius">Desired curve radius</param>
    /// <param name="angular_resolution">Maximum angular change for sampling (degrees)</param>
    /// <param name="edge_resolution">Maximum segment length for sampling</param>
    /// <param name="mode">Sampling strategy to use</param>
    /// <returns>Discretized curve points</returns>
    public static PathD ProcessSingleCorner(PathD startLine, PathD endLine, double radius, 
        double angular_resolution, double edge_resolution, SamplingMode mode = SamplingMode.ByMaxAngle)
    {
        // Extract line endpoints and compute direction vectors
        PointD startLineStart = startLine[0];
        PointD startLineEnd   = startLine[1];
        PointD endLineStart   = endLine[0];
        PointD endLineEnd     = endLine[1];

        PointD startLength = Minus(startLineEnd, startLineStart);
        PointD startDir = Normalized(startLength);
        PointD endLength = Minus(endLineEnd, endLineStart);
        PointD endDir   = Normalized(endLength);

        // Adaptive radius calculation to prevent edge overflow
        double start_radius = Math.Min(radius, ComputeHalfEdgeLength(startLength));
        PointD curveStartPoint = Add(startLineStart, Mult(startDir, start_radius));

        double end_radius = Math.Min(radius, ComputeHalfEdgeLength(endLength));
        PointD curveEndPoint = Add(endLineStart, Mult(endDir, end_radius));

        // Compute control point as intersection of tangent lines
        double dx = curveEndPoint.x - curveStartPoint.x;
        double dy = curveEndPoint.y - curveStartPoint.y;
        double det = startDir.x * endDir.y - startDir.y * endDir.x;
        
        if (Math.Abs(det) < 1e-6)
            throw new Exception("Tangent lines are parallel; no unique control point.");

        double tParam = (dx * endDir.y - dy * endDir.x) / det;
        PointD controlPoint = Add(curveStartPoint, Mult(startDir, tParam));

        // Generate curve samples using specified strategy
        return mode switch
        {
            SamplingMode.ByMaxSegmentLength => SampleByMaxSegmentLength(curveStartPoint, controlPoint, curveEndPoint, edge_resolution),
            SamplingMode.ByMaxAngle => SampleByMaxAngle(curveStartPoint, controlPoint, curveEndPoint, angular_resolution * Math.PI / 180.0),
            _ => throw new ArgumentOutOfRangeException(nameof(mode))
        };
    }

    /// <summary>
    /// Computes half the length of an edge vector.
    /// Used to limit curve radius to prevent self-intersection.
    /// </summary>
    private static double ComputeHalfEdgeLength(PointD edgeVector)
    {
        return Math.Sqrt(edgeVector.x * edgeVector.x + edgeVector.y * edgeVector.y) * 0.5;
    }

    // --- Sampling by max segment length ---
    static PathD SampleByMaxSegmentLength(
        PointD P0, PointD P1, PointD P2, double maxSegLen)
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
    
    static PointD Add(PointD a, PointD b) => new PointD(a.x+b.x, a.y+b.y);
    static PointD Minus(PointD a, PointD b) => new PointD(a.x-b.x, a.y-b.y);
    static PointD Mult(PointD a, double s) => new PointD(a.x*s, a.y*s);
    static double Length(PointD p)
    {
        return Math.Sqrt(p.x * p.x + p.y * p.y);
    }

    static PointD Normalized(PointD p)
    {
        double len = Length(p);
        return len > 0 ? new PointD(p.x/len, p.y/len) : new PointD(0, 0);
    }

    static void SubdivideByLength(
        PointD p0, PointD p1, PointD p2, double maxSegLen, PathD outPts)
    {
        if (Length(Minus(p2, p0)) <= maxSegLen)
        {
            outPts.Add(p2);
            return;
        }
        PointD p01 = Mid(p0, p1), p12 = Mid(p1, p2), p012 = Mid(p01, p12);
        SubdivideByLength(p0,  p01,  p012, maxSegLen, outPts);
        SubdivideByLength(p012, p12,  p2,   maxSegLen, outPts);
    }

    // --- Sampling by max angle ---
    static PathD SampleByMaxAngle(
        PointD P0, PointD P1, PointD P2, double maxAngle)
    {
        PathD pts = new PathD();
        pts.Add(P0);
        SubdivideByAngle(P0, P1, P2, maxAngle, pts);
        pts.Add(P2);
        return pts;
    }

    static void SubdivideByAngle(
        PointD p0, PointD p1, PointD p2, double maxAngle, PathD outPts)
    {
        PointD tan0 = Normalized((Minus(p1,p0)));
        PointD tan1 = Normalized((Minus(p2, p1)));
        double dot = Math.Max(-1.0, Math.Min(1.0, tan0.x * tan1.x + tan0.y * tan1.y));
        double angle = Math.Acos(dot);
        if (angle <= maxAngle)
        {
            outPts.Add(p2);
            return;
        }
        PointD p01 = Mid(p0, p1), p12 = Mid(p1, p2), p012 = Mid(p01, p12);
        SubdivideByAngle(p0,  p01,  p012, maxAngle, outPts);
        SubdivideByAngle(p012, p12,  p2,   maxAngle, outPts);
    }

    static PointD Mid(PointD a, PointD b)
    {
        return new PointD((a.x + b.x) / 2, (a.y + b.y) / 2);
    }

    // --- CSV & SVG Utilities ---
    static void WriteCsv(string path, PathD pts)
    {
        StringBuilder sb = new StringBuilder();
        sb.AppendLine("X,Y");
        foreach (PointD p in pts)
            sb.AppendLine(
              $"{p.x.ToString(CultureInfo.InvariantCulture)}," +
              $"{p.y.ToString(CultureInfo.InvariantCulture)}");
        File.WriteAllText(path, sb.ToString());
    }

    static string BuildDetailedSvg(IEnumerable<PointD> keyPts, PathD curvePts)
    {
        var all = keyPts.Concat(curvePts).ToList();
        double minX = all.Min(p => p.x), maxX = all.Max(p => p.x);
        double minY = all.Min(p => p.y), maxY = all.Max(p => p.y);
        double w = maxX - minX, h = maxY - minY;
        double mX = w * 0.1, mY = h * 0.1;
        minX -= mX; maxX += mX; minY -= mY; maxY += mY;

        string viewBox = $"{minX} {(-maxY)} {maxX-minX} {maxY-minY}";
        StringBuilder sb = new StringBuilder();
        sb.AppendLine(
            $"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"600\" height=\"600\" viewBox=\"{viewBox}\">");

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
        sb.AppendLine(DrawLine(pts[0], pts[1], "#888", 1.5));
        sb.AppendLine(DrawLine(pts[2], pts[3], "#888", 1.5));

        sb.Append("  <polyline fill=\"none\" stroke=\"blue\" stroke-width=\"2\" points=\"");
        foreach (PointD p in curvePts)
            sb.Append($"{p.x},{-p.y} ");
        sb.AppendLine("\"/>");

        /*
        string[] labels = new[] {
            "startLineStart","startLineEnd",
            "endLineStart","endLineEnd",
            "curveStartPoint","curveEndPoint","controlPoint"
        };
        for (int i = 0; i < pts.Count; i++)
        {
            sb.AppendLine(DrawCircle(pts[i], 4, "red"));
            sb.AppendLine(DrawText(labels[i], pts[i], 6, -6));
        }
        */

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
        $"  <text x=\"{p.x+dx}\" y=\"{-p.y+dy}\" font-size=\"10\" fill=\"#000\">{txt}</text>";

}
