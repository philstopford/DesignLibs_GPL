using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;

/// <summary>
/// Demonstrates two adaptive sampling strategies for quadratic Bézier curves.
/// 
/// This prototype implements and compares two different approaches for converting
/// smooth quadratic Bézier curves into discrete point sequences suitable for
/// rendering or further geometric processing.
/// 
/// Sampling Strategies:
/// 1. Maximum Segment Length: Recursively subdivides curve until all segments
///    are shorter than a specified threshold
/// 2. Maximum Angular Change: Subdivides based on tangent angle changes,
///    ensuring smooth curvature representation
/// 
/// Mathematical Foundation:
/// - Quadratic Bézier: B(t) = (1-t)²P₀ + 2(1-t)tP₁ + t²P₂
/// - De Casteljau subdivision for numerical stability
/// - Tangent vectors for angle-based sampling
/// 
/// Applications:
/// - Vector graphics rasterization
/// - CNC toolpath generation
/// - 3D printing path planning
/// - Font rendering
/// </summary>
public static class QuadraticBezierSamplingSwitcher
{
    /// <summary>
    /// Enumeration of available sampling modes for curve discretization.
    /// </summary>
    public enum SamplingMode
    {
        /// <summary>Subdivide until all segments are below maximum length threshold</summary>
        ByMaxSegmentLength,
        /// <summary>Subdivide based on maximum angular change between tangent vectors</summary>
        ByMaxAngle
    }

    static void Main()
    {
        const double scale = 100.0;

        // Define tangent lines that will intersect to form curve endpoints and control point
        var startLineStart = new Point(0, 0) * scale;      // Origin of first tangent line
        var startLineEnd   = new Point(3, 15) * scale;     // Direction point of first tangent
        var endLineStart   = new Point(0, 0) * scale;      // Origin of second tangent line  
        var endLineEnd     = new Point(20, 5) * scale;     // Direction point of second tangent

        // Extract normalized direction vectors from the tangent lines
        var startDir = (startLineEnd - startLineStart).Normalized();
        var endDir   = (endLineEnd   - endLineStart  ).Normalized();

        // Position curve endpoints along the tangent lines at specified distances
        var curveStartPoint = startLineStart + startDir * (5 * scale);
        var curveEndPoint   = endLineStart   + endDir   * (7 * scale);

        // Calculate control point as intersection of tangent lines
        var controlPoint = ComputeControlPoint(curveStartPoint, startDir, curveEndPoint, endDir);

        // Demonstrate different sampling strategies
        SamplingMode mode = SamplingMode.ByMaxSegmentLength;
        List<Point> samples = SampleCurve(curveStartPoint, controlPoint, curveEndPoint, mode, scale);

        // Export results for visualization and analysis
        ExportResults(samples, curveStartPoint, curveEndPoint, controlPoint, 
                     startLineStart, startLineEnd, endLineStart, endLineEnd);
    }

    /// <summary>
    /// Computes the control point for a quadratic Bézier curve given endpoints and tangent directions.
    /// </summary>
    /// <param name="start">Curve start point</param>
    /// <param name="startDir">Normalized tangent direction at start</param>
    /// <param name="end">Curve end point</param>
    /// <param name="endDir">Normalized tangent direction at end</param>
    /// <returns>Control point coordinates</returns>
    public static Point ComputeControlPoint(Point start, Point startDir, Point end, Point endDir)
    {
        double dx = end.X - start.X;
        double dy = end.Y - start.Y;
        double det = startDir.X * endDir.Y - startDir.Y * endDir.X;
        
        if (Math.Abs(det) < 1e-6)
            throw new Exception("Tangent lines are parallel; no unique control point.");

        double tParam = (dx * endDir.Y - dy * endDir.X) / det;
        return start + startDir * tParam;
    }

    /// <summary>
    /// Samples a quadratic Bézier curve using the specified sampling strategy.
    /// </summary>
    /// <param name="start">Curve start point</param>
    /// <param name="control">Curve control point</param>
    /// <param name="end">Curve end point</param>
    /// <param name="mode">Sampling strategy to use</param>
    /// <param name="scale">Scale factor for parameter adjustment</param>
    /// <returns>List of sampled points along the curve</returns>
    public static List<Point> SampleCurve(Point start, Point control, Point end, 
        SamplingMode mode, double scale)
    {
        return mode switch
        {
            SamplingMode.ByMaxSegmentLength => SampleByMaxSegmentLength(start, control, end, 0.05 * scale),
            SamplingMode.ByMaxAngle => SampleByMaxAngle(start, control, end, 5.0 * Math.PI / 180.0),
            _ => throw new ArgumentOutOfRangeException(nameof(mode))
        };
    }

    /// <summary>
    /// Samples curve by recursively subdividing until all segments are below length threshold.
    /// Uses De Casteljau's algorithm for numerically stable subdivision.
    /// </summary>
    /// <param name="P0">Start point</param>
    /// <param name="P1">Control point</param>
    /// <param name="P2">End point</param>
    /// <param name="maxSegLen">Maximum allowed segment length</param>
    /// <returns>List of points representing the curve</returns>
    public static List<Point> SampleByMaxSegmentLength(Point P0, Point P1, Point P2, double maxSegLen)
    {
        var pts = new List<Point> { P0 };
        SubdivideByLength(P0, P1, P2, maxSegLen, pts);
        pts.Add(P2);
        return pts;
    }

    /// <summary>
    /// Recursive helper for length-based subdivision using De Casteljau's algorithm.
    /// </summary>
    private static void SubdivideByLength(Point p0, Point p1, Point p2, double maxSegLen, List<Point> outPts)
    {
        // Check if current segment is small enough
        if ((p2 - p0).Length() <= maxSegLen)
        {
            outPts.Add(p2);
            return;
        }
        
        // Apply De Casteljau subdivision
        Point p01 = Mid(p0, p1), p12 = Mid(p1, p2), p012 = Mid(p01, p12);
        SubdivideByLength(p0,  p01,  p012, maxSegLen, outPts);
        SubdivideByLength(p012, p12,  p2,   maxSegLen, outPts);
    }

    /// <summary>
    /// Samples curve by subdividing based on angular changes between tangent vectors.
    /// Provides better curvature adaptation than length-based sampling.
    /// </summary>
    /// <param name="P0">Start point</param>
    /// <param name="P1">Control point</param>
    /// <param name="P2">End point</param>
    /// <param name="maxAngle">Maximum allowed angular change in radians</param>
    /// <returns>List of points representing the curve</returns>
    public static List<Point> SampleByMaxAngle(Point P0, Point P1, Point P2, double maxAngle)
    {
        var pts = new List<Point> { P0 };
        SubdivideByAngle(P0, P1, P2, maxAngle, pts);
        pts.Add(P2);
        return pts;
    }

    /// <summary>
    /// Recursive helper for angle-based subdivision.
    /// </summary>
    private static void SubdivideByAngle(Point p0, Point p1, Point p2, double maxAngle, List<Point> outPts)
    {
        // Calculate tangent vectors for angle comparison
        var tan0 = (p1 - p0).Normalized();
        var tan1 = (p2 - p1).Normalized();
        
        // Compute angle between tangents using dot product
        double dot = Math.Max(-1.0, Math.Min(1.0, tan0.X * tan1.X + tan0.Y * tan1.Y));
        double angle = Math.Acos(dot);
        
        if (angle <= maxAngle)
        {
            outPts.Add(p2);
            return;
        }
        
        // Subdivide using De Casteljau's algorithm
        Point p01 = Mid(p0, p1), p12 = Mid(p1, p2), p012 = Mid(p01, p12);
        SubdivideByAngle(p0,  p01,  p012, maxAngle, outPts);
        SubdivideByAngle(p012, p12,  p2,   maxAngle, outPts);
    }

    /// <summary>
    /// Calculates midpoint between two points (used in De Casteljau subdivision).
    /// </summary>
    private static Point Mid(Point a, Point b) => new Point((a.X + b.X) / 2, (a.Y + b.Y) / 2);

    /// <summary>
    /// Exports sampling results to CSV and SVG files for analysis and visualization.
    /// </summary>
    private static void ExportResults(List<Point> samples, Point curveStart, Point curveEnd, Point control,
        Point startLineStart, Point startLineEnd, Point endLineStart, Point endLineEnd)
    {
        WriteCsv("quad_curve_samples.csv", samples);
        Console.WriteLine($"<b>CSV saved ({samples.Count} pts)</b>");

        var keyPoints = new[] {
            startLineStart, startLineEnd,
            endLineStart,   endLineEnd,
            curveStart, curveEnd, control
        };
        string svg = BuildDetailedSvg(keyPoints, samples);
        File.WriteAllText("quad_curve_samples.svg", svg, Encoding.UTF8);
        Console.WriteLine("<b>SVG saved</b>");
    }

    /// <summary>
    /// Writes point data to CSV file for external analysis.
    /// </summary>
    public static void WriteCsv(string path, List<Point> pts)
    {
        var sb = new StringBuilder();
        sb.AppendLine("X,Y");
        foreach (var p in pts)
            sb.AppendLine(
              $"{p.X.ToString(CultureInfo.InvariantCulture)}," +
              $"{p.Y.ToString(CultureInfo.InvariantCulture)}");
        File.WriteAllText(path, sb.ToString());
    }

    /// <summary>
    /// Creates detailed SVG visualization showing construction lines, curve, and sample points.
    /// </summary>
    public static string BuildDetailedSvg(IEnumerable<Point> keyPts, List<Point> curvePts)
    {
        var all = keyPts.Concat(curvePts).ToList();
        double minX = all.Min(p => p.X), maxX = all.Max(p => p.X);
        double minY = all.Min(p => p.Y), maxY = all.Max(p => p.Y);
        double w = maxX - minX, h = maxY - minY;
        double mX = w * 0.1, mY = h * 0.1;
        minX -= mX; maxX += mX; minY -= mY; maxY += mY;

        string viewBox = $"{minX} {(-maxY)} {maxX-minX} {maxY-minY}";
        var sb = new StringBuilder();
        sb.AppendLine(
            $"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"600\" height=\"600\" viewBox=\"{viewBox}\">");

        // Add coordinate grid for reference
        double xs = NiceStep(w), ys = NiceStep(h);
        sb.AppendLine("  <g stroke=\"#ddd\" stroke-width=\"0.5\">");
        for (double x = Math.Ceiling(minX); x <= maxX; x += xs)
            sb.AppendLine($"    <line x1=\"{x}\" y1=\"{-maxY}\" x2=\"{x}\" y2=\"{-minY}\"/>");
        for (double y = Math.Floor(minY); y <= maxY; y += ys)
            sb.AppendLine($"    <line x1=\"{minX}\" y1=\"{-y}\" x2=\"{maxX}\" y2=\"{-y}\"/>");
        sb.AppendLine("  </g>");

        // Add coordinate axes
        sb.AppendLine("  <g stroke=\"#000\" stroke-width=\"1\">");
        sb.AppendLine($"    <line x1=\"{minX}\" y1=\"0\" x2=\"{maxX}\" y2=\"0\"/>");
        sb.AppendLine($"    <line x1=\"0\" y1=\"{-maxY}\" x2=\"0\" y2=\"{-minY}\"/>");
        sb.AppendLine("  </g>");

        // Draw construction lines and curve visualization
        var pts = keyPts.ToList();
        sb.AppendLine(DrawLine(pts[0], pts[1], "#888", 1.5));  // Start tangent line
        sb.AppendLine(DrawLine(pts[2], pts[3], "#888", 1.5));  // End tangent line

        // Draw sampled curve as polyline
        sb.Append("  <polyline fill=\"none\" stroke=\"blue\" stroke-width=\"2\" points=\"");
        foreach (var p in curvePts)
            sb.Append($"{p.X},{-p.Y} ");
        sb.AppendLine("\"/>");

        // Mark key points with labels
        var labels = new[] {
            "startLineStart","startLineEnd",
            "endLineStart","endLineEnd",
            "curveStartPoint","curveEndPoint","controlPoint"
        };
        for (int i = 0; i < pts.Count; i++)
        {
            sb.AppendLine(DrawCircle(pts[i], 4, "red"));
            sb.AppendLine(DrawText(labels[i], pts[i], 6, -6));
        }

        sb.AppendLine("</svg>");
        return sb.ToString();
    }

    /// <summary>Computes appropriate grid step size for visualization</summary>
    private static double NiceStep(double r)
    {
        double e = Math.Floor(Math.Log10(r));
        double f = r / Math.Pow(10, e);
        double bs = f <= 2 ? 0.5 : f <= 5 ? 1 : 2;
        return bs * Math.Pow(10, e);
    }

    /// <summary>SVG helper for drawing lines</summary>
    private static string DrawLine(Point a, Point b, string c, double w) =>
        $"  <line x1=\"{a.X}\" y1=\"{-a.Y}\" x2=\"{b.X}\" y2=\"{-b.Y}\" stroke=\"{c}\" stroke-width=\"{w}\"/>";

    /// <summary>SVG helper for drawing circles</summary>
    private static string DrawCircle(Point p, int r, string fill) =>
        $"  <circle cx=\"{p.X}\" cy=\"{-p.Y}\" r=\"{r}\" fill=\"{fill}\"/>";

    /// <summary>SVG helper for drawing text labels</summary>
    private static string DrawText(string txt, Point p, int dx, int dy) =>
        $"  <text x=\"{p.X+dx}\" y=\"{-p.Y+dy}\" font-size=\"10\" fill=\"#000\">{txt}</text>";

    /// <summary>
    /// Represents a 2D point with vector operations for geometric calculations.
    /// </summary>
    public struct Point
    {
        public double X, Y;
        
        public Point(double x, double y) { X = x; Y = y; }
        
        public static Point operator +(Point a, Point b) => new Point(a.X+b.X, a.Y+b.Y);
        public static Point operator -(Point a, Point b) => new Point(a.X-b.X, a.Y-b.Y);
        public static Point operator *(Point a, double s) => new Point(a.X*s, a.Y*s);
        
        /// <summary>Calculates Euclidean distance from origin</summary>
        public double Length() => Math.Sqrt(X*X + Y*Y);
        
        /// <summary>Returns unit vector in same direction</summary>
        public Point Normalized()
        {
            double len = Length();
            return len > 0 ? new Point(X/len, Y/len) : new Point(0, 0);
        }
    }
}
