using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;

class QuadraticBezierWithDetailedSvgScaled
{
    static void Main()
    {
        // 0. Uniform scale factor
        const double scale = 100.0;

        // 1. Define the two base lines (scaled)
        var startLineStart     = new Point(0, 0) * scale;
        var startLineEnd       = new Point(3, 15) * scale;
        var endLineStart       = new Point(0, 0) * scale;
        var endLineEnd         = new Point(20, 5) * scale;

        // 2. Compute curve start point at a scaled distance along the start line
        double distanceAlongStartLine = 5 * scale;
        var startDirection            = (startLineEnd - startLineStart).Normalized();
        var curveStartPoint           = startLineStart + startDirection * distanceAlongStartLine;

        // 3. Compute curve end point at a scaled distance along the end line
        double distanceAlongEndLine = 7 * scale;
        var endDirection             = (endLineEnd - endLineStart).Normalized();
        var curveEndPoint            = endLineStart + endDirection * distanceAlongEndLine;

        // 4. Compute the quadratic Bézier control point
        double deltaX      = curveEndPoint.X - curveStartPoint.X;
        double deltaY      = curveEndPoint.Y - curveStartPoint.Y;
        double determinant = startDirection.X * endDirection.Y
                           - startDirection.Y * endDirection.X;
        if (Math.Abs(determinant) < 1e-6)
            throw new Exception("Tangent lines are parallel; no unique control point.");

        double tParameter = (deltaX * endDirection.Y - deltaY * endDirection.X)
                          / determinant;
        var controlPoint = curveStartPoint + startDirection * tParameter;

        // 5. Sample the quadratic Bézier curve (still 'unitless' t-steps)
        int sampleCount = 100;
        var curveSamples = SampleQuadraticBezier(
            curveStartPoint, controlPoint, curveEndPoint, sampleCount
        );

        // 6. Export CSV (coordinates are already 100×)
        WriteCsv("quad_curve_points_scaled.csv", curveSamples);
        Console.WriteLine($"<b>CSV saved to quad_curve_points_scaled.csv ({curveSamples.Count} points)</b>");

        // 7. Build & write detailed SVG (all coords 100×)
        string svgContent = BuildDetailedSvg(
            new[] {
                startLineStart, startLineEnd,
                endLineStart,   endLineEnd,
                curveStartPoint, curveEndPoint, controlPoint
            },
            curveSamples
        );
        File.WriteAllText("quad_curve_scaled.svg", svgContent, Encoding.UTF8);
        Console.WriteLine("<b>SVG saved to quad_curve_scaled.svg</b>");
    }

    static List<Point> SampleQuadraticBezier(
        Point P0, Point P1, Point P2, int sampleCount)
    {
        var samples = new List<Point>(sampleCount + 1);
        for (int i = 0; i <= sampleCount; i++)
        {
            double t = i / (double)sampleCount;
            double u = 1 - t;
            double x = u*u*P0.X + 2*u*t*P1.X + t*t*P2.X;
            double y = u*u*P0.Y + 2*u*t*P1.Y + t*t*P2.Y;
            samples.Add(new Point(x, y));
        }
        return samples;
    }

    static void WriteCsv(string filePath, List<Point> points)
    {
        var sb = new StringBuilder();
        sb.AppendLine("X,Y");
        foreach (var pt in points)
            sb.AppendLine(
              $"{pt.X.ToString(CultureInfo.InvariantCulture)}," +
              $"{pt.Y.ToString(CultureInfo.InvariantCulture)}"
            );
        File.WriteAllText(filePath, sb.ToString());
    }

    static string BuildDetailedSvg(
        IEnumerable<Point> keyPts,
        List<Point> curvePts)
    {
        var allPoints = keyPts.Concat(curvePts).ToList();
        double minX = allPoints.Min(p => p.X), maxX = allPoints.Max(p => p.X);
        double minY = allPoints.Min(p => p.Y), maxY = allPoints.Max(p => p.Y);

        // 10% margin (scaled)
        double width  = maxX - minX,  height = maxY - minY;
        double marginX = width  * 0.1, marginY = height * 0.1;
        minX -= marginX; maxX += marginX;
        minY -= marginY; maxY += marginY;

        int svgW = 600, svgH = 600;
        string viewBox = $"{minX.ToString(CultureInfo.InvariantCulture)} " +
                         $"{(-maxY).ToString(CultureInfo.InvariantCulture)} " +
                         $"{(maxX - minX).ToString(CultureInfo.InvariantCulture)} " +
                         $"{(maxY - minY).ToString(CultureInfo.InvariantCulture)}";

        var sb = new StringBuilder();
        sb.AppendLine(
          $"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{svgW}\" height=\"{svgH}\" viewBox=\"{viewBox}\">"
        );

        // Grid
        double xStep = NiceStep(width), yStep = NiceStep(height);
        sb.AppendLine("  <g stroke=\"#ddd\" stroke-width=\"0.5\">");
        for (double x = Math.Ceiling(minX); x <= maxX; x += xStep)
            sb.AppendLine($"    <line x1=\"{x}\" y1=\"{-maxY}\" x2=\"{x}\" y2=\"{-minY}\"/>");
        for (double y = Math.Floor(minY); y <= maxY; y += yStep)
            sb.AppendLine($"    <line x1=\"{minX}\" y1=\"{-y}\" x2=\"{maxX}\" y2=\"{-y}\"/>");
        sb.AppendLine("  </g>");

        // Axes
        sb.AppendLine("  <g stroke=\"#000\" stroke-width=\"1\">");
        sb.AppendLine($"    <line x1=\"{minX}\" y1=\"0\" x2=\"{maxX}\" y2=\"0\"/>");
        sb.AppendLine($"    <line x1=\"0\" y1=\"{-maxY}\" x2=\"0\" y2=\"{-minY}\"/>");
        sb.AppendLine("  </g>");

        var pts = keyPts.ToList();
        // Base lines in gray
        sb.AppendLine(DrawLine(pts[0], pts[1], "#888", 1.5));
        sb.AppendLine(DrawLine(pts[2], pts[3], "#888", 1.5));

        // Quadratic curve in blue
        sb.Append("  <polyline fill=\"none\" stroke=\"blue\" stroke-width=\"2\" points=\"");
        foreach (var p in curvePts)
            sb.Append($"{p.X.ToString(CultureInfo.InvariantCulture)}," +
                      $"{(-p.Y).ToString(CultureInfo.InvariantCulture)} ");
        sb.AppendLine("\"/>");

        // Annotate key points
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

    static double NiceStep(double range)
    {
        double exp = Math.Floor(Math.Log10(range));
        double frac = range / Math.Pow(10, exp);
        double baseStep = frac <= 2 ? 0.5 : frac <= 5 ? 1 : 2;
        return baseStep * Math.Pow(10, exp);
    }

    static string DrawLine(Point p1, Point p2, string color, double w) =>
        $"  <line x1=\"{p1.X}\" y1=\"{-p1.Y}\" x2=\"{p2.X}\" y2=\"{-p2.Y}\" " +
        $"stroke=\"{color}\" stroke-width=\"{w}\"/>";

    static string DrawCircle(Point p, int r, string fill) =>
        $"  <circle cx=\"{p.X}\" cy=\"{-p.Y}\" r=\"{r}\" fill=\"{fill}\"/>";

    static string DrawText(string txt, Point p, int dx, int dy) =>
        $"  <text x=\"{p.X+dx}\" y=\"{-p.Y+dy}\" font-size=\"10\" fill=\"#000\">{txt}</text>";

    struct Point
    {
        public double X, Y;
        public Point(double x, double y) { X = x; Y = y; }
        public static Point operator +(Point a, Point b) =>
            new Point(a.X + b.X, a.Y + b.Y);
        public static Point operator -(Point a, Point b) =>
            new Point(a.X - b.X, a.Y - b.Y);
        public static Point operator *(Point a, double s) =>
            new Point(a.X * s, a.Y * s);
        public double Length() => Math.Sqrt(X*X + Y*Y);
        public Point Normalized()
        {
            var len = Length();
            return len > 0 ? new Point(X/len, Y/len) : new Point(0,0);
        }
    }
}
