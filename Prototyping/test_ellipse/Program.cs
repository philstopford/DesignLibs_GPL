using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;

class QuadraticBezierSamplingSwitcher
{
    enum SamplingMode
    {
        ByMaxSegmentLength,
        ByMaxAngle
    }

    static void Main()
    {
        const double scale = 100.0;

        // 1. Define base lines (scaled)
        var startLineStart = new Point(0, 0) * scale;
        var startLineEnd   = new Point(3, 15) * scale;
        var endLineStart   = new Point(0, 0) * scale;
        var endLineEnd     = new Point(20, 5) * scale;

        // 2. Compute curve start/end
        var startDir = (startLineEnd - startLineStart).Normalized();
        var endDir   = (endLineEnd   - endLineStart  ).Normalized();

        var curveStartPoint = startLineStart + startDir * (5 * scale);
        var curveEndPoint   = endLineStart   + endDir   * (7 * scale);

        // 3. Compute unique control point
        double dx = curveEndPoint.X - curveStartPoint.X;
        double dy = curveEndPoint.Y - curveStartPoint.Y;
        double det = startDir.X * endDir.Y - startDir.Y * endDir.X;
        if (Math.Abs(det) < 1e-6)
            throw new Exception("Tangent lines are parallel; no unique control point.");

        double tParam = (dx * endDir.Y - dy * endDir.X) / det;
        var controlPoint = curveStartPoint + startDir * tParam;

        // 4. Choose sampling mode
        SamplingMode mode = SamplingMode.ByMaxSegmentLength;
        List<Point> samples;
        switch (mode)
        {
            case SamplingMode.ByMaxSegmentLength:
                double maxSegLen = 0.05 * scale;
                samples = SampleByMaxSegmentLength(curveStartPoint, controlPoint, curveEndPoint, maxSegLen);
                break;

            case SamplingMode.ByMaxAngle:
                double maxAngleDeg = 5.0;
                double maxAngleRad = maxAngleDeg * Math.PI / 180.0;
                samples = SampleByMaxAngle(curveStartPoint, controlPoint, curveEndPoint, maxAngleRad);
                break;

            default:
                throw new ArgumentOutOfRangeException();
        }

        // 5. Export CSV & SVG
        WriteCsv("quad_curve_samples.csv", samples);
        Console.WriteLine($"<b>CSV saved ({samples.Count} pts)</b>");

        var keyPoints = new[] {
            startLineStart, startLineEnd,
            endLineStart,   endLineEnd,
            curveStartPoint, curveEndPoint, controlPoint
        };
        string svg = BuildDetailedSvg(keyPoints, samples);
        File.WriteAllText("quad_curve_samples.svg", svg, Encoding.UTF8);
        Console.WriteLine("<b>SVG saved</b>");
    }

    // --- Sampling by max segment length ---
    static List<Point> SampleByMaxSegmentLength(
        Point P0, Point P1, Point P2, double maxSegLen)
    {
        var pts = new List<Point> { P0 };
        SubdivideByLength(P0, P1, P2, maxSegLen, pts);
        pts.Add(P2);
        return pts;
    }

    static void SubdivideByLength(
        Point p0, Point p1, Point p2, double maxSegLen, List<Point> outPts)
    {
        if ((p2 - p0).Length() <= maxSegLen)
        {
            outPts.Add(p2);
            return;
        }
        Point p01 = Mid(p0, p1), p12 = Mid(p1, p2), p012 = Mid(p01, p12);
        SubdivideByLength(p0,  p01,  p012, maxSegLen, outPts);
        SubdivideByLength(p012, p12,  p2,   maxSegLen, outPts);
    }

    // --- Sampling by max angle ---
    static List<Point> SampleByMaxAngle(
        Point P0, Point P1, Point P2, double maxAngle)
    {
        var pts = new List<Point> { P0 };
        SubdivideByAngle(P0, P1, P2, maxAngle, pts);
        pts.Add(P2);
        return pts;
    }

    static void SubdivideByAngle(
        Point p0, Point p1, Point p2, double maxAngle, List<Point> outPts)
    {
        var tan0 = (p1 - p0).Normalized();
        var tan1 = (p2 - p1).Normalized();
        double dot = Math.Max(-1.0, Math.Min(1.0, tan0.X * tan1.X + tan0.Y * tan1.Y));
        double angle = Math.Acos(dot);
        if (angle <= maxAngle)
        {
            outPts.Add(p2);
            return;
        }
        Point p01 = Mid(p0, p1), p12 = Mid(p1, p2), p012 = Mid(p01, p12);
        SubdivideByAngle(p0,  p01,  p012, maxAngle, outPts);
        SubdivideByAngle(p012, p12,  p2,   maxAngle, outPts);
    }

    static Point Mid(Point a, Point b) =>
        new Point((a.X + b.X) / 2, (a.Y + b.Y) / 2);

    // --- CSV & SVG Utilities ---
    static void WriteCsv(string path, List<Point> pts)
    {
        var sb = new StringBuilder();
        sb.AppendLine("X,Y");
        foreach (var p in pts)
            sb.AppendLine(
              $"{p.X.ToString(CultureInfo.InvariantCulture)}," +
              $"{p.Y.ToString(CultureInfo.InvariantCulture)}");
        File.WriteAllText(path, sb.ToString());
    }

    static string BuildDetailedSvg(IEnumerable<Point> keyPts, List<Point> curvePts)
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
        foreach (var p in curvePts)
            sb.Append($"{p.X},{-p.Y} ");
        sb.AppendLine("\"/>");

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

    static double NiceStep(double r)
    {
        double e = Math.Floor(Math.Log10(r));
        double f = r / Math.Pow(10, e);
        double bs = f <= 2 ? 0.5 : f <= 5 ? 1 : 2;
        return bs * Math.Pow(10, e);
    }

    static string DrawLine(Point a, Point b, string c, double w) =>
        $"  <line x1=\"{a.X}\" y1=\"{-a.Y}\" x2=\"{b.X}\" y2=\"{-b.Y}\" stroke=\"{c}\" stroke-width=\"{w}\"/>";

    static string DrawCircle(Point p, int r, string fill) =>
        $"  <circle cx=\"{p.X}\" cy=\"{-p.Y}\" r=\"{r}\" fill=\"{fill}\"/>";

    static string DrawText(string txt, Point p, int dx, int dy) =>
        $"  <text x=\"{p.X+dx}\" y=\"{-p.Y+dy}\" font-size=\"10\" fill=\"#000\">{txt}</text>";

    struct Point
    {
        public double X, Y;
        public Point(double x, double y) { X = x; Y = y; }
        public static Point operator +(Point a, Point b) => new Point(a.X+b.X, a.Y+b.Y);
        public static Point operator -(Point a, Point b) => new Point(a.X-b.X, a.Y-b.Y);
        public static Point operator *(Point a, double s) => new Point(a.X*s, a.Y*s);
        public double Length() => Math.Sqrt(X*X + Y*Y);
        public Point Normalized()
        {
            double len = Length();
            return len > 0 ? new Point(X/len, Y/len) : new Point(0, 0);
        }
    }
}
