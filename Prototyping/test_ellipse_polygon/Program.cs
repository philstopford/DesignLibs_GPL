using System.Globalization;
using System.Text;
using Clipper2Lib;


class QuadraticBezierSamplingSwitcher_Polygon
{
    enum SamplingMode
    {
        ByMaxSegmentLength,
        ByMaxAngle
    }

    static void Main()
    {
        // Lowest X, lowest Y point first. Closed polygon, clockwise oriented for now.
        // Start simple.
        PathD original_path = [];
        original_path.Add(new(0,0));
        original_path.Add(new(0,50));
        original_path.Add(new(50,50));
        original_path.Add(new(50,0));
        original_path.Add(new(0,0));

        // Get our corner types
        int[] corner_types = categorizeCorners(original_path);
        
        double concave_radius = 20.0;
        double convex_radius = 50.0;
        double edge_resolution = 0.01;

        PathsD processed = new PathsD();
        // Need to now walk the corners.
        // Closed path expected.
        for (int i = 0; i < original_path.Count - 1; i++)
        {
            // Need to define our lines for the corner.
            // Start line...
            PathD startLine = new();
            startLine.Add(original_path[i]);
            PathD endLine = new();
            endLine.Add(original_path[i]);

            // First and last points are the same, so need to modify in case of first and last point.
            PointD midPoint;
            if (i == 0)
            {
                midPoint = new ((original_path[i].x + original_path[^2].x) * 0.5, (original_path[i].y + original_path[^2].y) * 0.5);
            }
            else
            {
                midPoint = new ((original_path[i].x + original_path[i-1].x) * 0.5, (original_path[i].y + original_path[i-1].y) * 0.5);
            }
            startLine.Add(midPoint);
            if (i == original_path.Count - 1)
            {
                midPoint = new ((original_path[i].x + original_path[2].x) * 0.5, (original_path[i].y + original_path[2].y) * 0.5);
            }
            else
            {
                midPoint = new ((original_path[i].x + original_path[i+1].x) * 0.5, (original_path[i].y + original_path[i+1].y) * 0.5);
            }
            endLine.Add(midPoint);
            
            // What's the nature of our corner?
            double radius = convex_radius;
            if (corner_types[i] == (int)types.concave)
            {
                radius = convex_radius;
            }

            PathD current_corner = processCorner(startLine, endLine, radius);
            processed.Add(current_corner);
        }
        
        // Processed now contains a list of corners that can be assembled.
        // Hopefully.
        PathD assembled = new();
        for (int i = 0; i < processed.Count; i++)
        {
            assembled.AddRange(processed[i]);
        }
        
        string svg = BuildDetailedSvg(original_path, assembled);
        File.WriteAllText("assembled.svg", svg, Encoding.UTF8);
    }
    
    private enum types {concave, convex, shortedge}
    
    static int[] categorizeCorners(PathD path)
    {
        // Determine polygon orientation: positive = CCW, negative = CW
        double area2 = 0;
        for (int i = 0; i < path.Count; i++)
        {
            var p1 = path[i];
            var p2 = path[(i + 1) % path.Count];
            area2 += p1.x * p2.y - p2.x * p1.y;
        }
        bool isCCW = area2 > 0;

        // Prepare status list
        var status = new int[path.Count];

        for (int i = 0; i < path.Count - 1; i++)
        {
            var prev = path[(i - 1 + path.Count) % path.Count];
            var curr = path[i];
            var next = path[(i + 1) % path.Count];

            // Vectors: prev→curr and curr→next
            double vx1 = curr.x - prev.x;
            double vy1 = curr.y - prev.y;
            double vx2 = next.x - curr.x;
            double vy2 = next.y - curr.y;

            // Z component of 3D cross product
            double crossZ = vx1 * vy2 - vy1 * vx2;

            // For CCW polygon, positive crossZ = convex. For CW, negative = convex.
            bool isVertexConvex = isCCW ? (crossZ > 0) : (crossZ < 0);
            status[i] = (int)(isVertexConvex ? types.convex : types.concave);
        }
        
        // Close path, first and last are the same
        status[path.Count - 1] = status[0];

        return status;
    }
    
    static PathD processCorner(PathD startLine, PathD endLine, double radius)
    {

        // 1. Define base lines (scaled)
        var startLineStart = startLine[0];
        var startLineEnd   = startLine[1];
        var endLineStart   = endLine[0];
        var endLineEnd     = endLine[1];

        // 2. Compute curve start/end
        var startLength = Minus(startLineEnd, startLineStart);
        var startDir = Normalized(startLength);
        var endLength = Minus(endLineEnd, endLineStart);
        var endDir   = Normalized(endLength);

        // Set the radius for the curve each side.
        // Is the radius larger than our midpoint?

        double start_radius = radius;
        double half_edge_length = Math.Abs(Math.Sqrt(startLength.x * startLength.x + startLength.y * startLength.y) * 0.5);
        if (start_radius > half_edge_length)
        {
            start_radius = half_edge_length;
        }
        var curveStartPoint = Add(startLineStart, Mult(startDir, start_radius));

        double end_radius = radius;
        half_edge_length = Math.Abs(Math.Sqrt(endLength.x * endLength.x + endLength.y * endLength.y) * 0.5);
        if (end_radius > half_edge_length)
        {
            end_radius = half_edge_length;
        }
        var curveEndPoint   = Add(endLineStart, Mult(endDir, end_radius));

        // 3. Compute unique control point
        double dx = curveEndPoint.x - curveStartPoint.x;
        double dy = curveEndPoint.y - curveStartPoint.y;
        double det = startDir.x * endDir.y - startDir.y * endDir.x;
        if (Math.Abs(det) < 1e-6)
            throw new Exception("Tangent lines are parallel; no unique control point.");

        double tParam = (dx * endDir.y - dy * endDir.x) / det;
        var controlPoint = Add(curveStartPoint, Mult(startDir, tParam));

        // 4. Choose sampling mode
        SamplingMode mode = SamplingMode.ByMaxSegmentLength;
        PathD samples;
        switch (mode)
        {
            case SamplingMode.ByMaxSegmentLength:
                double maxSegLen = 0.05;
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

        return samples;

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
    static PathD SampleByMaxSegmentLength(
        PointD P0, PointD P1, PointD P2, double maxSegLen)
    {
        var pts = new PathD();
        pts.Add(P0);
        SubdivideByLength(P0, P1, P2, maxSegLen, pts);
        pts.Add(P2);
        return pts;
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
        var pts = new PathD();
        pts.Add(P0);
        SubdivideByAngle(P0, P1, P2, maxAngle, pts);
        pts.Add(P2);
        return pts;
    }

    static void SubdivideByAngle(
        PointD p0, PointD p1, PointD p2, double maxAngle, PathD outPts)
    {
        var tan0 = Normalized((Minus(p1,p0)));
        var tan1 = Normalized((Minus(p2, p1)));
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
        var sb = new StringBuilder();
        sb.AppendLine("X,Y");
        foreach (var p in pts)
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
            sb.Append($"{p.x},{-p.y} ");
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

    static string DrawLine(PointD a, PointD b, string c, double w) =>
        $"  <line x1=\"{a.x}\" y1=\"{-a.y}\" x2=\"{b.x}\" y2=\"{-b.y}\" stroke=\"{c}\" stroke-width=\"{w}\"/>";

    static string DrawCircle(PointD p, int r, string fill) =>
        $"  <circle cx=\"{p.x}\" cy=\"{-p.y}\" r=\"{r}\" fill=\"{fill}\"/>";

    static string DrawText(string txt, PointD p, int dx, int dy) =>
        $"  <text x=\"{p.x+dx}\" y=\"{-p.y+dy}\" font-size=\"10\" fill=\"#000\">{txt}</text>";

    /*
    struct Point_
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
    */
}
