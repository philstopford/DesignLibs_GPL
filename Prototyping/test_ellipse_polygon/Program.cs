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
        double angular_resolution = 1;

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
                radius = concave_radius;
            }

            PathD current_corner = processCorner(startLine, endLine, radius, angular_resolution, edge_resolution, SamplingMode.ByMaxSegmentLength);
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
    
    static int[] categorizeCorners(PathD path_)
    {
        // Need to do some prep work here.
        // Prepare status list
        int[] status = new int[path_.Count];
        
        // Remove our last point in case the polygon is explicitly closed - this avoids some trouble.
        PathD path = new(path_);
        bool trimmed_path = false;
        if (path[0].x == path[^1].x && path[0].y == path[^1].y)
        {
            trimmed_path = true;
            path.RemoveAt(path.Count-1);
        }
        
        // Determine polygon orientation: positive = CCW, negative = CW
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

        if (trimmed_path)
        {
            // Close path, first and last are the same
            status[^1] = status[0];
        }

        return status;
    }
    
    static PathD processCorner(PathD startLine, PathD endLine, double radius, double angular_resolution, double edge_resolution, SamplingMode mode = SamplingMode.ByMaxAngle)
    {

        // 1. Define base lines (scaled)
        PointD startLineStart = startLine[0];
        PointD startLineEnd   = startLine[1];
        PointD endLineStart   = endLine[0];
        PointD endLineEnd     = endLine[1];

        // 2. Compute curve start/end
        PointD startLength = Minus(startLineEnd, startLineStart);
        PointD startDir = Normalized(startLength);
        PointD endLength = Minus(endLineEnd, endLineStart);
        PointD endDir   = Normalized(endLength);

        // Set the radius for the curve each side.
        // Is the radius larger than our midpoint?

        double start_radius = radius;
        double half_edge_length = Math.Abs(Math.Sqrt(startLength.x * startLength.x + startLength.y * startLength.y));
        if (start_radius > half_edge_length)
        {
            start_radius = half_edge_length;
        }
        PointD curveStartPoint = Add(startLineStart, Mult(startDir, start_radius));

        double end_radius = radius;
        half_edge_length = Math.Abs(Math.Sqrt(endLength.x * endLength.x + endLength.y * endLength.y));
        if (end_radius > half_edge_length)
        {
            end_radius = half_edge_length;
        }
        PointD curveEndPoint   = Add(endLineStart, Mult(endDir, end_radius));

        // 3. Compute unique control point
        double dx = curveEndPoint.x - curveStartPoint.x;
        double dy = curveEndPoint.y - curveStartPoint.y;
        double det = startDir.x * endDir.y - startDir.y * endDir.x;
        if (Math.Abs(det) < 1e-6)
            throw new Exception("Tangent lines are parallel; no unique control point.");

        double tParam = (dx * endDir.y - dy * endDir.x) / det;
        PointD controlPoint = Add(curveStartPoint, Mult(startDir, tParam));

        // 4. Choose sampling mode
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

    // --- Sampling by max segment length ---
    static PathD SampleByMaxSegmentLength(
        PointD P0, PointD P1, PointD P2, double maxSegLen)
    {
        PathD pts = new PathD();
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
