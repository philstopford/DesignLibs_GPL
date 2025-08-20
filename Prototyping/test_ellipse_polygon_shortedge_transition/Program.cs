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

    // Assemble, replacing runs of short-edge corners with diagonals
    PathD assembled = AssembleWithDiagonals(processed, corner_types, cornerMidpoints, cornerVertices,
        short_edge_length, max_short_edge_length, smoothMode);

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

static PathD AssembleWithDiagonals(PathsD processedCorners, int[] corner_types, List<PointD> cornerMidpoints, List<PointD> cornerVertices,
    double short_edge_length, double max_short_edge_length, SmoothStepMode smoothMode)
{
    int n = processedCorners.Count;
    if (cornerMidpoints == null || cornerMidpoints.Count != n)
        throw new ArgumentException("cornerMidpoints must be same length as processedCorners");
    if (cornerVertices == null || cornerVertices.Count != n)
        throw new ArgumentException("cornerVertices must be same length as processedCorners");

    PathD outPts = new PathD();

    // helper smoothstep functions
    static double CubicSmoothstep(double t) => t <= 0 ? 0 : t >= 1 ? 1 : (3.0 * t * t - 2.0 * t * t * t);
    static double QuinticSmoothstep(double t) => t <= 0 ? 0 : t >= 1 ? 1 : (6.0 * Math.Pow(t, 5) - 15.0 * Math.Pow(t, 4) + 10.0 * Math.Pow(t, 3));
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
        PointD diagEndMid = cornerMidpoints[nextIdx];  // midpoint on edge leaving nextIdx

        // Corner vertex positions adjacent to the short run:
        PointD vertexPrev = cornerVertices[prevIdx];
        PointD vertexNext = cornerVertices[nextIdx];

        // Distances from vertex to its adjacent midpoint (per-side lengths)
        double dPrev = Length(Minus(diagStartMid, vertexPrev)); // distance from prev vertex to its outgoing midpoint
        double dNext = Length(Minus(diagEndMid, vertexNext));  // distance from next vertex to its outgoing midpoint

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
            outPts.Add(processedEndPt);
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
        PointD diagControl = new PointD((processedStartPt.x + processedEndPt.x) / 2.0, (processedStartPt.y + processedEndPt.y) / 2.0);
        PointD bezierControl = new PointD((bezierStartEndpoint.x + bezierEndEndpoint.x) / 2.0, (bezierStartEndpoint.y + bezierEndEndpoint.y) / 2.0);
        double avgBlend = (blendPrev + blendNext) * 0.5;
        PointD blendedControl = new PointD(
            diagControl.x * (1 - avgBlend) + bezierControl.x * avgBlend,
            diagControl.y * (1 - avgBlend) + bezierControl.y * avgBlend
        );

        // Sample curve between blendedStart -> blendedControl -> blendedEnd
        PathD blendedCurve = SampleByMaxSegmentLength(blendedStart, blendedControl, blendedEnd, 0.5);
        outPts.AddRange(blendedCurve);
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

    static PointD Add(PointD a, PointD b) => new PointD(a.x + b.x, a.y + b.y);
    static PointD Minus(PointD a, PointD b) => new PointD(a.x - b.x, a.y - b.y);
    static PointD Mult(PointD a, double s) => new PointD(a.x * s, a.y * s);
    static double Length(PointD p) => Math.Sqrt(p.x * p.x + p.y * p.y);

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
    public double Length() => Math.Sqrt(x * x + y * y);

    public PointD Normalized()
    {
        double len = Length();
        return len > 0 ? new PointD(x / len, y / len) : new PointD(0, 0);
    }
}