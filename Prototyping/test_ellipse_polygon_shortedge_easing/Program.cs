using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;

using PathD = System.Collections.Generic.List<PointD>;
using PathsD = System.Collections.Generic.List<System.Collections.Generic.List<PointD>>;

struct PointD
{
    public double x, y;
    public PointD(double x_, double y_) { x = x_; y = y_; }
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

class QuadraticBezierSamplingSwitcher_Polygon_Easing
{
    enum SamplingMode { ByMaxSegmentLength, ByMaxAngle }
    private enum CornerType { concave, convex, shortedge }

    enum EasingMode { None, QuadraticEase, CubicHermite }

    static void Main()
    {
        // --- Configurable input & parameters ---
        List<PointD> original_path = new List<PointD> {
            new PointD(0,0),
            new PointD(0,100),
            new PointD(20,100),
            new PointD(20,80),
            new PointD(40,80),
            new PointD(40,60),
            new PointD(60,60),
            new PointD(60,40),
            new PointD(80,40),
            new PointD(80,0),
            new PointD(0,0)
        };

        double concave_radius = 20.0;
        double convex_radius = 50.0;
        double edge_resolution = 0.5;
        double angular_resolution = 1.0;
        double short_edge_length = 20.0;

        // Easing configuration:
        EasingMode easingMode = EasingMode.CubicHermite; // None | QuadraticEase | CubicHermite
        double insetFraction = 0.12;   // fraction of diagonal length to inset at each end (0..0.5)
        double minInset = 2.0;         // minimum inset in world units
        double maxInset = 40.0;        // maximum inset
        int diagStraightSample = 0;    // samples for straight center section (0 = just endpoints)

        // --- Processing ---
        int[] corner_types = categorizeCorners(original_path, short_edge_length);

        PathsD processed = new PathsD();
        for (int i = 0; i < original_path.Count - 1; i++)
        {
            PathD startLine = new PathD { original_path[i] };
            PathD endLine = new PathD { original_path[i] };

            PointD midPrev = (i == 0) ? new PointD((original_path[i].x + original_path[^2].x) * 0.5, (original_path[i].y + original_path[^2].y) * 0.5)
                                      : new PointD((original_path[i].x + original_path[i - 1].x) * 0.5, (original_path[i].y + original_path[i - 1].y) * 0.5);
            startLine.Add(midPrev);

            PointD midNext = (i == original_path.Count - 2) ? new PointD((original_path[i].x + original_path[0].x) * 0.5, (original_path[i].y + original_path[0].y) * 0.5)
                                                            : new PointD((original_path[i].x + original_path[i + 1].x) * 0.5, (original_path[i].y + original_path[i + 1].y) * 0.5);
            endLine.Add(midNext);

            double radius = convex_radius;
            if (corner_types[i] == (int)CornerType.concave) radius = concave_radius;

            processed.Add(processCorner(startLine, endLine, radius, angular_resolution, edge_resolution));
        }

        PathD assembled = AssembleWithEasing(processed, corner_types, easingMode, insetFraction, minInset, maxInset, diagStraightSample);

        string svg = BuildDetailedSvg(original_path, assembled);
        File.WriteAllText("assembled_with_easing.svg", svg, Encoding.UTF8);
        Console.WriteLine("Wrote assembled_with_easing.svg");
    }

    // --- categorizeCorners (unchanged logic) ---
    static int[] categorizeCorners(PathD path_, double short_edge_length)
    {
        int[] status = new int[path_.Count];
        PathD path = new PathD(path_);
        bool trimmed = false;
        if (path.Count > 1 && path[0].x == path[^1].x && path[0].y == path[^1].y)
        {
            trimmed = true;
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

            double len1 = Math.Sqrt(vx1 * vx1 + vy1 * vy1);
            double len2 = Math.Sqrt(vx2 * vx2 + vy2 * vy2);

            if (len1 <= short_edge_length && len2 <= short_edge_length)
            {
                status[i] = (int)CornerType.shortedge;
                continue;
            }

            double crossZ = vx1 * vy2 - vy1 * vx2;
            bool isVertexConvex = isCCW ? (crossZ > 0) : (crossZ < 0);
            status[i] = (int)(isVertexConvex ? CornerType.convex : CornerType.concave);
        }

        if (trimmed)
        {
            Array.Resize(ref status, status.Length + 1);
            status[^1] = status[0];
        }

        return status;
    }

    // --- corner processing / sampling (unchanged) ---
    static PathD processCorner(PathD startLine, PathD endLine, double radius, double angular_resolution, double edge_resolution, SamplingMode mode = SamplingMode.ByMaxSegmentLength)
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
            PathD straight = new PathD { curveStartPoint, curveEndPoint };
            return straight;
        }

        double tParam = (dx * endDir.y - dy * endDir.x) / det;
        PointD controlPoint = Add(curveStartPoint, Mult(startDir, tParam));

        // use fixed mode ByMaxSegmentLength for now
        return SampleByMaxSegmentLength(curveStartPoint, controlPoint, curveEndPoint, edge_resolution);
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
    
    // --- Easing assembly ---
    static PathD AssembleWithEasing(
        PathsD processedCorners,
        int[] corner_types,
        EasingMode easingMode,
        double insetFraction,
        double minInset,
        double maxInset,
        int diagStraightSample)
    {
        int n = processedCorners.Count;
        PathD outPts = new PathD();

        bool[] isShort = new bool[n];
        for (int i = 0; i < n; i++) isShort[i] = (i < corner_types.Length && corner_types[i] == (int)CornerType.shortedge);

        int idx = 0;
        while (idx < n)
        {
            if (!isShort[idx])
            {
                // append corner polyline avoiding duplicate first point
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

            // found short run
            int runStart = idx;
            while (idx < n && isShort[idx]) idx++;
            int runEnd = idx - 1;

            // find prev non-short index
            int prevIdx = (runStart - 1 + n) % n;
            while (isShort[prevIdx]) prevIdx = (prevIdx - 1 + n) % n;

            // find next non-short index
            int nextIdx = (runEnd + 1) % n;
            while (isShort[nextIdx]) nextIdx = (nextIdx + 1) % n;

            PointD startPt = processedCorners[prevIdx].Last();
            PointD endPt = processedCorners[nextIdx].First();

            // Ensure startPt present
            if (outPts.Count == 0) outPts.Add(startPt);
            else
            {
                PointD last = outPts[^1];
                if (!(Math.Abs(last.x - startPt.x) < 1e-9 && Math.Abs(last.y - startPt.y) < 1e-9))
                    outPts.Add(startPt);
            }

            // compute diagonal vector and inset
            PointD diag = Minus(endPt, startPt);
            double diagLen = Length(diag);
            if (diagLen <= 1e-9)
            {
                // degenerate: just append end
                outPts.Add(endPt);
                continue;
            }
            PointD diagDir = Normalized(diag);
            double inset = Math.Max(minInset, Math.Min(maxInset, diagLen * insetFraction));
            if (inset * 2.0 > diagLen) inset = diagLen * 0.5 * 0.999; // keep positive center

            PointD S = startPt + diagDir * inset; // inward from start
            PointD E = endPt - diagDir * inset;   // inward from end

            if (easingMode == EasingMode.None)
            {
                // direct diagonal; optionally sample center
                if (diagStraightSample <= 0)
                {
                    outPts.Add(endPt);
                }
                else
                {
                    // split into diagStraightSample equal segments (excluding endpoints already present)
                    for (int s = 1; s <= diagStraightSample; s++)
                    {
                        double t = (double)s / (diagStraightSample + 1);
                        PointD p = new PointD(startPt.x + diag.x * t, startPt.y + diag.y * t);
                        outPts.Add(p);
                    }
                    outPts.Add(endPt);
                }
                continue;
            }

            // we will insert easing segments at both ends and a (possibly sampled) straight center between S and E
            // estimate tangents from processed polylines to blend with
            PointD prevTangent = EstimateOutgoingTangent(processedCorners[prevIdx]);
            PointD nextTangent = EstimateIncomingTangent(processedCorners[nextIdx]);

            if (easingMode == EasingMode.QuadraticEase)
            {
                // Build small quadratic from startPt -> S, control near a blend of diagDir and prevTangent
                PointD startControl = startPt + (diagDir * (inset * 0.5)) + (prevTangent * (inset * 0.25));
                PathD startQuad = SampleQuadratic(startPt, startControl, S, Math.Min(0.5, inset * 0.25)); // small sample length

                // Build small quadratic from E -> endPt
                PointD endControl = endPt - (diagDir * (inset * 0.5)) - (nextTangent * (inset * 0.25));
                PathD endQuad = SampleQuadratic(E, endControl, endPt, Math.Min(0.5, inset * 0.25));

                // append startQuad (avoid duplicate of startPt)
                if (outPts.Count > 0 && PointsEqual(outPts[^1], startQuad[0])) outPts.AddRange(startQuad.Skip(1));
                else outPts.AddRange(startQuad);

                // center straight S->E
                if (diagStraightSample <= 0)
                {
                    outPts.Add(E);
                }
                else
                {
                    for (int s = 1; s <= diagStraightSample; s++)
                    {
                        double t = (double)s / (diagStraightSample + 1);
                        PointD p = S + (diagDir * (diagLen - 2 * inset) * t);
                        outPts.Add(p);
                    }
                    outPts.Add(E);
                }

                // append endQuad (avoid duplicate E)
                if (PointsEqual(outPts[^1], endQuad[0])) outPts.AddRange(endQuad.Skip(1));
                else outPts.AddRange(endQuad);
            }
            else if (easingMode == EasingMode.CubicHermite)
            {
                // Cubic Hermite-style easing segments:
                // start cubic: startPt -> S, H0 tangent along diagDir, H1 tangent = prevTangent (scaled)
                double hLen = inset * 0.6; // handle length scale; tweak as desired
                PointD m0 = diagDir * hLen;
                PointD m1 = prevTangent * hLen;

                // convert Hermite (p0,p1,m0,m1) to cubic Bezier (p0, b1, b2, p1)
                PathD startCubic = HermiteToCubicAndSample(startPt, S, m0, m1, Math.Min(0.5, inset * 0.25));

                // end cubic: E -> endPt, tangent at E = -diagDir * hLen, tangent at endPt = nextTangent * hLen
                PointD m0e = (-1.0 * diagDir) * hLen;
                PointD m1e = nextTangent * hLen;
                PathD endCubic = HermiteToCubicAndSample(E, endPt, m0e, m1e, Math.Min(0.5, inset * 0.25));

                // append start cubic
                if (outPts.Count > 0 && PointsEqual(outPts[^1], startCubic[0])) outPts.AddRange(startCubic.Skip(1));
                else outPts.AddRange(startCubic);

                // center straight
                if (diagStraightSample <= 0)
                {
                    outPts.Add(E);
                }
                else
                {
                    for (int s = 1; s <= diagStraightSample; s++)
                    {
                        double t = (double)s / (diagStraightSample + 1);
                        PointD p = S + (diagDir * (diagLen - 2 * inset) * t);
                        outPts.Add(p);
                    }
                    outPts.Add(E);
                }

                // append end cubic
                if (PointsEqual(outPts[^1], endCubic[0])) outPts.AddRange(endCubic.Skip(1));
                else outPts.AddRange(endCubic);
            }
        }

        return outPts;
    }

    static bool PointsEqual(PointD a, PointD b)
    {
        return Math.Abs(a.x - b.x) < 1e-9 && Math.Abs(a.y - b.y) < 1e-9;
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

    // --- small math & svg utilities ---
    static PointD Mid(PointD a, PointD b) => new PointD((a.x + b.x) / 2.0, (a.y + b.y) / 2.0);
    static PointD Add(PointD a, PointD b) => new PointD(a.x + b.x, a.y + b.y);
    static PointD Minus(PointD a, PointD b) => new PointD(a.x - b.x, a.y - b.y);
    static PointD Mult(PointD a, double s) => new PointD(a.x * s, a.y * s);
    static double Length(PointD p) => Math.Sqrt(p.x * p.x + p.y * p.y);
    static PointD Normalized(PointD p)
    {
        double l = Length(p);
        return l > 0 ? new PointD(p.x / l, p.y / l) : new PointD(0, 0);
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
        sb.AppendLine($"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"800\" height=\"800\" viewBox=\"{viewBox}\">");

        // grid
        double xs = NiceStep(w), ys = NiceStep(h);
        sb.AppendLine("  <g stroke=\"#eee\" stroke-width=\"0.5\">");
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
        // draw original polygon key segments (first two and next two) as guides
        if (pts.Count >= 4)
        {
            sb.AppendLine(DrawLine(pts[0], pts[1], "#999", 1.0));
            sb.AppendLine(DrawLine(pts[2], pts[3], "#999", 1.0));
        }

        // draw assembled polyline
        sb.Append("  <polyline fill=\"none\" stroke=\"blue\" stroke-width=\"2\" points=\"");
        foreach (PointD p in curvePts)
            sb.Append($"{p.x},{-p.y} ");
        sb.AppendLine("\"/>");

        // draw small circles at key points
        foreach (PointD p in keyPts)
            sb.AppendLine(DrawCircle(p, 3, "red"));

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

}
