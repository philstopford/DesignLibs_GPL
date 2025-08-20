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
    enum SamplingMode
    {
        ByMaxSegmentLength,
        ByMaxAngle
    }

    private enum CornerType
    {
        concave,
        convex,
        shortedge
    }
    
    static void Main()
    {
        // --- Configurable input & parameters ---
        List<PointD> original_path = new List<PointD>
        {
            new PointD(0, 0),
            new PointD(0, 100),
            new PointD(20, 100),
            new PointD(20, 80),
            new PointD(40, 80),
            new PointD(40, 60),
            new PointD(60, 60),
            new PointD(60, 40),
            new PointD(80, 40),
            new PointD(80, 0),
            new PointD(0, 0)
        };

        double concave_radius = 20.0;
        double convex_radius = 50.0;
        double edge_resolution = 0.5;
        double angular_resolution = 1.0;
        double short_edge_length = 20.0;

        // Easing configuration:
        EasingStrategy easingMode = EasingStrategy.SmoothBlend;
        double insetFraction = 0.12; // fraction of diagonal length to inset at each end (0..0.5)
        double minInset = 2.0; // minimum inset in world units
        double maxInset = 40.0; // maximum inset
        int diagStraightSample = 0; // samples for straight center section (0 = just endpoints)

        // --- Processing ---
        int[] corner_types = categorizeCorners(original_path, short_edge_length);

        // Need to now walk the corners.
        // Closed path expected.
        PathsD processed = new PathsD();
        for (int i = 0; i < original_path.Count - 1; i++)
        {
            PathD startLine = new PathD { original_path[i] };
            PathD endLine = new PathD { original_path[i] };

            // First and last points are the same, so need to modify in case of first and last point.
            PointD midPrev = (i == 0)
                ? new PointD((original_path[i].x + original_path[^2].x) * 0.5,
                    (original_path[i].y + original_path[^2].y) * 0.5)
                : new PointD((original_path[i].x + original_path[i - 1].x) * 0.5,
                    (original_path[i].y + original_path[i - 1].y) * 0.5);
            startLine.Add(midPrev);

            PointD midNext = (i == original_path.Count - 2)
                ? new PointD((original_path[i].x + original_path[0].x) * 0.5,
                    (original_path[i].y + original_path[0].y) * 0.5)
                : new PointD((original_path[i].x + original_path[i + 1].x) * 0.5,
                    (original_path[i].y + original_path[i + 1].y) * 0.5);
            endLine.Add(midNext);

            // What's the nature of our corner?
            double radius = convex_radius;
            if (corner_types[i] == (int)CornerType.concave) radius = concave_radius;

            processed.Add(processCorner(startLine, endLine, radius, angular_resolution, edge_resolution));
        }

        // Processed now contains a list of corners that can be assembled.
        // Hopefully.
        // Easing here softens the transition between the bezier corner and any short-edge derived diagonals.
        PathD assembled = AssembleWithEasing(processed, corner_types, easingMode, insetFraction, minInset, maxInset,
            diagStraightSample);

        string svg = BuildDetailedSvg(original_path, assembled);
        File.WriteAllText("assembled_with_easing.svg", svg, Encoding.UTF8);
        Console.WriteLine("Wrote assembled_with_easing.svg");
    }

    // --- categorizeCorners (unchanged logic) ---
    static int[] categorizeCorners(PathD path_, double short_edge_length)
    {
        // Need to do some prep work here.
        // Prepare status list
        int[] status = new int[path_.Count];

        // Remove our last point in case the polygon is explicitly closed - this avoids some trouble.
        PathD path = new PathD(path_);
        bool trimmed = false;
        // Remove our last point in case the polygon is explicitly closed - this avoids some trouble.
        if (path.Count > 1 && path[0].x == path[^1].x && path[0].y == path[^1].y)
        {
            trimmed = true;
            path.RemoveAt(path.Count - 1);
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

            double len1 = Math.Sqrt(vx1 * vx1 + vy1 * vy1);
            double len2 = Math.Sqrt(vx2 * vx2 + vy2 * vy2);

            if (len1 <= short_edge_length && len2 <= short_edge_length)
            {
                status[i] = (int)CornerType.shortedge;
                continue;
            }

            // Z component of 3D cross product
            double crossZ = vx1 * vy2 - vy1 * vx2;

            // For CCW polygon, positive crossZ = convex. For CW, negative = convex.
            bool isVertexConvex = isCCW ? (crossZ > 0) : (crossZ < 0);
            status[i] = (int)(isVertexConvex ? CornerType.convex : CornerType.concave);
        }

        if (trimmed)
        {
            // Close path, first and last are the same
            Array.Resize(ref status, status.Length + 1);
            status[^1] = status[0];
        }

        return status;
    }

    // --- corner processing / sampling (unchanged) ---
    static PathD processCorner(PathD startLine, PathD endLine, double radius, double angular_resolution,
        double edge_resolution, SamplingMode mode = SamplingMode.ByMaxSegmentLength)
    {
        // 1. Define base lines (scaled)
        PointD startLineStart = startLine[0];
        PointD startLineEnd = startLine[1];
        PointD endLineStart = endLine[0];
        PointD endLineEnd = endLine[1];

        // 2. Compute curve start/end
        PointD startLength = Minus(startLineEnd, startLineStart);
        PointD startDir = Normalized(startLength);
        PointD endLength = Minus(endLineEnd, endLineStart);
        PointD endDir = Normalized(endLength);

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

        PointD curveEndPoint = Add(endLineStart, Mult(endDir, end_radius));

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

    enum EasingStrategy
    {
        None,
        CubicBezierOptimized,
        QuinticHermite,
        SmoothBlend,
        CircularArc
    }

    static PathD AssembleWithEasing(
        PathsD processedCorners,
        int[] corner_types,
        EasingStrategy strategy, // new: choose implementation
        double insetFraction,
        double minInset,
        double maxInset,
        int diagStraightSample)
    {
        int n = processedCorners.Count;
        PathD outPts = new PathD();

        bool[] isShort = new bool[n];
        for (int i = 0; i < n; i++)
            isShort[i] = (i < corner_types.Length && corner_types[i] == (int)CornerType.shortedge);

        // --- Local helper math ---
        static PointD Minus(PointD a, PointD b) => new PointD(a.x - b.x, a.y - b.y);
        static PointD Add(PointD a, PointD b) => new PointD(a.x + b.x, a.y + b.y);
        static PointD Mul(PointD a, double s) => new PointD(a.x * s, a.y * s);
        static double Length(PointD a) => Math.Sqrt(a.x * a.x + a.y * a.y);

        static PointD Normalized(PointD a)
        {
            double L = Length(a);
            return L < 1e-12 ? new PointD(0, 0) : new PointD(a.x / L, a.y / L);
        }

        static double Dot(PointD a, PointD b) => a.x * b.x + a.y * b.y;
        static PointD Neg(PointD a) => new PointD(-a.x, -a.y);
        static bool PointsEqual(PointD a, PointD b) => Math.Abs(a.x - b.x) < 1e-9 && Math.Abs(a.y - b.y) < 1e-9;

        // Sample cubic bezier uniformly in t (inclusive endpoints)
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

        // Quintic hermite sampler (C2): positions p0,p1; first derivatives m0,m1; second derivatives s0,s1
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

        // Quintic smoothstep (C2)
        static double SmoothStep5(double t) =>
            t <= 0 ? 0 : t >= 1 ? 1 : (6 * Math.Pow(t, 5) - 15 * Math.Pow(t, 4) + 10 * Math.Pow(t, 3));

        // Attempt to build a circle that is tangent to direction t0 at p0 and t1 at p1.
        // Returns false if the geometry is degenerate/unbuildable.
        static bool TryBuildTangentArc(PointD p0, PointD t0Dir, PointD p1, PointD t1Dir, out PointD center,
            out double radius, out double startAng, out double sweep)
        {
            center = new PointD(0, 0);
            radius = 0;
            startAng = 0;
            sweep = 0;
            // Solve for center C such that (p0 - C) dot t0Dir = 0 and (p1 - C) dot t1Dir = 0
            // This is two linear equations for C.x and C.y.
            double a11 = -t0Dir.x;
            double a12 = -t0Dir.y;
            double b1 = -(-t0Dir.x * p0.x - t0Dir.y * p0.y);
            double a21 = -t1Dir.x;
            double a22 = -t1Dir.y;
            double b2 = -(-t1Dir.x * p1.x - t1Dir.y * p1.y);
            double det = a11 * a22 - a12 * a21;
            if (Math.Abs(det) < 1e-12) return false;
            double cx = (b1 * a22 - a12 * b2) / det;
            double cy = (a11 * b2 - b1 * a21) / det;
            center = new PointD(cx, cy);
            PointD v0 = new PointD(p0.x - cx, p0.y - cy);
            PointD v1 = new PointD(p1.x - cx, p1.y - cy);
            double r0 = Length(v0);
            double r1 = Length(v1);
            if (r0 < 1e-9 || r1 < 1e-9) return false;
            // radii must be similar
            if (Math.Abs(r0 - r1) > Math.Max(r0, r1) * 1e-3) return false;
            radius = 0.5 * (r0 + r1);
            startAng = Math.Atan2(v0.y, v0.x);
            double endAng = Math.Atan2(v1.y, v1.x);
            // choose sweep so that tangent directions agree: compute sign from cross of radius->point and tangent
            double rawSweep = endAng - startAng;
            // normalize to [-pi, pi]
            while (rawSweep <= -Math.PI) rawSweep += 2 * Math.PI;
            while (rawSweep > Math.PI) rawSweep -= 2 * Math.PI;
            // verify tangency orientation: for a circle, tangent direction at angle a is ( -sin(a), cos(a) ) times radius sign
            // Compute tangent at start: perp(v0) normalized
            PointD tangentAtStart = new PointD(-v0.y, v0.x);
            double ssign = Math.Sign(Dot(tangentAtStart, t0Dir));
            if (ssign == 0) ssign = 1;
            // Ensure sweep sign matches tangent direction; if not flip
            if (ssign < 0 && rawSweep > 0) rawSweep -= 2 * Math.PI;
            if (ssign > 0 && rawSweep < 0) rawSweep += 2 * Math.PI;
            sweep = rawSweep;
            // final sanity: sweep length not huge
            if (Math.Abs(sweep) > Math.PI * 1.5) return false;
            return true;
        }

        // adaptive sample count helper
        int SampleCountForInset(double ins)
        {
            int sc = Math.Max(8, (int)Math.Ceiling(ins * 40)); // baseline: more samples for larger insets
            return sc;
        }

        // build cubic-bezier optimized segment from a->b using incoming/outgoing tangents
        PathD BuildCubicBezierOptimized(PointD a, PointD b, PointD incomingTangent, PointD outgoingDir, double inset)
        {
            double h = inset * 0.6; // base handle length
            PointD tIn = Normalized(incomingTangent);
            PointD dir = Normalized(outgoingDir);
            if (Length(tIn) < 1e-9) tIn = Neg(dir); // fallback
            PointD b1 = Add(a, Mul(tIn, h)); // leave along poly tangent
            PointD b2 = Minus(b, Mul(dir, h * 0.5)); // approach along diagDir
            int samples = SampleCountForInset(inset);
            return SampleCubicBezier(a, b1, b2, b, samples);
        }

        // build quintic-hermite segment from a->b using first derivatives and zero second derivatives
        PathD BuildQuinticHermite(PointD a, PointD b, PointD incomingTangent, PointD outgoingDir, double inset)
        {
            double h = inset * 0.6;
            PointD m0 = Mul(Normalized(incomingTangent), h); // first derivative at a
            PointD m1 = Mul(Normalized(outgoingDir), h * 0.5); // derivative at b along diag
            PointD s0 = new PointD(0, 0);
            PointD s1 = new PointD(0, 0);
            int samples = SampleCountForInset(inset);
            return SampleQuinticHermite(a, b, m0, m1, s0, s1, samples);
        }

        // build smooth blend between a->b: blend poly path (short straight along incomingTangent) and diag path
        PathD BuildSmoothBlend(PointD a, PointD b, PointD incomingTangent, PointD diagDir, double inset)
        {
            int samples = SampleCountForInset(inset);
            PathD outp = new PathD();
            // T(t): short polyline approximated by moving a small fraction along incomingTangent then to b
            PointD approxEnd = Add(a, Mul(Normalized(incomingTangent), inset * 0.6));
            for (int i = 0; i <= samples; i++)
            {
                double t = (double)i / samples;
                double f = SmoothStep5(t);
                // diag position D(t)
                PointD D = Add(a, Mul(Minus(b, a), t));
                // tangent approx path T(t): linear between a and approxEnd for first half, then to b
                PointD T;
                if (t < 0.5)
                {
                    double tt = t / 0.5;
                    T = Add(Mul(a, 1 - tt), Mul(approxEnd, tt));
                }
                else
                {
                    double tt = (t - 0.5) / 0.5;
                    T = Add(Mul(approxEnd, 1 - tt), Mul(b, tt));
                }

                PointD P = Add(Mul(T, 1 - f), Mul(D, f));
                outp.Add(P);
            }

            return outp;
        }

        // build circular arc tangent blend or return null if failed
        PathD BuildCircularArcOrNull(PointD a, PointD b, PointD incomingTangent, PointD diagDir, double inset)
        {
            if (!TryBuildTangentArc(a, incomingTangent, b, Neg(diagDir), out PointD center, out double radius,
                    out double startAng, out double sweep))
            {
                return null;
            }

            int samples = SampleCountForInset(inset);
            PathD outp = new PathD();
            for (int i = 0; i <= samples; i++)
            {
                double t = (double)i / samples;
                double ang = startAng + sweep * t;
                PointD p = new PointD(center.x + radius * Math.Cos(ang), center.y + radius * Math.Sin(ang));
                outp.Add(p);
            }

            return outp;
        }

        // --- Main loop (kept your original structure) ---
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

            PointD startPt = processedCorners[prevIdx].Last();
            PointD endPt = processedCorners[nextIdx].First();

            if (outPts.Count == 0) outPts.Add(startPt);
            else
            {
                PointD last = outPts[^1];
                if (!(Math.Abs(last.x - startPt.x) < 1e-9 && Math.Abs(last.y - startPt.y) < 1e-9))
                    outPts.Add(startPt);
            }

            // diagonal and inset
            PointD diag = Minus(endPt, startPt);
            double diagLen = Length(diag);
            if (diagLen <= 1e-9)
            {
                outPts.Add(endPt);
                continue;
            }

            PointD diagDir = Normalized(diag);
            double inset = Math.Max(minInset, Math.Min(maxInset, diagLen * insetFraction));
            if (inset * 2.0 > diagLen) inset = diagLen * 0.5 * 0.999;

            PointD S = Add(startPt, Mul(diagDir, inset));
            PointD E = Minus(endPt, Mul(diagDir, inset));

            if (strategy == EasingStrategy.None)
            {
                if (diagStraightSample <= 0)
                {
                    outPts.Add(endPt);
                }
                else
                {
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

            PointD prevTangent = EstimateOutgoingTangent(processedCorners[prevIdx]);
            PointD nextTangent = EstimateIncomingTangent(processedCorners[nextIdx]);

            // choose builder based on strategy
            PathD startSeg = null;
            PathD endSeg = null;

            switch (strategy)
            {
                case EasingStrategy.CubicBezierOptimized:
                    startSeg = BuildCubicBezierOptimized(startPt, S, prevTangent, diagDir, inset);
                    endSeg = BuildCubicBezierOptimized(E, endPt, Neg(diagDir), nextTangent, inset);
                    break;

                case EasingStrategy.QuinticHermite:
                    startSeg = BuildQuinticHermite(startPt, S, prevTangent, diagDir, inset);
                    endSeg = BuildQuinticHermite(E, endPt, Neg(diagDir), nextTangent, inset);
                    break;

                case EasingStrategy.SmoothBlend:
                    startSeg = BuildSmoothBlend(startPt, S, prevTangent, diagDir, inset);
                    endSeg = BuildSmoothBlend(E, endPt, Neg(diagDir), nextTangent, inset);
                    break;

                case EasingStrategy.CircularArc:
                    startSeg = BuildCircularArcOrNull(startPt, S, prevTangent, diagDir, inset);
                    endSeg = BuildCircularArcOrNull(E, endPt, Neg(diagDir), nextTangent, inset);
                    // fallback to cubic if either arc fails
                    if (startSeg == null) startSeg = BuildCubicBezierOptimized(startPt, S, prevTangent, diagDir, inset);
                    if (endSeg == null) endSeg = BuildCubicBezierOptimized(E, endPt, Neg(diagDir), nextTangent, inset);
                    break;
            }

            // append startSeg avoiding duplicate start point
            if (startSeg != null)
            {
                if (outPts.Count > 0 && PointsEqual(outPts[^1], startSeg[0])) outPts.AddRange(startSeg.Skip(1));
                else outPts.AddRange(startSeg);
            }
            else
            {
                // fallback: straight to S
                if (!(PointsEqual(outPts[^1], S))) outPts.Add(S);
            }

            // center straight S->E (sampled)
            if (diagStraightSample <= 0)
            {
                outPts.Add(E);
            }
            else
            {
                for (int s = 1; s <= diagStraightSample; s++)
                {
                    double t = (double)s / (diagStraightSample + 1);
                    PointD p = Add(S, Mul(diagDir, (diagLen - 2 * inset) * t));
                    outPts.Add(p);
                }

                outPts.Add(E);
            }

            // append endSeg avoiding duplicate E
            if (endSeg != null)
            {
                if (PointsEqual(outPts[^1], endSeg[0])) outPts.AddRange(endSeg.Skip(1));
                else outPts.AddRange(endSeg);
            }
            else
            {
                if (!(PointsEqual(outPts[^1], endPt))) outPts.Add(endPt);
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
        minX -= mX;
        maxX += mX;
        minY -= mY;
        maxY += mY;

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
