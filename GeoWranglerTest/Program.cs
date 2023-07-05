using System;
using Clipper2Lib;
using geoWrangler;

namespace GeoWranglerTest;

internal static class Program
{
    private static void Main()
    {
        grassfire_test();
        unidirectional_bias();
        test_strip_collinear();
        test_fragmentPath();
    }

    private static void grassfire_test()
    {
        PathD original = new()
        {
            new(0, 0),
            new(0, 80),
            new(80, 80),
            new(80, 60),
            new(40, 60),
            new(40, 0),
            new(0, 0)
        };

        PathD median = GeoWrangler.skeleton(original);

        SvgWriter svg = new();
        svg.FillRule = FillRule.EvenOdd;
        
        SvgUtils.AddSubject(svg, original);
        SvgUtils.AddOpenSolution(svg, new PathsD() {median}, true);
        SvgUtils.SaveToFile(svg, "median.svg", FillRule.NonZero, 5500,5500, 10);
    }

    private static void unidirectional_bias()
    {
        SvgWriter svg = new();
        svg.FillRule = FillRule.EvenOdd;
        
        Path64 test = Clipper.MakePath(new[] { 0, 0, 4500, 0, 4500, 4500, 0, 4500 });

        Path64 testvector = Clipper.MakePath(new[] { 0, 0, 500, 0 });

        Paths64 res = Minkowski.Sum(test, testvector, true);
        
        SvgUtils.AddSubject(svg, test);
        SvgUtils.AddSolution(svg, Clipper.Union(res, new() {test}, FillRule.Positive), true);
        SvgUtils.SaveToFile(svg, "curious.svg", FillRule.NonZero, 5500,5500, 10);

        PathD original = new()
        {
            new(0, 0),
            new (0, 50),
            new (20,50),
            new (20, 80),
            new (0,80),
            new(0, 100),
            new(40, 100),
            new(40, 60),
            new(60, 30),
            new(40, 0),
            new(0, 0)
        };

        PathsD original_ = new();
        original_.Add(original);
        
        PathD vector = new()
        {
            new(0.0, 0.0),
            new PointD(0.0, 10.0)
        };

        PathD vector2 = new()
        {
            new(0.0, 0.0),
            new PointD(0.0, -10.0)
        };
        
        // Need both transforms to get the full result. Not sure if implementation-related or a bug in C2.
        // No big deal either way.

        PathsD result = Minkowski.Sum(original, vector, true);

        PathsD result_i = Minkowski.Sum(original, vector2, true);

        result.AddRange(result_i);
        
        SvgUtils.AddSubject(svg, original_);
        SvgUtils.AddSolution(svg, Clipper.Union(result, original_, FillRule.Positive, 2), true);
        SvgUtils.SaveToFile(svg, "tmp.svg", FillRule.NonZero, 150,150, 10);

        PathsD subject = new () { Clipper.Ellipse(new PointD (50.0,50.0),20,20)};

        PathsD result2 = Minkowski.Sum(subject[0], vector, true);
        
        SvgWriter svg2 = new();
        svg2.FillRule = FillRule.EvenOdd;
        SvgUtils.AddSubject(svg2, subject);
        SvgUtils.AddSolution(svg2, Clipper.Union(result2, subject, FillRule.Positive, 2), true);
        SvgUtils.SaveToFile(svg2, "tmp2.svg", FillRule.NonZero, 150,150, 10);
        
        // System("tmp.svg");
        //
    }

    internal static double CrossProduct(Point64 pt1, Point64 pt2, Point64 pt3)
    {
        // typecast to double to avoid potential int overflow
        return ((double) (pt2.X - pt1.X) * (pt3.Y - pt2.Y) -
                (double) (pt2.Y - pt1.Y) * (pt3.X - pt2.X));
    }

    internal static double DotProduct(Point64 pt1, Point64 pt2, Point64 pt3)
    {
        // typecast to double to avoid potential int overflow
        return ((double) (pt2.X - pt1.X) * (pt3.X - pt2.X) +
                (double) (pt2.Y - pt1.Y) * (pt3.Y - pt2.Y));
    }

    internal static double CrossProduct(PointD vec1, PointD vec2)
    {
        return (vec1.y * vec2.x - vec2.y * vec1.x);
    }

    internal static double DotProduct(PointD vec1, PointD vec2)
    {
        return (vec1.x * vec2.x + vec1.y * vec2.y);
    }

    internal static double CrossProduct(PointD pt1, PointD pt2, PointD pt3)
    {
        // typecast to double to avoid potential int overflow
        return ((pt2.x - pt1.x) * (pt3.y - pt2.y) -
                (pt2.y - pt1.y) * (pt3.x - pt2.x));
    }

    internal static double DotProduct(PointD pt1, PointD pt2, PointD pt3)
    {
        // typecast to double to avoid potential int overflow
        return ((pt2.x - pt1.x) * (pt3.x - pt2.x) +
                (pt2.y - pt1.y) * (pt3.y - pt2.y));
    }

    private static double varOffset(Path64 path, PathD norms, int pt, int prevPt)
    {
        // sin(A) < 0: right turning
        // cos(A) < 0: change in angle is more than 90 degree
        int nextPt = pt + 1;
        if (nextPt == path.Count)
        {
            nextPt = 0;
        }
        double sinA = CrossProduct(norms[nextPt], norms[pt], norms[prevPt]);
        double cosA = DotProduct(norms[nextPt], norms[pt], norms[prevPt]);

        return Math.Pow(norms[pt].y + norms[prevPt].y, 2) * 10;
        return 0.0;
    }
    
    private static void test_fragmentPath()
    {
        PathD original = new()
        {
            new(0, 0),
            new(0, 10),
            new(10, 10),
            new(10, 0),
            new(0, 0)
        };

        Fragmenter f = new();

        PathD fragmented_10 = f.fragmentPath(original, 1.0);
        PathD fragmented_05 = f.fragmentPath(original, 0.5);

        // Inject a point into the path that is off-grid for the fragmentation. This has to be collinear to avoid changing the actual shape.
        original.Insert(1, new(0,5.2));
        
        // The fragmentation below should change the result on the first edge due to the injected point
        PathD fragmented_b_10 = f.fragmentPath(original, 1.0);
        PathD fragmented_b_05 = f.fragmentPath(original, 0.5);

        // The below should give consistent results with the original fragmentation - the injected collinear point should have been stripped.
        PathD refragmented_10 = f.refragmentPath(original, 1.0);
        PathD refragmented_05 = f.refragmentPath(original, 0.5);
    }

    private static void test_strip_collinear()
    {
        PathD source = new () {
            new(0.02985, 0.18999),
            new(0.00864, 0.21120),
            new(0.01217, 0.21474),
            new(0.01571, 0.21827),
            new(0.02278, 0.21120),
            new(0.02631, 0.21474),
            new(0.02985, 0.21827),
            new(0.00156, 0.24656),
            new(-0.00196, 0.24302),
            new(-0.00550, 0.23949),
            new(-0.01257, 0.24656),
            new(-0.04121, 0.21791),
            new(-0.04829, 0.21084),
            new(-0.05500, 0.20413),
            new(-0.06207, 0.21120),
            new(-0.05500, 0.21827),
            new(-0.06207, 0.22534),
            new(-0.06560, 0.22181),
            new(-0.06914, 0.21828),
            new(-0.07621, 0.22535),
            new(-0.07267, 0.22888),
            new(-0.06914, 0.23242),
            new(-0.06560, 0.23595),
            new(-0.06207, 0.23949),
            new(-0.05853, 0.24302),
            new(-0.05500, 0.24656),
            new(-0.04793, 0.23949),
            new(-0.05499, 0.23242),
            new(-0.05146, 0.22888),
            new(-0.04792, 0.22535),
            new(-0.01964, 0.25363),
            new(-0.02671, 0.26070),
            new(-0.02318, 0.26424),
            new(-0.01964, 0.26777),
            new(-0.06277, 0.31091),
            new(-0.05570, 0.31798),
            new(-0.04792, 0.31020),
            new(-0.04085, 0.31727),
            new(-0.04792, 0.32434),
            new(-0.04085, 0.33141),
            new(-0.01964, 0.31020),
            new(-0.02671, 0.30313),
            new(-0.03378, 0.31020),
            new(-0.04085, 0.30313),
            new(-0.01257, 0.27484),
            new(-0.00550, 0.28191),
            new(0.00156, 0.27484),
            new(0.04399, 0.31727),
            new(0.05106, 0.31020),
            new(0.04753, 0.30666),
            new(0.04399, 0.30313),
            new(0.05106, 0.29606),
            new(0.05813, 0.30313),
            new(0.06520, 0.29606),
            new(0.06167, 0.29252),
            new(0.05813, 0.28899),
            new(0.05460, 0.28545),
            new(0.04399, 0.27484),
            new(0.03692, 0.28191),
            new(0.04399, 0.28899),
            new(0.04045, 0.29252),
            new(0.03692, 0.29606),
            new(0.00864, 0.26777),
            new(0.01571, 0.26070),
            new(0.00864, 0.25363),
            new(0.05177, 0.21050),
            new(0.04470, 0.20343),
            new(0.03692, 0.21120),
            new(0.03338, 0.20767),
            new(0.02985, 0.20413),
            new(0.03692, 0.19706),
            new(0.03692, 0.19706)
        };

        PathD cleaned = GeoWrangler.stripCollinear(source, precision:6);
    }
}