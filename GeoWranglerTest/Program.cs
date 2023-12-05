using System;
using Clipper2Lib;
using geoWrangler;
using NUnit.Framework;

namespace GeoWranglerTest;

internal static class Program
{
    private static string root_loc = "/d/development/DesignLibs_GPL/geowrangler_test/";

    private static void Main()
    {
        grassfire_test();
        unidirectional_bias();
        test_strip_collinear();
        test_fragmentPath();
        proximity();
        proximity2();
        customBoolean();
        customBoolean2();
        customBoolean3();
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
        SvgUtils.SaveToFile(svg, root_loc + "median.svg", FillRule.NonZero, 5500,5500, 10);

        Assert.AreEqual(6, median.Count);
        Assert.AreEqual(77.5, median[0].x);
        Assert.AreEqual(70, median[0].y);
        Assert.AreEqual(75, median[1].x);
        Assert.AreEqual(70, median[1].y);
        Assert.AreEqual(25, median[2].x);
        Assert.AreEqual(70, median[2].y);
        Assert.AreEqual(20, median[3].x);
        Assert.AreEqual(60, median[3].y);
        Assert.AreEqual(20, median[4].x);
        Assert.AreEqual(13.33, median[4].y);
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
        SvgUtils.SaveToFile(svg, root_loc + "unidirectional1.svg", FillRule.NonZero, 5500,5500, 10);
        
        Assert.AreEqual(4500000, Clipper.Area(res));

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
        SvgUtils.SaveToFile(svg, root_loc + "unidirectional2.svg", FillRule.NonZero, 150,150, 10);

        Assert.AreEqual(3167, Clipper.Area(result));

        PathsD subject = new () { Clipper.Ellipse(new PointD (50.0,50.0),20,20)};

        PathsD result2 = Minkowski.Sum(subject[0], vector, true);
        
        SvgWriter svg2 = new();
        svg2.FillRule = FillRule.EvenOdd;
        SvgUtils.AddSubject(svg2, subject);
        SvgUtils.AddSolution(svg2, Clipper.Union(result2, subject, FillRule.Positive, 2), true);
        SvgUtils.SaveToFile(svg2, root_loc + "unidirectional3.svg", FillRule.NonZero, 150,150, 10);
        
        Assert.LessOrEqual(Clipper.Area(result2) - 785.5610, 0.001);
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

        Assert.AreEqual(81, fragmented_05.Count);
        Assert.AreEqual(41, fragmented_10.Count);

        // Inject a point into the path that is off-grid for the fragmentation. This has to be collinear to avoid changing the actual shape.
        original.Insert(1, new(0,5.2));
        
        // The fragmentation below should change the result on the first edge due to the injected point
        PathD fragmented_b_10 = f.fragmentPath(original, 1.0);
        PathD fragmented_b_05 = f.fragmentPath(original, 0.5);

        Assert.AreEqual(80, fragmented_b_05.Count);
        Assert.AreEqual(40, fragmented_b_10.Count);

        // The below should give consistent results with the original fragmentation - the injected collinear point should have been stripped.
        PathD refragmented_10 = f.refragmentPath(original, 1.0);
        PathD refragmented_05 = f.refragmentPath(original, 0.5);
        
        Assert.AreEqual(81, refragmented_05.Count);
        Assert.AreEqual(41, refragmented_10.Count);
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
        Assert.AreEqual(71, cleaned.Count);
    }

    private static void proximity()
    {
        PathsD input = new ();
        input.Add(Clipper.MakePath(new [] {
        -5, -15,
        -4.97, -13.950000000000001,
        -4.89, -12.91,
        -4.75, -11.870000000000001,
        -4.5600000000000005, -10.84,
        -4.32, -9.82,
        -4.0200000000000005, -8.82,
        -3.67, -7.83,
        -3.27, -6.87,
        -2.82, -5.92,
        -2.32, -5,
        -1.77, -4.11,
        -1.18, -3.24,
        -0.54, -2.41,
        0.14, -1.62,
        0.86, -0.86,
        1.62, -0.14,
        2.41, 0.54,
        3.24, 1.18,
        4.11, 1.77,
        5, 2.32,
        5.92, 2.82,
        6.87, 3.27,
        7.83, 3.67,
        8.82, 4.0200000000000005,
        9.82, 4.32,
        10.84, 4.5600000000000005,
        11.870000000000001, 4.75,
        12.91, 4.89,
        13.950000000000001, 4.97,
        15, 5,
        16.05, 4.97,
        17.09, 4.89,
        18.13, 4.75,
        19.16, 4.5600000000000005,
        20.18, 4.32,
        21.18, 4.0200000000000005,
        22.17, 3.67,
        23.13, 3.27,
        24.080000000000002, 2.82,
        25, 2.32,
        25.89, 1.77,
        26.76, 1.18,
        27.59, 0.54,
        28.38, -0.14,
        29.14, -0.86,
        29.86, -1.62,
        30.54, -2.41,
        31.18, -3.24,
        31.77, -4.11,
        32.32, -5,
        32.82, -5.92,
        33.27, -6.87,
        33.67, -7.83,
        34.02, -8.82,
        34.32, -9.82,
        34.56, -10.84,
        34.75, -11.870000000000001,
        34.89, -12.91,
        34.97, -13.950000000000001,
        35, -15,
        34.99, -16.05,
        34.96, -17.09,
        34.910000000000004, -18.13,
        34.84, -19.16,
        34.74, -20.18,
        34.63, -21.18,
        34.45, -22.490000000000002,
        34.24, -23.77,
        34, -25,
        33.72, -26.18,
        33.410000000000004, -27.310000000000002,
        33.07, -28.38,
        32.71, -29.39,
        32.32, -30.32,
        31.8, -31.38,
        31.25, -32.32,
        30.55, -33.27,
        29.82, -34.02,
        28.93, -34.63,
        27.89, -34.97,
        27.5, -35,
        26.46, -34.93,
        25.43, -34.71,
        24.45, -34.35,
        23.53, -33.86,
        22.68, -33.25,
        21.93, -32.52,
        21.28, -31.69,
        20.76, -30.79,
        20.37, -29.82,
        20.11, -28.8,
        20, -27.76,
        20, -27.5,
        20.05, -26.46,
        20.19, -25.43,
        20.43, -24.45,
        20.81, -23.42,
        21.28, -22.48,
        21.92, -21.59,
        22.650000000000002, -20.88,
        23.54, -20.330000000000002,
        24.560000000000002, -20.03,
        25, -20,
        26.04, -19.89,
        27.03, -19.57,
        27.94, -19.05,
        28.72, -18.35,
        29.330000000000002, -17.5,
        29.76, -16.55,
        29.97, -15.52,
        30, -15,
        29.89, -13.96,
        29.57, -12.97,
        29.05, -12.06,
        28.35, -11.28,
        27.5, -10.67,
        26.55, -10.24,
        25.52, -10.03,
        25, -10,
        23.96, -9.89,
        22.97, -9.57,
        22.06, -9.05,
        21.28, -8.35,
        20.67, -7.5,
        20.240000000000002, -6.55,
        20.03, -5.5200000000000005,
        20, -5,
        19.89, -3.96,
        19.57, -2.97,
        19.05, -2.06,
        18.35, -1.28,
        17.5, -0.67,
        16.55, -0.24,
        15.52, -0.03,
        15, 0,
        13.96, -0.11,
        12.97, -0.43,
        12.06, -0.9500000000000001,
        11.28, -1.6500000000000001,
        10.67, -2.5,
        10.24, -3.45,
        10.03, -4.48,
        10, -5,
        9.89, -6.04,
        9.57, -7.03,
        9.05, -7.94,
        8.35, -8.72,
        7.5, -9.33,
        6.55, -9.76,
        5.5200000000000005, -9.97,
        5, -10,
        3.96, -10.11,
        2.97, -10.43,
        2.06, -10.950000000000001,
        1.28, -11.65,
        0.67, -12.5,
        0.24, -13.450000000000001,
        0.03, -14.48,
        0, -15,
        0.11, -16.04,
        0.43, -17.03,
        0.9500000000000001, -17.94,
        1.6500000000000001, -18.72,
        2.5, -19.330000000000002,
        3.45, -19.76,
        4.48, -19.97,
        5, -20,
        6.04, -20.16,
        6.95, -20.6,
        7.8, -21.28,
        8.47, -22.1,
        8.99, -22.990000000000002,
        9.41, -23.98,
        9.73, -25.060000000000002,
        9.91, -26.07,
        9.99, -27.11,
        10, -27.5,
        9.93, -28.54,
        9.71, -29.57,
        9.35, -30.55,
        8.86, -31.470000000000002,
        8.25, -32.32,
        7.5200000000000005, -33.07,
        6.69, -33.72,
        5.79, -34.24,
        4.82, -34.63,
        3.8000000000000003, -34.89,
        2.7600000000000002, -35,
        2.5, -35,
        1.46, -34.81,
        0.56, -34.32,
        -0.31, -33.54,
        -1.02, -32.660000000000004,
        -1.58, -31.77,
        -2.12, -30.76,
        -2.61, -29.63,
        -2.99, -28.64,
        -3.33, -27.59,
        -3.64, -26.47,
        -3.93, -25.3,
        -4.18, -24.080000000000002,
        -4.4, -22.81,
        -4.59, -21.51,
        -4.71, -20.51,
        -4.8100000000000005, -19.5,
        -4.89, -18.47,
        -4.94, -17.44,
        -4.98, -16.4,
        -5, -15.35,
        -5, -15,
        }));

        GeometryResult res = Proximity.proximityBias(input, new() { false }, 5, 60, 64, 0, 1, 1.03m, 1, false, 1000);
        
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSubject(svgSrc, input);
        SvgUtils.AddSolution(svgSrc, res.geometry, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "proximity.svg", FillRule.NonZero, 800, 800, 10);
        
        Assert.AreEqual(1, res.geometry.Count);
        double area = Clipper.Area(res.geometry[0]);
        Assert.AreEqual(211, res.geometry[0].Count);
        Assert.LessOrEqual(area - 1539.759, 0.001);
    }
    
    private static void proximity2()
    {
        PathsD input = new ();
        input.Add(Clipper.MakePath(new double[] {
            0, -25,
            0, -24,
            0, -23,
            0, -22,
            0, -21,
            0, -20,
            0, -19,
            0, -18,
            0, -17,
            0, -16,
            0, -15,
            0, -14,
            0, -13,
            0, -12,
            0, -11,
            0, -10,
            0, -9,
            0, -8,
            0, -7,
            0, -6,
            0, -5,
            0.11, -3.96,
            0.43, -2.97,
            0.9500000000000001, -2.06,
            1.6500000000000001, -1.28,
            2.5, -0.67,
            3.45, -0.24,
            4.48, -0.03,
            5, 0,
            6.04, -0.11,
            7.03, -0.43,
            7.94, -0.9500000000000001,
            8.72, -1.6500000000000001,
            9.33, -2.5,
            9.76, -3.45,
            9.97, -4.48,
            10, -5,
            10, -6,
            10, -7,
            10, -8,
            10, -9,
            10, -10,
            10, -11,
            10, -12,
            10, -13,
            10, -14,
            10, -15,
            10, -16,
            10, -17,
            10, -18,
            10, -19,
            10, -20,
            10, -21,
            10, -22,
            10, -23,
            10, -24,
            10, -25,
            9.89, -26.04,
            9.57, -27.03,
            9.05, -27.94,
            8.35, -28.72,
            7.5, -29.330000000000002,
            6.55, -29.76,
            5.5200000000000005, -29.97,
            5, -30,
            3.96, -29.89,
            2.97, -29.57,
            2.06, -29.05,
            1.28, -28.35,
            0.67, -27.5,
            0.24, -26.55,
            0.03, -25.52,
            0, -25,
        }));
        input.Add(Clipper.MakePath(new double[] {
            20, -25,
            20, -24,
            20, -23,
            20, -22,
            20, -21,
            20, -20,
            20, -19,
            20, -18,
            20, -17,
            20, -16,
            20, -15,
            20, -14,
            20, -13,
            20, -12,
            20, -11,
            20, -10,
            20, -9,
            20, -8,
            20, -7,
            20, -6,
            20, -5,
            20.11, -3.96,
            20.43, -2.97,
            20.95, -2.06,
            21.650000000000002, -1.28,
            22.5, -0.67,
            23.45, -0.24,
            24.48, -0.03,
            25, 0,
            26.04, -0.11,
            27.03, -0.43,
            27.94, -0.9500000000000001,
            28.72, -1.6500000000000001,
            29.330000000000002, -2.5,
            29.76, -3.45,
            29.97, -4.48,
            30, -5,
            30, -6,
            30, -7,
            30, -8,
            30, -9,
            30, -10,
            30, -11,
            30, -12,
            30, -13,
            30, -14,
            30, -15,
            30, -16,
            30, -17,
            30, -18,
            30, -19,
            30, -20,
            30, -21,
            30, -22,
            30, -23,
            30, -24,
            30, -25,
            29.89, -26.04,
            29.57, -27.03,
            29.05, -27.94,
            28.35, -28.72,
            27.5, -29.330000000000002,
            26.55, -29.76,
            25.52, -29.97,
            25, -30,
            23.96, -29.89,
            22.97, -29.57,
            22.06, -29.05,
            21.28, -28.35,
            20.67, -27.5,
            20.240000000000002, -26.55,
            20.03, -25.52,
            20, -25,
        }));

        GeometryResult res = Proximity.proximityBias(input, new() { false, false }, 5, 30, 64, 0, 1, 1.03m, 1, false, 1000);
        
        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSubject(svgSrc, input);
        SvgUtils.AddSolution(svgSrc, res.geometry, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "proximity2.svg", FillRule.NonZero, 800, 800, 10);
        
        Assert.AreEqual(2, res.geometry.Count);
        double area = Clipper.Area(res.geometry);
        Assert.AreEqual(73, res.geometry[0].Count);
        Assert.AreEqual(73, res.geometry[1].Count);
        Assert.LessOrEqual(area - 1179.626, 0.001);
    }

    // This tests the complex boolean engine in GeoWrangler, which should deliver a keyholed representation of the 
    // boolean output, if needed.
    // This is a little more sophisticated than just calling Clipper directly - it pulls in various GeoWrangler 
    // systems to work its magic, so this is something of an integrated test.
    // Note that some area is lost due to the keyholed representation.
    private static void customBoolean()
    {
        PathsD layerAPaths = new();
        layerAPaths.Add(Clipper.MakePath(new double []
        {
            0,0,
            0,80,
            80,80,
            80,0
        }));
        layerAPaths.Add(Clipper.MakePath(new double []
        {
            90,0,
            90,80,
            100,80,
            100,0
        }));
        double aArea = Clipper.Area(layerAPaths);

        PathsD layerBPaths = new();
        layerBPaths.Add(Clipper.MakePath(new double[]
        {
            10,10,
            10,70,
            20,70,
            20, 10
        }));
        layerBPaths.Add(Clipper.MakePath(new double[]
        {
            60, 10,
            60, 70,
            70, 70,
            70, 10
        }));
        double bArea = Clipper.Area(layerBPaths);
        
        // Subtract layerBPaths from layerAPaths
        PathsD booleanPaths = GeoWrangler.customBoolean(
            firstLayerOperator: (int)GeoWrangler.LayerFlag.none,
            firstLayer: layerAPaths, 
            secondLayerOperator: (int)GeoWrangler.LayerFlag.NOT, 
            secondLayer: layerBPaths, 
            booleanFlag: (int)GeoWrangler.booleanOperation.AND,
            resolution: 1.0,
            extension: 1.03
        );

        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSubject(svgSrc, layerAPaths);
        SvgUtils.AddClip(svgSrc, layerBPaths);
        SvgUtils.AddSolution(svgSrc, booleanPaths, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "customboolean_not.svg", FillRule.NonZero, 800, 800, 10);
        
        double notArea_expected = aArea - bArea;
        double notArea = Clipper.Area(booleanPaths);
        Assert.AreEqual(2, booleanPaths.Count);
        Assert.LessOrEqual(notArea_expected - notArea, 1);
    }
    
    private static void customBoolean2()
    {
        PathsD layerAPaths = new();
        layerAPaths.Add(Clipper.MakePath(new double []
        {
            0,0,
            0,80,
            80,80,
            80,0
        }));
        layerAPaths.Add(Clipper.MakePath(new double []
        {
            -10,30,
            -10,50,
            90,50,
            90,30
        }));

        PathsD layerBPaths = new();
        layerBPaths.Add(Clipper.MakePath(new double[]
        {
            10,10,
            10,70,
            20,70,
            20, 10
        }));
        layerBPaths.Add(Clipper.MakePath(new double[]
        {
            60, 10,
            60, 70,
            70, 70,
            70, 10
        }));
        
        // Subtract layerBPaths from layerAPaths
        PathsD booleanPaths = GeoWrangler.customBoolean(
            firstLayerOperator: (int)GeoWrangler.LayerFlag.none,
            firstLayer: layerAPaths, 
            secondLayerOperator: (int)GeoWrangler.LayerFlag.NOT, 
            secondLayer: layerBPaths, 
            booleanFlag: (int)GeoWrangler.booleanOperation.AND,
            resolution: 1.0,
            extension: 1.03
        );

        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSubject(svgSrc, layerAPaths);
        SvgUtils.AddClip(svgSrc, layerBPaths);
        SvgUtils.AddSolution(svgSrc, booleanPaths, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "customboolean_not2.svg", FillRule.NonZero, 800, 800, 10);
        
        double notArea_expected = -5598;
        double notArea = Clipper.Area(booleanPaths);
        Assert.AreEqual(1, booleanPaths.Count);
        Assert.LessOrEqual(notArea_expected - notArea, 0.001);
    }
    
    private static void customBoolean3()
    {
        PathsD layerAPaths = new();
        layerAPaths.Add(Clipper.MakePath(new double []
        {
            0,0,
            0,80,
            80,80,
            80,0
        }));
        layerAPaths.Add(Clipper.MakePath(new double []
        {
            90,0,
            90,80,
            100,80,
            100,0
        }));
        double aArea = Clipper.Area(layerAPaths);

        PathsD layerBPaths = new();
        layerBPaths.Add(Clipper.MakePath(new double[]
        {
            10,10,
            10,70,
            20,70,
            20, 10
        }));
        layerBPaths.Add(Clipper.MakePath(new double[]
        {
            60, 10,
            60, 70,
            70, 70,
            70, 10
        }));
        double bArea = Clipper.Area(layerBPaths);
        
        // Subtract layerBPaths from layerAPaths
        PathsD tmpPaths = GeoWrangler.customBoolean(
            firstLayerOperator: (int)GeoWrangler.LayerFlag.none,
            firstLayer: layerAPaths, 
            secondLayerOperator: (int)GeoWrangler.LayerFlag.NOT, 
            secondLayer: layerBPaths, 
            booleanFlag: (int)GeoWrangler.booleanOperation.AND,
            resolution: 1.0,
            extension: 1.03
        );

        PathsD cPaths = new();
        cPaths.Add(Clipper.MakePath(new double []
        {
            -10,30,
            -10,35,
            90,35,
            90,30
        }));

        PathsD booleanPaths = GeoWrangler.customBoolean(
            firstLayerOperator: (int)GeoWrangler.LayerFlag.none,
            firstLayer: tmpPaths, 
            secondLayerOperator: (int)GeoWrangler.LayerFlag.none, 
            secondLayer: cPaths, 
            booleanFlag: (int)GeoWrangler.booleanOperation.OR,
            resolution: 1.0,
            extension: 1.03
        );

        SvgWriter svgSrc = new SvgWriter();
        SvgUtils.AddSubject(svgSrc, layerAPaths);
        SvgUtils.AddClip(svgSrc, layerBPaths);
        SvgUtils.AddSolution(svgSrc, booleanPaths, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "customboolean_not3.svg", FillRule.NonZero, 800, 800, 10);
        
        double notArea_expected = -6796;
        double notArea = Clipper.Area(booleanPaths);
        Assert.AreEqual(1, booleanPaths.Count);
        Assert.LessOrEqual(notArea_expected - notArea, 0.001);
    }
}