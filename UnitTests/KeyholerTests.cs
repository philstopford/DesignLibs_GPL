using Clipper2Lib;
using geoWrangler;
using TestHelpers;
using utility;

namespace UnitTests;

public class KeyholerTests
{
    /*
     *
     * There are multiple tests here.
     * singleTest : creates a single keyholed polygon from input. Uses that to create a single slivered polygon. Then run gap removal and sliver removal.
     *              The result should be gR_a having two polygons. sR_a should have one polygon.
     *
     * multiTest : creates a single keyholed polygon with two keyholes. Uses this to also create a single polygon with two slivers.
     *              Result should be that gR_a has three polygons (one outer, two inner). sR_a should have one polygon.
     *
     * multiCutTest : creates a keyholed polygon using a '-' shape, and then cuts across that with a '|'. This tests that the keyholer only generates one keyhole
     *                  Results are that kH has one polygon (the '-' result) and that kH2 should also have one polygon (the '+' from the additional '|' cut).
     *                  gR should have two polygons (each of 5 points, for the outer and the '-' cutter).
     *                  gR2 should have two polygons. The second polygon should have 13 points as it is the '+' from the combined '-' and '|' shapes.
     *                  sL should have one polygon of eight points.
     *                  sL2 should have one polygon of sixteen points.
     *                  sR should have one polygon of five points.
     *                  sR2 should have one polygon of thirteen points.
     *
     * selfOverlapTest : Manually created an overlap within the same Path. Tests how the overlap region is handled for different approaches, to also ensure a keyhole is defined
     *                  Results (x means the correct keyhole is inserted) :
     *
                                Normal		Reversed
                        uR
                        uRP 	x
                        uRNZ 	x			x
                        simp
                        iR
                        iRP
                        iRNZ
                        iR2
                        iR2P	x
                        iR2NZ	x			x
     *
     * comboTest : creates a polygon with a keyhole and a sliver. The sliver is manually defined at this time.
     *              Result should be that gR_a has two polygons. sR_a should have one polygon.
     *
     * simple_islandTest : two not-overlapping instances of singleTest.
     *                      Result should be gR_a having four polygons. sR_a should have two polygons.
     *
     * complex_islandTest : not-overlapping instances of singleTest, multiTest and comboTest.
     *                      Result should be gR_a having seven polygons. sR_a should have three polygons.
     *
     * NOTES:
     *          The keyholer requires orientation rules to be followed to avoid trouble.
     *          _f variants are fragmented forms due to clean-up required when removing vertices. Closing shape and using GeoWrangler might be an option here.
     *          A basic attempt is made to check the orientation of the results; might not work well for complex cases.
     *
     * TODO:
     */

    static string root_loc = TestPath.Get("keyhole_out") + Path.DirectorySeparatorChar;

    // [SetUp]
    public static void KeyholerSetup()
    {
        singleTest();
        spiralTest();
        multiTest();
        multiCutTest();
        selfOverlapTest();
        selfOverlapTest_reversed();
        comboTest();
        simple_islandTest();
        complex_islandTest();
        multiHoleTest();
    }

    [Test]
    public static void singleTest()
    {
        PathD outer = Clipper.MakePath(new double[]
        {
            -200, -200,
            200, -200,
            200, 200,
            -200, 200,
            -200, -200
        });

        PathD inner1 = Clipper.MakePath(new double[]
        {
            -100, -100,
            -100, 100,
            100, 100,
            100, -100,
            -100, -100
        });

        PathsD kHSource = new()
        {
            outer,
            inner1
        };

        // Generate keyholed geometry
        PathsD kH = GeoWrangler.makeKeyHole(kHSource, false, true);
        SvgWriter svgSrc;
        svgSrc = new SvgWriter();
        SvgUtils.AddClip(svgSrc, kHSource);
        SvgUtils.AddSolution(svgSrc, kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "singletest_kh.svg", FillRule.NonZero, 800, 800, 10);

        Assert.That(kH.Count, Is.EqualTo(1));
        // Expected area should be 120000
        double expectedArea = Math.Abs(Clipper.Area(outer)) - Math.Abs(Clipper.Area(inner1));
        // There is a small delta due to the keyhole, but it should be negligible.
        Assert.That(Math.Abs(Clipper.Area(kH)), Is.EqualTo(expectedArea).Within(10.01));
        // Let's make sure our keyhole didn't move too much as well, that could be annoying.
        Assert.That(kH[0][2].x, Is.EqualTo(-100.05).Within(0.001));
        Assert.That(kH[0][3].x, Is.EqualTo(-100.05).Within(0.001));
        Assert.That(kH[0][8].x, Is.EqualTo(-99.95).Within(0.001));
        Assert.That(kH[0][9].x, Is.EqualTo(-99.95).Within(0.001));

        // Generate sliver geometry.
        PathsD sL = new();
        ClipperD c = new(Constants.roundingDecimalPrecision);
        c.AddSubject(outer);
        c.AddClip(kH);
        c.Execute(ClipType.Difference, FillRule.EvenOdd, sL);
        Assert.That(sL.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(sL), Is.EqualTo(40010.0025).Within(0.0001));
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, kH);
        SvgUtils.AddSolution(svgSrc, sL, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "singletest_sl.svg", FillRule.NonZero, 800, 800, 10);

        // Gap removal test
        PathsD gR = GeoWrangler.gapRemoval(kH, 100);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, kH);
        SvgUtils.AddSolution(svgSrc, gR, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "keyholetest_singletest_gr.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(gR.Count, Is.EqualTo(2));
        Assert.That(Math.Abs(Clipper.Area(outer)), Is.EqualTo(Math.Abs(Clipper.Area(gR[0]))));
        Assert.That(Math.Abs(Clipper.Area(inner1)), Is.EqualTo(Math.Abs(Clipper.Area(gR[1]))));
        Assert.That(Clipper.IsPositive(gR[0]), Is.False);
        Assert.That(Clipper.IsPositive(gR[1]), Is.True);

        // Sliver removal test
        PathsD sR = GeoWrangler.gapRemoval(sL, -100);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, sL);
        SvgUtils.AddSolution(svgSrc, sR, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "singletest_sr.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(Clipper.Area(sR), Is.EqualTo(40000));
    }

    [Test]
    public static void spiralTest()
    {
        PathD spiral = Clipper.MakePath(new double[]
        {
            -0.06300, -0.06300,
            -0.06300, 0.06300,
            0.00200, 0.06300,
            0.00200, 0.05700,
            -0.05700, 0.05700,
            -0.05700, -0.05700,
            0.05700, -0.05700,
            0.05700, 0.03700,
            -0.03700, 0.03700,
            -0.03700, -0.03700,
            0.03700, -0.03700,
            0.03700, 0.01700,
            -0.01700, 0.01700,
            -0.01700, -0.01700,
            0.01700, -0.01700,
            0.01700, -0.00300,
            -0.00200, -0.00300,
            -0.00200, 0.00300,
            0.02300, 0.00300,
            0.02300, -0.02300,
            -0.02300, -0.02300,
            -0.02300, 0.02300,
            0.04300, 0.02300,
            0.04300, -0.04300,
            -0.04300, -0.04300,
            -0.04300, 0.04300,
            0.06300, 0.04300,
            0.06300, -0.06300
        });
        spiral = GeoWrangler.resize(spiral, 1000);
        PathsD kH = GeoWrangler.makeKeyHole(spiral, false, false);
        PathD test = GeoWrangler.stripTerminators(kH[0], false);
        Assert.That(Utils.GetSHA256Hash(test), Is.EqualTo(Utils.GetSHA256Hash(spiral)));
        spiral = GeoWrangler.resize(spiral, 0.001);
        kH = GeoWrangler.resize(kH, 0.001);
        SvgWriter svgSrc;
        svgSrc = new SvgWriter();
        SvgUtils.AddClip(svgSrc, spiral);
        SvgUtils.AddSolution(svgSrc, kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "spiral.svg", FillRule.NonZero, 800, 800, 10);
    }

    [Test]
    public static void multiTest()
    {
        PathD outer = Clipper.MakePath(new double[]
        {
            -300, -200,
            300, -200,
            300, 200,
            -300, 200,
            -300, -200
        });

        PathD inner1 = Clipper.MakePath(new double[]
        {
            -200, -100,
            -200, 100,
            -100, 100,
            -100, -100,
            -200, -100
        });

        PathD inner2 = Clipper.MakePath(new double[]
        {
            100, -100,
            100, 100,
            200, 100,
            200, -100,
            100, -100
        });

        PathsD kHSource = new()
        {
            outer,
            inner1,
            inner2
        };


        // Generate keyholed geometry
        PathsD kH = GeoWrangler.makeKeyHole(new PathsD(kHSource), false, true);
        SvgWriter svgSrc;
        svgSrc = new SvgWriter();
        SvgUtils.AddClip(svgSrc, kHSource);
        SvgUtils.AddSolution(svgSrc, kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "multitest_kh.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(kH.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(kH), Is.EqualTo(-199979.995d).Within(0.001));

        // Generate sliver geometry.
        PathsD sL = new();
        ClipperD c = new(Constants.roundingDecimalPrecision);
        c.AddSubject(outer);
        c.AddClip(kH);
        c.Execute(ClipType.Difference, FillRule.EvenOdd, sL);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, kH);
        SvgUtils.AddSolution(svgSrc, sL, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "multitest_sl.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(sL.Count, Is.EqualTo(2));
        Assert.That(Clipper.Area(sL), Is.EqualTo(40020.005).Within(0.001));

        // Gap removal test
        PathsD gR = GeoWrangler.gapRemoval(kH, 100);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, gR, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "multitest_gr.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(gR.Count, Is.EqualTo(3));
        Assert.That(Clipper.Area(gR), Is.EqualTo(-200000));

        // Sliver removal test
        PathsD sR = GeoWrangler.gapRemoval(sL, -100);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, sL);
        SvgUtils.AddSolution(svgSrc, sR, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "multitest_sr.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(sR.Count, Is.EqualTo(2));
        Assert.That(Clipper.Area(sR), Is.EqualTo(40000));
    }

    [Test]
    public static void multiCutTest()
    {
        PathD outer = new()
        {
            new PointD(0, 0),
            new PointD(400, 0),
            new PointD(400, 400),
            new PointD(0, 400),
            new PointD(0, 0)
        };

        PathD inner1 = new()
        {
            new PointD(50, 150),
            new PointD(50, 250),
            new PointD(350, 250),
            new PointD(350, 150),
            new PointD(50, 150)
        };

        PathD inner2 = new()
        {
            new PointD(150, 50),
            new PointD(150, 350),
            new PointD(250, 350),
            new PointD(250, 50),
            new PointD(150, 50)
        };

        PathsD kHSource = new()
        {
            outer,
            inner1
        };

        // Generate keyholed geometry
        PathsD kH = GeoWrangler.makeKeyHole(new PathsD(kHSource), false, true);
        SvgWriter svgSrc;
        svgSrc = new SvgWriter();
        SvgUtils.AddClip(svgSrc, kHSource);
        SvgUtils.AddSolution(svgSrc, kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "multicuttest_kh.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(kH.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(kH), Is.EqualTo(-129994.9975d).Within(0.001));

        PathsD kHSource2 = new();
        kHSource2.AddRange(kH);
        kHSource2.Add(inner2);

        // Generate keyholed geometry
        PathsD kH2 = GeoWrangler.makeKeyHole(new PathsD(kHSource2), false, true);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, kH2, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "multicuttest_kh2.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(kH2.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(kH2), Is.EqualTo(-139989.995).Within(0.001));

        // Generate sliver geometry.
        ClipperD c = new(Constants.roundingDecimalPrecision);
        PathsD sL = new();
        c.AddSubject(outer);
        c.AddClip(kH);
        c.Execute(ClipType.Difference, FillRule.EvenOdd, sL);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, sL, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "multicuttest_sl.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(sL.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(sL), Is.EqualTo(30005.0025).Within(0.001));

        c.Clear();
        PathsD sL2 = new();
        c.AddSubject(outer);
        c.AddClip(kH2);
        c.Execute(ClipType.Difference, FillRule.EvenOdd, sL2);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, sL2, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "multicuttest_sl2.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(sL2.Count, Is.EqualTo(2));
        Assert.That(Clipper.Area(sL2), Is.EqualTo(20010.005).Within(0.001));

        // Gap removal test
        PathsD gR = GeoWrangler.gapRemoval(kH, 100);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, gR, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "multicuttest_gr.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(gR.Count, Is.EqualTo(2));
        Assert.That(Clipper.Area(gR), Is.EqualTo(-130000).Within(0.001));

        PathsD gR2 = GeoWrangler.gapRemoval(kH2, 100);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, gR2, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "multicuttest_gr2.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(gR2.Count, Is.EqualTo(3));
        Assert.That(Clipper.Area(gR2), Is.EqualTo(-140000).Within(0.001));

        PathsD kHSource3 = new();
        kHSource3.AddRange(gR);

        PathsD kHSource4 = new();
        kHSource4.AddRange(gR2);

        // Generate keyholed geometry
        PathsD kH3 = GeoWrangler.makeKeyHole(new PathsD(kHSource3), false, true);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, kH3, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "multicuttest_kh3.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(kH3.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(kH3), Is.EqualTo(-129994.9975).Within(0.001));

        PathsD kH4 = GeoWrangler.makeKeyHole(new PathsD(kHSource4), false, true);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, kH4, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "multicuttest_kh4.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(kH4.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(kH4), Is.EqualTo(-139989.995).Within(0.001));

        // Sliver removal test
        PathsD sR = GeoWrangler.gapRemoval(sL, -100);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, sR, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "multicuttest_sr.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(sR.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(sR), Is.EqualTo(30000).Within(0.001));

        PathsD sR2 = GeoWrangler.gapRemoval(sL2, -100);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, sR2, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "multicuttest_sr2.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(sR2.Count, Is.EqualTo(2));
        Assert.That(Clipper.Area(sR2), Is.EqualTo(20000).Within(0.001));
    }

    [Test]
    public static void selfOverlapTest()
    {
        PathD outer = new()
        {
            new PointD(0, 0),
            new PointD(110, 0),
            new PointD(110, 50),
            new PointD(50, 50),
            new PointD(50, 150),
            new PointD(150, 150),
            new PointD(150, 50),
            new PointD(90, 50),
            new PointD(90, 0),
            new PointD(200, 0),
            new PointD(200, 200),
            new PointD(0, 200),
            new PointD(0, 0)
        };

        // decomposer test
        PathsD dSource = new() { outer };
        PathsD decomp = GeoWrangler.decompose(dSource);
        SvgWriter svgSrc;
        svgSrc = new SvgWriter();
        SvgUtils.AddClip(svgSrc, dSource);
        SvgUtils.AddSolution(svgSrc, decomp, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(decomp.Count, Is.EqualTo(2));
        Assert.That(Clipper.Area(decomp), Is.EqualTo(30000).Within(0.001));

        PathsD kHD = GeoWrangler.makeKeyHole(dSource, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, kHD, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_khd.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(kHD.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(kHD), Is.EqualTo(-29994.9975).Within(0.001));

        // keyholer test
        PathsD kHSource = new() { outer };
        PathsD kH = GeoWrangler.makeKeyHole(kHSource, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, kHSource);
        SvgUtils.AddSolution(svgSrc, kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_kh.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(kH.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(kH), Is.EqualTo(-29994.9975).Within(0.001));

        ClipperD c = new(Constants.roundingDecimalPrecision);
        c.AddSubject(outer);

        PathsD unionRes = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, unionRes);

        PathsD unionRes_kH = GeoWrangler.makeKeyHole(unionRes, false, true);

        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, unionRes);
        SvgUtils.AddSolution(svgSrc, unionRes_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_unionres.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(unionRes_kH.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(unionRes_kH), Is.EqualTo(29000).Within(0.001));

        PathsD unionResc = GeoWrangler.close(unionRes);
        PathsD unionResc_kH = GeoWrangler.makeKeyHole(unionResc, false, true);

        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, unionResc);
        SvgUtils.AddSolution(svgSrc, unionResc_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_unionresc.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(unionResc_kH.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(unionResc_kH), Is.EqualTo(29000).Within(0.001));

        PathsD unionResP = new();
        c.Execute(ClipType.Union, FillRule.Positive, unionResP);

        PathsD unionResP_kH = GeoWrangler.makeKeyHole(unionResP, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, unionResP);
        SvgUtils.AddSolution(svgSrc, unionResP_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_unionresp.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(unionResP_kH.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(unionResP_kH), Is.EqualTo(-29994.9975).Within(0.001));

        PathsD unionResPc = GeoWrangler.close(unionResP);
        PathsD unionResPc_kH = GeoWrangler.makeKeyHole(unionResPc, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, unionResPc);
        SvgUtils.AddSolution(svgSrc, unionResPc_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_unionrespc.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(unionResPc_kH.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(unionResPc_kH), Is.EqualTo(-29994.9975).Within(0.001));

        // seems good - get keyhole
        PathsD unionResNZ = new();
        c.Execute(ClipType.Union, FillRule.NonZero, unionResNZ);

        PathsD unionResNZ_kH = GeoWrangler.makeKeyHole(unionResNZ, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, unionResNZ);
        SvgUtils.AddSolution(svgSrc, unionResNZ_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_unionresnz.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(unionResNZ_kH.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(unionResNZ_kH), Is.EqualTo(-29994.9975).Within(0.001));

        PathsD unionResNZc = GeoWrangler.close(unionResNZ);
        PathsD unionResNZc_kH = GeoWrangler.makeKeyHole(unionResNZc, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, unionResNZc);
        SvgUtils.AddSolution(svgSrc, unionResNZc_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_unionresnzc.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(unionResNZc_kH.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(unionResNZc_kH), Is.EqualTo(-29994.9975).Within(0.001));

        PathsD simplifyRes = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, simplifyRes);
        simplifyRes = GeoWrangler.stripCollinear(simplifyRes);

        PathsD simplifyRes_kH = GeoWrangler.makeKeyHole(simplifyRes, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, simplifyRes);
        SvgUtils.AddSolution(svgSrc, simplifyRes_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_simplifyres.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(simplifyRes_kH.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(simplifyRes_kH), Is.EqualTo(29000).Within(0.001));

        PathsD simplifyResc = GeoWrangler.close(simplifyRes);
        PathsD simplifyResc_kH = GeoWrangler.makeKeyHole(simplifyResc, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, simplifyResc);
        SvgUtils.AddSolution(svgSrc, simplifyResc_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_simplifyresc.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(simplifyResc_kH.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(simplifyResc_kH), Is.EqualTo(29000).Within(0.001));

        PathsD simplifyRes2 = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, simplifyRes2);
        PathsD simplifyRes2_kH = GeoWrangler.makeKeyHole(simplifyRes2, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, simplifyRes2);
        SvgUtils.AddSolution(svgSrc, simplifyRes2_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_simplifyres2.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(simplifyRes2_kH.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(simplifyRes2_kH), Is.EqualTo(29000).Within(0.001));

        PathsD simplifyRes2c = GeoWrangler.close(simplifyRes2);
        PathsD simplifyRes2c_kH = GeoWrangler.makeKeyHole(simplifyRes2c, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, simplifyRes2c);
        SvgUtils.AddSolution(svgSrc, simplifyRes2c_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_simplifyres2c.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(simplifyRes2c_kH.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(simplifyRes2c_kH), Is.EqualTo(29000).Within(0.001));

        // no good - no result
        PathsD intRes = new();
        c.Execute(ClipType.Intersection, FillRule.EvenOdd, intRes);
        Assert.That(intRes.Count, Is.EqualTo(0));

        PathsD intRes_kH = GeoWrangler.makeKeyHole(intRes, false, true);
        Assert.That(intRes_kH.Count, Is.EqualTo(0));

        PathsD intResc = GeoWrangler.close(intRes);
        PathsD intResc_kH = GeoWrangler.makeKeyHole(intResc, false, true);
        Assert.That(intResc_kH.Count, Is.EqualTo(0));

        // no good - no result
        PathsD intResP = new();
        c.Execute(ClipType.Intersection, FillRule.Positive, intResP);
        PathsD intResP_kH = GeoWrangler.makeKeyHole(intResP, false, true);
        Assert.That(intResP_kH.Count, Is.EqualTo(0));
        PathsD intResPc = GeoWrangler.close(intResP);
        PathsD intResPc_kH = GeoWrangler.makeKeyHole(intResPc, false, true);
        Assert.That(intResPc_kH.Count, Is.EqualTo(0));

        // no good - no result
        PathsD intResNZ = new();
        c.Execute(ClipType.Intersection, FillRule.NonZero, intResNZ);
        PathsD intResNZ_kH = GeoWrangler.makeKeyHole(intResNZ, false, true);
        Assert.That(intResNZ_kH.Count, Is.EqualTo(0));
        PathsD intResNZc = GeoWrangler.close(intResNZ);
        PathsD intResNZc_kH = GeoWrangler.makeKeyHole(intResNZc, false, true);
        Assert.That(intResNZc_kH.Count, Is.EqualTo(0));

        RectD bounds = Clipper.GetBounds(new PathsD { outer });
        PathD bb = new()
        {
            new PointD(bounds.left, bounds.bottom),
            new PointD(bounds.left, bounds.top),
            new PointD(bounds.right, bounds.top),
            new PointD(bounds.right, bounds.bottom)
        };

        c.Clear();
        c.AddSubject(bb);
        c.AddClip(outer);

        PathsD intRes2 = new();
        c.Execute(ClipType.Intersection, FillRule.EvenOdd, intRes2);

        PathsD intRes2_kH = GeoWrangler.makeKeyHole(intRes2, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, intRes2);
        SvgUtils.AddSolution(svgSrc, intRes2_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_intres2.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(intRes2_kH.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(intRes2_kH), Is.EqualTo(29000).Within(0.001));

        PathsD intRes2c = GeoWrangler.close(intRes2);
        PathsD intRes2c_kH = GeoWrangler.makeKeyHole(intRes2c, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, intRes2c);
        SvgUtils.AddSolution(svgSrc, intRes2c_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_intres2c.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(intRes2c_kH.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(intRes2c_kH), Is.EqualTo(29000).Within(0.001));

        PathsD intRes2P = new();
        c.Execute(ClipType.Intersection, FillRule.Positive, intRes2P);
        Assert.That(intRes2P.Count, Is.EqualTo(2));

        PathsD intRes2P_kH = GeoWrangler.makeKeyHole(intRes2P, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, intRes2P);
        SvgUtils.AddSolution(svgSrc, intRes2P_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_intres2p.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(intRes2P_kH.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(intRes2P_kH), Is.EqualTo(-29994.9975).Within(0.001));

        PathsD intRes2Pc = GeoWrangler.close(intRes2P);
        PathsD intRes2Pc_kH = GeoWrangler.makeKeyHole(intRes2Pc, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, intRes2Pc);
        SvgUtils.AddSolution(svgSrc, intRes2Pc_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_intres2pc.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(intRes2Pc_kH.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(intRes2Pc_kH), Is.EqualTo(-29994.9975).Within(0.001));

        // seems good - get keyhole
        PathsD intRes2NZ = new();
        c.Execute(ClipType.Intersection, FillRule.NonZero, intRes2NZ);

        PathsD intRes2NZ_kH = GeoWrangler.makeKeyHole(intRes2NZ, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, intRes2NZ);
        SvgUtils.AddSolution(svgSrc, intRes2NZ_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_intres2nz.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(intRes2NZ_kH.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(intRes2NZ_kH), Is.EqualTo(-29994.9975).Within(0.001));

        PathsD intRes2NZc = GeoWrangler.close(intRes2NZ);
        PathsD intRes2NZc_kH = GeoWrangler.makeKeyHole(intRes2NZc, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, intRes2NZc);
        SvgUtils.AddSolution(svgSrc, intRes2NZc_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_intres2nzc.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(intRes2NZc_kH.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(intRes2NZc_kH), Is.EqualTo(-29994.9975).Within(0.001));
    }

    [Test]
    public static void selfOverlapTest_reversed()
    {
        PathD outer = Clipper.MakePath(new double[]
        {
            0, 0,
            110, 0,
            110, 50,
            50, 50,
            50, 150,
            150, 150,
            150, 50,
            90, 50,
            90, 0,
            200, 0,
            200, 200,
            0, 200,
            0, 0
        });

        outer.Reverse();

        PathsD dSource = new() { outer };
        PathsD decomp = GeoWrangler.decompose(dSource);

        SvgWriter svgSrc;
        svgSrc = new SvgWriter();
        SvgUtils.AddClip(svgSrc, dSource);
        SvgUtils.AddSolution(svgSrc, decomp, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_reversed.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(decomp.Count, Is.EqualTo(2));
        // Delta expected. Sign for dSource is opposite to decomp.
        Assert.That(Clipper.Area(decomp), Is.EqualTo(30000));
        Assert.That(Clipper.Area(dSource), Is.EqualTo(-31000));

        PathsD kHD = GeoWrangler.makeKeyHole(dSource, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, dSource);
        SvgUtils.AddSolution(svgSrc, kHD, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_reversed_khd.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(kHD.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(kHD), Is.EqualTo(-29994.9975).Within(0.001));

        // keyholer test
        PathsD kHSource = new() { outer };
        PathsD kH = GeoWrangler.makeKeyHole(kHSource, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, kHSource);
        SvgUtils.AddSolution(svgSrc, kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_reversed_kh.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(kH.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(kH), Is.EqualTo(-29994.9975).Within(0.001));

        ClipperD c = new(Constants.roundingDecimalPrecision);
        c.AddSubject(outer);

        PathsD unionRes = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, unionRes);

        PathsD unionRes_kH = GeoWrangler.makeKeyHole(unionRes, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, unionRes);
        SvgUtils.AddSolution(svgSrc, unionRes_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_reversed_khunion.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(unionRes_kH.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(unionRes_kH), Is.EqualTo(Clipper.Area(unionRes)));

        PathsD unionResc = GeoWrangler.close(unionRes);
        PathsD unionResc_kH = GeoWrangler.makeKeyHole(unionResc, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, unionResc);
        SvgUtils.AddSolution(svgSrc, unionResc_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_reversed_khunionc.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(unionResc_kH.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(unionResc_kH), Is.EqualTo(Clipper.Area(unionResc)));

        PathsD unionResP = new();
        c.Clear();
        c.AddSubject(unionResc);
        c.Execute(ClipType.Union, FillRule.Positive, unionResP);

        PathsD unionResP_kH = GeoWrangler.makeKeyHole(unionResP, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, unionResP);
        SvgUtils.AddSolution(svgSrc, unionResP_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_reversed_khunionp.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(unionResP_kH.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(unionResP_kH), Is.EqualTo(Clipper.Area(unionResP)));

        PathsD unionResPc = GeoWrangler.close(unionResP);
        PathsD unionResPc_kH = GeoWrangler.makeKeyHole(unionResPc, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, unionResPc);
        SvgUtils.AddSolution(svgSrc, unionResPc_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_reversed_khunionpc.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(unionResPc_kH.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(unionResPc_kH), Is.EqualTo(Clipper.Area(unionResPc)));

        PathsD unionResNZ = new();
        c.Execute(ClipType.Union, FillRule.NonZero, unionResNZ);

        PathsD unionResNZ_kH = GeoWrangler.makeKeyHole(unionResNZ, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, unionResNZ);
        SvgUtils.AddSolution(svgSrc, unionResNZ_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_reversed_khunionnz.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(unionResNZ_kH.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(unionResNZ_kH), Is.EqualTo(Clipper.Area(unionResNZ)));

        PathsD unionResNZc = GeoWrangler.close(unionResNZ);
        PathsD unionResNZc_kH = GeoWrangler.makeKeyHole(unionResNZc, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, unionResNZc);
        SvgUtils.AddSolution(svgSrc, unionResNZc_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_reversed_khunionnzc.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(unionResNZc_kH.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(unionResNZc_kH), Is.EqualTo(Clipper.Area(unionResNZc)));

        // no good - overlap region is a gap.
        PathsD simplifyRes = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, simplifyRes);
        simplifyRes = GeoWrangler.stripCollinear(simplifyRes);
        PathsD simplifyRes_kH = GeoWrangler.makeKeyHole(simplifyRes, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, simplifyRes);
        SvgUtils.AddSolution(svgSrc, simplifyRes_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_reversed_khsimplify.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(simplifyRes_kH.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(simplifyRes_kH), Is.EqualTo(Clipper.Area(simplifyRes)));

        PathsD simplifyResc = GeoWrangler.close(simplifyRes);
        PathsD simplifyResc_kH = GeoWrangler.makeKeyHole(simplifyResc, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, simplifyResc);
        SvgUtils.AddSolution(svgSrc, simplifyRes_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_reversed_khsimplifyc.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(simplifyResc_kH.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(simplifyResc_kH), Is.EqualTo(Clipper.Area(simplifyResc)));

        PathsD simplifyRes2 = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, simplifyRes2);
        PathsD simplifyRes2_kH = GeoWrangler.makeKeyHole(simplifyRes2, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, simplifyRes2);
        SvgUtils.AddSolution(svgSrc, simplifyRes2_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_reversed_khsimplify2.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(simplifyRes2_kH.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(simplifyRes2_kH), Is.EqualTo(Clipper.Area(simplifyRes2)));

        PathsD simplifyRes2c = GeoWrangler.close(simplifyRes2);
        PathsD simplifyRes2c_kH = GeoWrangler.makeKeyHole(simplifyRes2c, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, simplifyRes2c);
        SvgUtils.AddSolution(svgSrc, simplifyRes2c_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_reversed_khsimplify2c.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(simplifyRes2c_kH.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(simplifyRes2c_kH), Is.EqualTo(Clipper.Area(simplifyRes2c)));

        // no good - no result
        PathsD intRes = new();
        c.Execute(ClipType.Intersection, FillRule.EvenOdd, intRes);
        Assert.That(intRes.Count, Is.EqualTo(0));
        PathsD intRes_kH = GeoWrangler.makeKeyHole(intRes, false, true);
        Assert.That(intRes_kH.Count, Is.EqualTo(0));

        PathsD intResc = GeoWrangler.close(intRes);
        PathsD intResc_kH = GeoWrangler.makeKeyHole(intResc, false, true);
        Assert.That(intResc_kH.Count, Is.EqualTo(0));

        // no good - no result
        PathsD intResP = new();
        c.Execute(ClipType.Intersection, FillRule.Positive, intResP);
        PathsD intResP_kH = GeoWrangler.makeKeyHole(intResP, false, true);
        Assert.That(intResP_kH.Count, Is.EqualTo(0));
        PathsD intResPc = GeoWrangler.close(intResP);
        PathsD intResPc_kH = GeoWrangler.makeKeyHole(intResPc, false, true);
        Assert.That(0, Is.EqualTo(0).Within(intResPc_kH.Count));

        // no good - no result
        PathsD intResNZ = new();
        c.Execute(ClipType.Intersection, FillRule.NonZero, intResNZ);
        PathsD intResNZ_kH = GeoWrangler.makeKeyHole(intResNZ, false, true);
        Assert.That(intResNZ_kH.Count, Is.EqualTo(0));
        PathsD intResNZc = GeoWrangler.close(intResNZ);
        PathsD intResNZc_kH = GeoWrangler.makeKeyHole(intResNZc, false, true);
        Assert.That(intResNZc_kH.Count, Is.EqualTo(0));

        RectD bounds = Clipper.GetBounds(new PathsD { outer });
        PathD bb = new()
        {
            new PointD(bounds.left, bounds.bottom),
            new PointD(bounds.left, bounds.top),
            new PointD(bounds.right, bounds.top),
            new PointD(bounds.right, bounds.bottom)
        };

        c.Clear();
        c.AddSubject(bb);
        c.AddClip(outer);

        PathsD intRes2 = new();
        c.Execute(ClipType.Intersection, FillRule.EvenOdd, intRes2);

        PathsD intRes2_kH = GeoWrangler.makeKeyHole(intRes2, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, intRes2);
        SvgUtils.AddSolution(svgSrc, intRes2_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_reversed_khint2.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(intRes2_kH.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(intRes2_kH), Is.EqualTo(Clipper.Area(intRes2)));

        PathsD intRes2c = GeoWrangler.close(intRes2);
        PathsD intRes2c_kH = GeoWrangler.makeKeyHole(intRes2c, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, intRes2c);
        SvgUtils.AddSolution(svgSrc, intRes2c_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_reversed_khint2c.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(intRes2c_kH.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(intRes2c_kH), Is.EqualTo(Clipper.Area(intRes2c)));

        // no good - no result
        PathsD intRes2P = new();
        c.Execute(ClipType.Intersection, FillRule.Positive, intRes2P);
        Assert.That(intRes2P.Count, Is.EqualTo(0));

        // No results - no geometry
        PathsD intRes2P_kH = GeoWrangler.makeKeyHole(intRes2P, false, true);
        Assert.That(intRes2P_kH.Count, Is.EqualTo(0));

        PathsD intRes2Pc = GeoWrangler.close(intRes2P);
        PathsD intRes2Pc_kH = GeoWrangler.makeKeyHole(intRes2Pc, false, true);
        Assert.That(intRes2Pc_kH.Count, Is.EqualTo(0));

        // seems good - get keyhole
        PathsD intRes2NZ = new();
        c.Execute(ClipType.Intersection, FillRule.NonZero, intRes2NZ);
        Assert.That(intRes2NZ.Count, Is.EqualTo(2));

        PathsD intRes2NZ_kH = GeoWrangler.makeKeyHole(intRes2NZ, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, intRes2NZ);
        SvgUtils.AddSolution(svgSrc, intRes2NZ_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_reversed_khint2nz.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(intRes2NZ_kH.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(intRes2NZ_kH), Is.EqualTo(-29994.9975).Within(0.001));

        PathsD intRes2NZc = GeoWrangler.close(intRes2NZ);
        PathsD intRes2NZc_kH = GeoWrangler.makeKeyHole(intRes2NZc, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, intRes2NZc);
        SvgUtils.AddSolution(svgSrc, intRes2NZc_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_reversed_khint2nzc.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(intRes2NZc_kH.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(intRes2NZc_kH), Is.EqualTo(-29994.9975).Within(0.001));
    }

    [Test]
    public static void comboTest()
    {
        // Manually create the sliver.
        PathD outer = Clipper.MakePath(new double[]
        {
            -200, -200,
            300, -200,
            300, 200,
            -200, 200,
            -200, -99.9,
            -300, -99.9,
            -300, -100.1,
            -200, -100.1,
            -200, -200
        });

        PathD inner1 = Clipper.MakePath(new double[]
        {
            100, -100,
            100, 100,
            200, 100,
            200, -100,
            100, -100
        });

        PathsD kHSource = new()
        {
            outer,
            inner1
        };

        PathsD kH = GeoWrangler.makeKeyHole(kHSource, false, true);
        SvgWriter svgSrc;
        svgSrc = new SvgWriter();
        SvgUtils.AddClip(svgSrc, kHSource);
        SvgUtils.AddSolution(svgSrc, kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "combotest_kh.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(kH.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(kH), Is.EqualTo(-180009.9975).Within(0.001));
        Assert.That(Clipper.Area(kHSource), Is.EqualTo(180020));

        // Gap removal test
        PathsD gR = GeoWrangler.gapRemoval(kH, 100);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, kH);
        SvgUtils.AddSolution(svgSrc, gR, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "combotest_gr.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(gR.Count, Is.EqualTo(2));
        Assert.That(Math.Abs(Clipper.Area(gR)), Is.EqualTo(Math.Abs(Clipper.Area(kHSource))));

        // Sliver removal test
        PathsD sR = GeoWrangler.gapRemoval(kH, -100);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, kH);
        SvgUtils.AddSolution(svgSrc, sR, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "combotest_sr.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(sR.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(sR), Is.EqualTo(-179989.9975).Within(0.001));
    }

    [Test]
    public static void simple_islandTest()
    {
        PathD outer1 = Clipper.MakePath(new double[]
        {
            -200, -200,
            200, -200,
            200, 200,
            -200, 200,
            -200, -200
        });

        PathD inner1 = Clipper.MakePath(new double[]
        {
            -100, -100,
            -100, 100,
            100, 100,
            100, -100,
            -100, -100
        });

        PathD outer2 = Clipper.MakePath(new double[]
        {
            -200, 400,
            200, 400,
            200, 800,
            -200, 800,
            -200, 400
        });

        PathD inner2 = Clipper.MakePath(new double[]
        {
            -100, 500,
            -100, 700,
            100, 700,
            100, 500,
            -100, 500
        });

        PathsD kHSource = new()
        {
            outer1,
            inner1,
            outer2,
            inner2
        };

        // Generate keyholed geometry
        PathsD kH = GeoWrangler.makeKeyHole(kHSource, false, true);
        SvgWriter svgSrc;
        svgSrc = new SvgWriter();
        SvgUtils.AddClip(svgSrc, kHSource);
        SvgUtils.AddSolution(svgSrc, kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "simpleislandtest_kh.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(kH.Count, Is.EqualTo(2));
        Assert.That(Clipper.Area(kH), Is.EqualTo(-239979.995).Within(0.001));
        Assert.That(Clipper.Area(kHSource), Is.EqualTo(240000));

        // Gap removal test
        PathsD gR = GeoWrangler.gapRemoval(kH, 100);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, kH);
        SvgUtils.AddSolution(svgSrc, gR, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "simpleislandtest_gr.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(gR.Count, Is.EqualTo(4));
        Assert.That(Math.Abs(Clipper.Area(gR)), Is.EqualTo(Math.Abs(Clipper.Area(kHSource))));

        // Generate sliver geometry.
        PathsD sL = new();
        ClipperD c = new(Constants.roundingDecimalPrecision);
        c.AddSubject(outer1);
        c.AddSubject(outer2);
        c.AddClip(kH);
        c.Execute(ClipType.Difference, FillRule.EvenOdd, sL);

        // Sliver removal test
        PathsD sR = GeoWrangler.gapRemoval(sL, -100);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, sL);
        SvgUtils.AddSolution(svgSrc, sR, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "simpleislandtest_sr.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(sR.Count, Is.EqualTo(2));
        Assert.That(Clipper.Area(sR), Is.EqualTo(80000));
    }

    [Test]
    public static void complex_islandTest()
    {
        // Island 1 - mix of sliver and gap at the end.
        // Manually create the sliver.
        PathD outer1 = Clipper.MakePath(new double[]
        {
            -200, -200,
            300, -200,
            300, 200,
            -200, 200,
            -200, -99.9,
            -300, -99.9,
            -300, -100.1,
            -200, -100.1,
            -200, -200
        });

        PathD inner1 = Clipper.MakePath(new double[]
        {
            100, -100,
            100, 100,
            200, 100,
            200, -100,
            100, -100
        });

        PathsD kHSource = new()
        {
            outer1,
            inner1
        };

        PathsD kH = GeoWrangler.makeKeyHole(kHSource, false, true);
        SvgWriter svgSrc;
        svgSrc = new SvgWriter();
        SvgUtils.AddClip(svgSrc, kHSource);
        SvgUtils.AddSolution(svgSrc, kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "complexislandtest_kh_1.svg", FillRule.NonZero, 800, 800, 10);
        double area_kh = Clipper.Area(kH);
        Assert.That(kH.Count, Is.EqualTo(1));
        Assert.That(area_kh, Is.EqualTo(-180009.9975).Within(0.001));
        Assert.That(Clipper.Area(kHSource), Is.EqualTo(180020));

        // Island 2 - simple single hole
        PathD outer2 = Clipper.MakePath(new double[]
        {
            -200, 400,
            200, 400,
            200, 800,
            -200, 800,
            -200, 400
        });

        PathD inner2 = Clipper.MakePath(new double[]
        {
            -100, 500,
            -100, 700,
            100, 700,
            100, 500,
            -100, 500
        });

        kHSource.Clear();
        kHSource.Add(outer2);
        kHSource.Add(inner2);

        PathsD kH2 = GeoWrangler.makeKeyHole(kHSource, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, kHSource);
        SvgUtils.AddSolution(svgSrc, kH2, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "complexislandtest_kh_2.svg", FillRule.NonZero, 800, 800, 10);
        double area_kh2 = Clipper.Area(kH2);
        Assert.That(kH2.Count, Is.EqualTo(1));
        Assert.That(area_kh2, Is.EqualTo(-119989.9975).Within(0.001));
        Assert.That(Clipper.Area(kHSource), Is.EqualTo(120000));

        // Island 3 - dual hole hole
        PathD outer3 = Clipper.MakePath(new double[]
        {
            -300, 1200,
            300, 1200,
            300, 1600,
            -300, 1600,
            -300, 1200
        });

        PathD inner3_1 = Clipper.MakePath(new double[]
        {
            -200, 1300,
            -200, 1500,
            -100, 1500,
            -100, 1300,
            -200, 1300
        });

        PathD inner3_2 = Clipper.MakePath(new double[]
        {
            100, 1300,
            100, 1500,
            200, 1500,
            200, 1300,
            100, 1300
        });

        kHSource.Clear();
        kHSource.Add(outer3);
        kHSource.Add(inner3_1);
        kHSource.Add(inner3_2);

        PathsD kH3 = GeoWrangler.makeKeyHole(kHSource, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, kHSource);
        SvgUtils.AddSolution(svgSrc, kH3, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "complexislandtest_kh_3.svg", FillRule.NonZero, 800, 800, 10);
        double area_kh3 = Clipper.Area(kH3);
        Assert.That(kH3.Count, Is.EqualTo(1));
        Assert.That(area_kh3, Is.EqualTo(-199979.995).Within(0.001));
        Assert.That(Clipper.Area(kHSource), Is.EqualTo(200000));

        kHSource.Clear();
        kHSource.Add(outer1);
        kHSource.Add(inner1);

        kHSource.Add(outer2);
        kHSource.Add(inner2);

        PathsD kH12 = GeoWrangler.makeKeyHole(kHSource, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, kHSource);
        SvgUtils.AddSolution(svgSrc, kH12, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "complexislandtest_kh_12.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(kH12.Count, Is.EqualTo(2));
        Assert.That(Clipper.Area(kH12), Is.EqualTo(area_kh + area_kh2));
        Assert.That(Clipper.Area(kHSource), Is.EqualTo(300020));

        kHSource.Clear();
        kHSource.Add(outer2);
        kHSource.Add(inner2);

        kHSource.Add(outer3);
        kHSource.Add(inner3_1);

        PathsD kH23h1 = GeoWrangler.makeKeyHole(kHSource, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, kHSource);
        SvgUtils.AddSolution(svgSrc, kH23h1, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "complexislandtest_kh_23h1.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(kH23h1.Count, Is.EqualTo(2));
        double area_ = Clipper.Area(kH23h1);
        Assert.That(Clipper.Area(kH23h1), Is.EqualTo(-339979.995).Within(0.001));
        Assert.That(Clipper.Area(kHSource), Is.EqualTo(340000));

        kHSource.Clear();
        kHSource.Add(outer2);
        kHSource.Add(inner2);

        kHSource.Add(outer3);
        kHSource.Add(inner3_2);

        PathsD kH23h2 = GeoWrangler.makeKeyHole(kHSource, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, kHSource);
        SvgUtils.AddSolution(svgSrc, kH23h2, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "complexislandtest_kh_23h2.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(kH23h2.Count, Is.EqualTo(2));
        Assert.That(Clipper.Area(kH23h2), Is.EqualTo(-339979.995).Within(0.001));
        Assert.That(Clipper.Area(kHSource), Is.EqualTo(340000));

        kHSource.Clear();
        kHSource.Add(outer1);
        kHSource.Add(inner1);

        kHSource.Add(outer3);
        kHSource.Add(inner3_1);
        kHSource.Add(inner3_2);

        PathsD kH13 = GeoWrangler.makeKeyHole(kHSource, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, kHSource);
        SvgUtils.AddSolution(svgSrc, kH13, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "complexislandtest_kh_13.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(kH13.Count, Is.EqualTo(2));
        Assert.That(Clipper.Area(kH13), Is.EqualTo(area_kh + area_kh3));
        Assert.That(Clipper.Area(kHSource), Is.EqualTo(380020));

        kHSource.Clear();
        kHSource.Add(outer1);
        kHSource.Add(inner1);

        kHSource.Add(outer2);
        kHSource.Add(inner2);

        kHSource.Add(outer3);
        kHSource.Add(inner3_1);
        kHSource.Add(inner3_2);

        PathsD kH123 = GeoWrangler.makeKeyHole(kHSource, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, kHSource);
        SvgUtils.AddSolution(svgSrc, kH123, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "complexislandtest_kh_123.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(kH123.Count, Is.EqualTo(3));
        Assert.That(Math.Abs(area_kh) + Math.Abs(area_kh2) + Math.Abs(area_kh3), Is.EqualTo(Math.Abs(Clipper.Area(kH123))));
        Assert.That(Clipper.Area(kHSource), Is.EqualTo(500020));

        // Gap removal test
        PathsD gR = GeoWrangler.gapRemoval(kH123, 100);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, kH123);
        SvgUtils.AddSolution(svgSrc, gR, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "complexislandtest_gr.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(gR.Count, Is.EqualTo(7));
        Assert.That(-Clipper.Area(gR), Is.EqualTo(Clipper.Area(kHSource)));

        // Sliver removal test
        PathsD sR = GeoWrangler.gapRemoval(kH123, -100);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, kH123);
        SvgUtils.AddSolution(svgSrc, sR, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "complexislandtest_sr.svg", FillRule.NonZero, 800, 800, 10);
        Assert.That(sR.Count, Is.EqualTo(3));
        Assert.That(Clipper.Area(sR), Is.EqualTo(-499959.989).Within(0.001));
    }

    [Test]
    public static void multiHoleTest()
    {
        PathsD kHSource = new();
        PathD path_1 = Clipper.MakePath(new double[]
        {
            50, 0,
            0, 0,
            0, 155,
            50, 155,
            50, 0
        });
        kHSource.Add(path_1);

        PathD path_2 = Clipper.MakePath(new double[]
        {
            35, 135,
            35, 150,
            5, 150,
            5, 135,
            35, 135
        });
        kHSource.Add(path_2);

        PathD path_3 = Clipper.MakePath(new double[]
        {
            22, 95,
            22, 125,
            5, 125,
            5, 95,
            22, 95
        });
        kHSource.Add(path_3);

        PathD path_4 = Clipper.MakePath(new double[]
        {
            35, 45,
            35, 75,
            5, 75,
            5, 45,
            35, 45
        });
        kHSource.Add(path_4);

        PathD path_5 = Clipper.MakePath(new double[]
        {
            35, 5,
            35, 35,
            5, 35,
            5, 5,
            35, 5
        });
        kHSource.Add(path_5);

        // Generate keyholed geometry
        PathsD kH = GeoWrangler.makeKeyHole(kHSource, false, true);
        SvgWriter svgSrc;
        svgSrc = new SvgWriter();
        SvgUtils.AddClip(svgSrc, kHSource);
        SvgUtils.AddSolution(svgSrc, kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "multiholetest_1.svg", FillRule.NonZero, 800, 800, 10);
        double area_kh = Clipper.Area(kH);
        Assert.That(kH.Count, Is.EqualTo(1));
        Assert.That(area_kh, Is.EqualTo(-4987.99).Within(0.001));
        Assert.That(Clipper.Area(kHSource), Is.EqualTo(-4990));
    }
}