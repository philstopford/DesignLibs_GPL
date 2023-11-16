using System;
using Clipper2Lib;
using geoWrangler;
using NUnit.Framework;
using utility;

namespace KeyHolerTest;

internal class Program
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

    private static string root_loc = "/d/development/DesignLibs_GPL/keyhole_out/";
    
    private static void Main(string[] args)
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

    private static void singleTest()
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

        Assert.AreEqual(1, kH.Count);
        // Expected area should be 120000
        double expectedArea = Math.Abs(Clipper.Area(outer)) - Math.Abs(Clipper.Area(inner1));
        // There is a small delta due to the keyhole, but it should be negligible.
        Assert.AreEqual(expectedArea,Math.Abs(Clipper.Area(kH)), 10.01);
        // Let's make sure our keyhole didn't move too much as well, that could be annoying.
        Assert.AreEqual(kH[0][2].x, -100.05, 0.001);
        Assert.AreEqual(kH[0][3].x, -100.05, 0.001);
        Assert.AreEqual(kH[0][8].x, -99.95, 0.001);
        Assert.AreEqual(kH[0][9].x, -99.95, 0.001);

        // Generate sliver geometry.
        PathsD sL = new();
        ClipperD c = new(Constants.roundingDecimalPrecision);
        c.AddSubject(outer);
        c.AddClip(kH);
        c.Execute(ClipType.Difference, FillRule.EvenOdd, sL);
        Assert.AreEqual(1, sL.Count);
        Assert.AreEqual(Clipper.Area(sL), 40010.0025, 0.0001);
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
        Assert.AreEqual(2, gR.Count);
        Assert.AreEqual(Math.Abs(Clipper.Area(gR[0])), Math.Abs(Clipper.Area(outer)));
        Assert.AreEqual(Math.Abs(Clipper.Area(gR[1])), Math.Abs(Clipper.Area(inner1)));
        Assert.False(Clipper.IsPositive(gR[0]));
        Assert.True(Clipper.IsPositive(gR[1]));

        // Sliver removal test
        PathsD sR = GeoWrangler.gapRemoval(sL, -100);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, sL);
        SvgUtils.AddSolution(svgSrc, sR, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "singletest_sr.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(40000, Clipper.Area(sR));
    }

    private static void spiralTest()
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
        Assert.AreEqual(Utils.GetSHA256Hash(spiral), Utils.GetSHA256Hash(test));
        spiral = GeoWrangler.resize(spiral, 0.001);
        kH = GeoWrangler.resize(kH, 0.001);
        SvgWriter svgSrc;
        svgSrc = new SvgWriter();
        SvgUtils.AddClip(svgSrc, spiral);
        SvgUtils.AddSolution(svgSrc, kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "spiral.svg", FillRule.NonZero, 800, 800, 10);
    }

    private static void multiTest()
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
        Assert.AreEqual(1, kH.Count);
        Assert.AreEqual(-199979.995d, Clipper.Area(kH), 0.001);

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
        Assert.AreEqual(2, sL.Count);
        Assert.AreEqual(40020.005, Clipper.Area(sL), 0.001);
        
        // Gap removal test
        PathsD gR = GeoWrangler.gapRemoval(kH, 100);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, gR, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "multitest_gr.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(3, gR.Count);
        Assert.AreEqual(-200000, Clipper.Area(gR));
        
        // Sliver removal test
        PathsD sR = GeoWrangler.gapRemoval(sL, -100);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, sL);
        SvgUtils.AddSolution(svgSrc, sR, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "multitest_sr.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(2, sR.Count);
        Assert.AreEqual(40000, Clipper.Area(sR));
    }

    private static void multiCutTest()
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
        Assert.AreEqual(1, kH.Count);
        Assert.AreEqual(-129994.9975d, Clipper.Area(kH), 0.001);
        
        PathsD kHSource2 = new();
        kHSource2.AddRange(kH);
        kHSource2.Add(inner2);

        // Generate keyholed geometry
        PathsD kH2 = GeoWrangler.makeKeyHole(new PathsD(kHSource2), false, true);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, kH2, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "multicuttest_kh2.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(1, kH2.Count);
        Assert.AreEqual(-139989.995, Clipper.Area(kH2), 0.001);
        
        // Generate sliver geometry.
        ClipperD c = new(Constants.roundingDecimalPrecision);
        PathsD sL = new();
        c.AddSubject(outer);
        c.AddClip(kH);
        c.Execute(ClipType.Difference, FillRule.EvenOdd, sL);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, sL, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "multicuttest_sl.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(1, sL.Count);
        Assert.AreEqual(30005.0025, Clipper.Area(sL), 0.001);

        c.Clear();
        PathsD sL2 = new();
        c.AddSubject(outer);
        c.AddClip(kH2);
        c.Execute(ClipType.Difference, FillRule.EvenOdd, sL2);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, sL2, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "multicuttest_sl2.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(2, sL2.Count);
        Assert.AreEqual(20010.005, Clipper.Area(sL2), 0.001);
        
        // Gap removal test
        PathsD gR = GeoWrangler.gapRemoval(kH, 100);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, gR, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "multicuttest_gr.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(2, gR.Count);
        Assert.AreEqual(-130000, Clipper.Area(gR), 0.001);

        PathsD gR2 = GeoWrangler.gapRemoval(kH2, 100);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, gR2, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "multicuttest_gr2.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(3, gR2.Count);
        Assert.AreEqual(-140000, Clipper.Area(gR2), 0.001);
        
        PathsD kHSource3 = new();
        kHSource3.AddRange(gR);

        PathsD kHSource4 = new();
        kHSource4.AddRange(gR2);

        // Generate keyholed geometry
        PathsD kH3 = GeoWrangler.makeKeyHole(new PathsD(kHSource3), false, true);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, kH3, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "multicuttest_kh3.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(1, kH3.Count);
        Assert.AreEqual(-129994.9975, Clipper.Area(kH3), 0.001);

        PathsD kH4 = GeoWrangler.makeKeyHole(new PathsD(kHSource4), false, true);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, kH4, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "multicuttest_kh4.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(1, kH4.Count);
        Assert.AreEqual(-139989.995, Clipper.Area(kH4), 0.001);

        // Sliver removal test
        PathsD sR = GeoWrangler.gapRemoval(sL, -100);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, sR, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "multicuttest_sr.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(1, sR.Count);
        Assert.AreEqual(30000, Clipper.Area(sR), 0.001);
        
        PathsD sR2 = GeoWrangler.gapRemoval(sL2, -100);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, sR2, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "multicuttest_sr2.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(2, sR2.Count);
        Assert.AreEqual(20000, Clipper.Area(sR2), 0.001);
    }

    private static void selfOverlapTest()
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
        Assert.AreEqual(2, decomp.Count);
        Assert.AreEqual(30000, Clipper.Area(decomp), 0.001);

        PathsD kHD = GeoWrangler.makeKeyHole(dSource, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddSolution(svgSrc, kHD, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_khd.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(1, kHD.Count);
        Assert.AreEqual(-29994.9975, Clipper.Area(kHD), 0.001);

        // keyholer test
        PathsD kHSource = new() { outer };
        PathsD kH = GeoWrangler.makeKeyHole(kHSource, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, kHSource);
        SvgUtils.AddSolution(svgSrc, kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_kh.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(1, kH.Count);
        Assert.AreEqual(-29994.9975, Clipper.Area(kH), 0.001);

        ClipperD c = new(Constants.roundingDecimalPrecision);
        c.AddSubject(outer);

        PathsD unionRes = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, unionRes);

        PathsD unionRes_kH = GeoWrangler.makeKeyHole(unionRes, false, true);

        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, unionRes);
        SvgUtils.AddSolution(svgSrc, unionRes_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_unionres.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(1, unionRes_kH.Count);
        Assert.AreEqual(29000, Clipper.Area(unionRes_kH), 0.001);

        PathsD unionResc = GeoWrangler.close(unionRes);
        PathsD unionResc_kH = GeoWrangler.makeKeyHole(unionResc, false, true);

        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, unionResc);
        SvgUtils.AddSolution(svgSrc, unionResc_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_unionresc.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(1, unionResc_kH.Count);
        Assert.AreEqual(29000, Clipper.Area(unionResc_kH), 0.001);

        PathsD unionResP = new();
        c.Execute(ClipType.Union, FillRule.Positive, unionResP);

        PathsD unionResP_kH = GeoWrangler.makeKeyHole(unionResP, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, unionResP);
        SvgUtils.AddSolution(svgSrc, unionResP_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_unionresp.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(1, unionResP_kH.Count);
        Assert.AreEqual(-29994.9975, Clipper.Area(unionResP_kH), 0.001);

        PathsD unionResPc = GeoWrangler.close(unionResP);
        PathsD unionResPc_kH = GeoWrangler.makeKeyHole(unionResPc, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, unionResPc);
        SvgUtils.AddSolution(svgSrc, unionResPc_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_unionrespc.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(1, unionResPc_kH.Count);
        Assert.AreEqual(-29994.9975, Clipper.Area(unionResPc_kH), 0.001);

        // seems good - get keyhole
        PathsD unionResNZ = new();
        c.Execute(ClipType.Union, FillRule.NonZero, unionResNZ);

        PathsD unionResNZ_kH = GeoWrangler.makeKeyHole(unionResNZ, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, unionResNZ);
        SvgUtils.AddSolution(svgSrc, unionResNZ_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_unionresnz.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(1, unionResNZ_kH.Count);
        Assert.AreEqual(-29994.9975, Clipper.Area(unionResNZ_kH), 0.001);

        PathsD unionResNZc = GeoWrangler.close(unionResNZ);
        PathsD unionResNZc_kH = GeoWrangler.makeKeyHole(unionResNZc, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, unionResNZc);
        SvgUtils.AddSolution(svgSrc, unionResNZc_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_unionresnzc.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(1, unionResNZc_kH.Count);
        Assert.AreEqual(-29994.9975, Clipper.Area(unionResNZc_kH), 0.001);

        PathsD simplifyRes = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, simplifyRes);
        simplifyRes = GeoWrangler.stripCollinear(simplifyRes);

        PathsD simplifyRes_kH = GeoWrangler.makeKeyHole(simplifyRes, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, simplifyRes);
        SvgUtils.AddSolution(svgSrc, simplifyRes_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_simplifyres.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(1, simplifyRes_kH.Count);
        Assert.AreEqual(29000, Clipper.Area(simplifyRes_kH), 0.001);

        PathsD simplifyResc = GeoWrangler.close(simplifyRes);
        PathsD simplifyResc_kH = GeoWrangler.makeKeyHole(simplifyResc, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, simplifyResc);
        SvgUtils.AddSolution(svgSrc, simplifyResc_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_simplifyresc.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(1, simplifyResc_kH.Count);
        Assert.AreEqual(29000, Clipper.Area(simplifyResc_kH), 0.001);

        PathsD simplifyRes2 = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, simplifyRes2);
        PathsD simplifyRes2_kH = GeoWrangler.makeKeyHole(simplifyRes2, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, simplifyRes2);
        SvgUtils.AddSolution(svgSrc, simplifyRes2_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_simplifyres2.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(1, simplifyRes2_kH.Count);
        Assert.AreEqual(29000, Clipper.Area(simplifyRes2_kH), 0.001);

        PathsD simplifyRes2c = GeoWrangler.close(simplifyRes2);
        PathsD simplifyRes2c_kH = GeoWrangler.makeKeyHole(simplifyRes2c, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, simplifyRes2c);
        SvgUtils.AddSolution(svgSrc, simplifyRes2c_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_simplifyres2c.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(1, simplifyRes2c_kH.Count);
        Assert.AreEqual(29000, Clipper.Area(simplifyRes2c_kH), 0.001);

        // no good - no result
        PathsD intRes = new();
        c.Execute(ClipType.Intersection, FillRule.EvenOdd, intRes);
        Assert.AreEqual(0, intRes.Count);

        PathsD intRes_kH = GeoWrangler.makeKeyHole(intRes, false, true);
        Assert.AreEqual(0, intRes_kH.Count);

        PathsD intResc = GeoWrangler.close(intRes);
        PathsD intResc_kH = GeoWrangler.makeKeyHole(intResc, false, true);
        Assert.AreEqual(0, intResc_kH.Count);

        // no good - no result
        PathsD intResP = new();
        c.Execute(ClipType.Intersection, FillRule.Positive, intResP);
        PathsD intResP_kH = GeoWrangler.makeKeyHole(intResP, false, true);
        Assert.AreEqual(0, intResP_kH.Count);
        PathsD intResPc = GeoWrangler.close(intResP);
        PathsD intResPc_kH = GeoWrangler.makeKeyHole(intResPc, false, true);
        Assert.AreEqual(0, intResPc_kH.Count);

        // no good - no result
        PathsD intResNZ = new();
        c.Execute(ClipType.Intersection, FillRule.NonZero, intResNZ);
        PathsD intResNZ_kH = GeoWrangler.makeKeyHole(intResNZ, false, true);
        Assert.AreEqual(0, intResNZ_kH.Count);
        PathsD intResNZc = GeoWrangler.close(intResNZ);
        PathsD intResNZc_kH = GeoWrangler.makeKeyHole(intResNZc, false, true);
        Assert.AreEqual(0, intResNZc_kH.Count);

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
        Assert.AreEqual(1, intRes2_kH.Count);
        Assert.AreEqual(29000, Clipper.Area(intRes2_kH), 0.001);

        PathsD intRes2c = GeoWrangler.close(intRes2);
        PathsD intRes2c_kH = GeoWrangler.makeKeyHole(intRes2c, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, intRes2c);
        SvgUtils.AddSolution(svgSrc, intRes2c_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_intres2c.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(1, intRes2c_kH.Count);
        Assert.AreEqual(29000, Clipper.Area(intRes2c_kH), 0.001);
        
        PathsD intRes2P = new();
        c.Execute(ClipType.Intersection, FillRule.Positive, intRes2P);
        Assert.AreEqual(intRes2P.Count, 2);

        PathsD intRes2P_kH = GeoWrangler.makeKeyHole(intRes2P, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, intRes2P);
        SvgUtils.AddSolution(svgSrc, intRes2P_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_intres2p.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(1, intRes2P_kH.Count);
        Assert.AreEqual(-29994.9975, Clipper.Area(intRes2P_kH), 0.001);

        PathsD intRes2Pc = GeoWrangler.close(intRes2P);
        PathsD intRes2Pc_kH = GeoWrangler.makeKeyHole(intRes2Pc, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, intRes2Pc);
        SvgUtils.AddSolution(svgSrc, intRes2Pc_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_intres2pc.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(1, intRes2Pc_kH.Count);
        Assert.AreEqual(-29994.9975, Clipper.Area(intRes2Pc_kH), 0.001);

        // seems good - get keyhole
        PathsD intRes2NZ = new();
        c.Execute(ClipType.Intersection, FillRule.NonZero, intRes2NZ);
        
        PathsD intRes2NZ_kH = GeoWrangler.makeKeyHole(intRes2NZ, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, intRes2NZ);
        SvgUtils.AddSolution(svgSrc, intRes2NZ_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_intres2nz.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(1, intRes2NZ_kH.Count);
        Assert.AreEqual(-29994.9975, Clipper.Area(intRes2NZ_kH), 0.001);

        PathsD intRes2NZc = GeoWrangler.close(intRes2NZ);
        PathsD intRes2NZc_kH = GeoWrangler.makeKeyHole(intRes2NZc, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, intRes2NZc);
        SvgUtils.AddSolution(svgSrc, intRes2NZc_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_intres2nzc.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(1, intRes2NZc_kH.Count);
        Assert.AreEqual(-29994.9975, Clipper.Area(intRes2NZc_kH), 0.001);
    }

    private static void selfOverlapTest_reversed()
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
        Assert.AreEqual(2, decomp.Count);
        // Delta expected. Sign for dSource is opposite to decomp.
        Assert.AreEqual(30000, Clipper.Area(decomp));
        Assert.AreEqual(-31000, Clipper.Area(dSource));

        PathsD kHD = GeoWrangler.makeKeyHole(dSource, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, dSource);
        SvgUtils.AddSolution(svgSrc, kHD, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_reversed_khd.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(1, kHD.Count);
        Assert.AreEqual(-29994.9975, Clipper.Area(kHD), 0.001);

        // keyholer test
        PathsD kHSource = new() { outer };
        PathsD kH = GeoWrangler.makeKeyHole(kHSource, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, kHSource);
        SvgUtils.AddSolution(svgSrc, kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_reversed_kh.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(1, kH.Count);
        Assert.AreEqual(-29994.9975, Clipper.Area(kH), 0.001);

        ClipperD c = new(Constants.roundingDecimalPrecision);
        c.AddSubject(outer);

        PathsD unionRes = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, unionRes);

        PathsD unionRes_kH = GeoWrangler.makeKeyHole(unionRes, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, unionRes);
        SvgUtils.AddSolution(svgSrc, unionRes_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_reversed_khunion.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(1, unionRes_kH.Count);
        Assert.AreEqual(Clipper.Area(unionRes),Clipper.Area(unionRes_kH));

        PathsD unionResc = GeoWrangler.close(unionRes);
        PathsD unionResc_kH = GeoWrangler.makeKeyHole(unionResc, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, unionResc);
        SvgUtils.AddSolution(svgSrc, unionResc_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_reversed_khunionc.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(1, unionResc_kH.Count);
        Assert.AreEqual(Clipper.Area(unionResc), Clipper.Area(unionResc_kH));

        PathsD unionResP = new();
        c.Clear();
        c.AddSubject(unionResc);
        c.Execute(ClipType.Union, FillRule.Positive, unionResP);

        PathsD unionResP_kH = GeoWrangler.makeKeyHole(unionResP, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, unionResP);
        SvgUtils.AddSolution(svgSrc, unionResP_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_reversed_khunionp.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(1, unionResP_kH.Count);
        Assert.AreEqual(Clipper.Area(unionResP), Clipper.Area(unionResP_kH));

        PathsD unionResPc = GeoWrangler.close(unionResP);
        PathsD unionResPc_kH = GeoWrangler.makeKeyHole(unionResPc, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, unionResPc);
        SvgUtils.AddSolution(svgSrc, unionResPc_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_reversed_khunionpc.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(1, unionResPc_kH.Count);
        Assert.AreEqual(Clipper.Area(unionResPc), Clipper.Area(unionResPc_kH));

        PathsD unionResNZ = new();
        c.Execute(ClipType.Union, FillRule.NonZero, unionResNZ);

        PathsD unionResNZ_kH = GeoWrangler.makeKeyHole(unionResNZ, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, unionResNZ);
        SvgUtils.AddSolution(svgSrc, unionResNZ_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_reversed_khunionnz.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(1, unionResNZ_kH.Count);
        Assert.AreEqual(Clipper.Area(unionResNZ), Clipper.Area(unionResNZ_kH));

        PathsD unionResNZc = GeoWrangler.close(unionResNZ);
        PathsD unionResNZc_kH = GeoWrangler.makeKeyHole(unionResNZc, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, unionResNZc);
        SvgUtils.AddSolution(svgSrc, unionResNZc_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_reversed_khunionnzc.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(1, unionResNZc_kH.Count);
        Assert.AreEqual(Clipper.Area(unionResNZc), Clipper.Area(unionResNZc_kH));

        // no good - overlap region is a gap.
        PathsD simplifyRes = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, simplifyRes);
        simplifyRes = GeoWrangler.stripCollinear(simplifyRes);
        PathsD simplifyRes_kH = GeoWrangler.makeKeyHole(simplifyRes, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, simplifyRes);
        SvgUtils.AddSolution(svgSrc, simplifyRes_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_reversed_khsimplify.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(1, simplifyRes_kH.Count);
        Assert.AreEqual(Clipper.Area(simplifyRes), Clipper.Area(simplifyRes_kH));

        PathsD simplifyResc = GeoWrangler.close(simplifyRes);
        PathsD simplifyResc_kH = GeoWrangler.makeKeyHole(simplifyResc, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, simplifyResc);
        SvgUtils.AddSolution(svgSrc, simplifyRes_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_reversed_khsimplifyc.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(1, simplifyResc_kH.Count);
        Assert.AreEqual(Clipper.Area(simplifyResc), Clipper.Area(simplifyResc_kH));

        PathsD simplifyRes2 = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, simplifyRes2);
        PathsD simplifyRes2_kH = GeoWrangler.makeKeyHole(simplifyRes2, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, simplifyRes2);
        SvgUtils.AddSolution(svgSrc, simplifyRes2_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_reversed_khsimplify2.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(1, simplifyRes2_kH.Count);
        Assert.AreEqual(Clipper.Area(simplifyRes2), Clipper.Area(simplifyRes2_kH));

        PathsD simplifyRes2c = GeoWrangler.close(simplifyRes2);
        PathsD simplifyRes2c_kH = GeoWrangler.makeKeyHole(simplifyRes2c, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, simplifyRes2c);
        SvgUtils.AddSolution(svgSrc, simplifyRes2c_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_reversed_khsimplify2c.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(1, simplifyRes2c_kH.Count);
        Assert.AreEqual(Clipper.Area(simplifyRes2c), Clipper.Area(simplifyRes2c_kH));

        // no good - no result
        PathsD intRes = new();
        c.Execute(ClipType.Intersection, FillRule.EvenOdd, intRes);
        Assert.AreEqual(0, intRes.Count);
        PathsD intRes_kH = GeoWrangler.makeKeyHole(intRes, false, true);
        Assert.AreEqual(0, intRes_kH.Count);

        PathsD intResc = GeoWrangler.close(intRes);
        PathsD intResc_kH = GeoWrangler.makeKeyHole(intResc, false, true);
        Assert.AreEqual(0, intResc_kH.Count);

        // no good - no result
        PathsD intResP = new();
        c.Execute(ClipType.Intersection, FillRule.Positive, intResP);
        PathsD intResP_kH = GeoWrangler.makeKeyHole(intResP, false, true);
        Assert.AreEqual(0, intResP_kH.Count);
        PathsD intResPc = GeoWrangler.close(intResP);
        PathsD intResPc_kH = GeoWrangler.makeKeyHole(intResPc, false, true);
        Assert.AreEqual(0, 0, intResPc_kH.Count);

        // no good - no result
        PathsD intResNZ = new();
        c.Execute(ClipType.Intersection, FillRule.NonZero, intResNZ);
        PathsD intResNZ_kH = GeoWrangler.makeKeyHole(intResNZ, false, true);
        Assert.AreEqual(0, intResNZ_kH.Count);
        PathsD intResNZc = GeoWrangler.close(intResNZ);
        PathsD intResNZc_kH = GeoWrangler.makeKeyHole(intResNZc, false, true);
        Assert.AreEqual(0, intResNZc_kH.Count);

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
        Assert.AreEqual(1, intRes2_kH.Count);
        Assert.AreEqual(Clipper.Area(intRes2), Clipper.Area(intRes2_kH));

        PathsD intRes2c = GeoWrangler.close(intRes2);
        PathsD intRes2c_kH = GeoWrangler.makeKeyHole(intRes2c, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, intRes2c);
        SvgUtils.AddSolution(svgSrc, intRes2c_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_reversed_khint2c.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(1, intRes2c_kH.Count);
        Assert.AreEqual(Clipper.Area(intRes2c), Clipper.Area(intRes2c_kH));

        // no good - no result
        PathsD intRes2P = new();
        c.Execute(ClipType.Intersection, FillRule.Positive, intRes2P);
        Assert.AreEqual(0, intRes2P.Count);

        // No results - no geometry
        PathsD intRes2P_kH = GeoWrangler.makeKeyHole(intRes2P, false, true);
        Assert.AreEqual(0, intRes2P_kH.Count);

        PathsD intRes2Pc = GeoWrangler.close(intRes2P);
        PathsD intRes2Pc_kH = GeoWrangler.makeKeyHole(intRes2Pc, false, true);
        Assert.AreEqual(0, intRes2Pc_kH.Count);

        // seems good - get keyhole
        PathsD intRes2NZ = new();
        c.Execute(ClipType.Intersection, FillRule.NonZero, intRes2NZ);
        Assert.AreEqual(2, intRes2NZ.Count);

        PathsD intRes2NZ_kH = GeoWrangler.makeKeyHole(intRes2NZ, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, intRes2NZ);
        SvgUtils.AddSolution(svgSrc, intRes2NZ_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_reversed_khint2nz.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(1, intRes2NZ_kH.Count);
        Assert.AreEqual(-29994.9975, Clipper.Area(intRes2NZ_kH), 0.001);

        PathsD intRes2NZc = GeoWrangler.close(intRes2NZ);
        PathsD intRes2NZc_kH = GeoWrangler.makeKeyHole(intRes2NZc, false, true);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, intRes2NZc);
        SvgUtils.AddSolution(svgSrc, intRes2NZc_kH, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "selfoverlaptest_reversed_khint2nzc.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(1, intRes2NZc_kH.Count);
        Assert.AreEqual(-29994.9975, Clipper.Area(intRes2NZc_kH), 0.001);
    }

    private static void comboTest()
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
        Assert.AreEqual(1, kH.Count);
        Assert.AreEqual(-180009.9975, Clipper.Area(kH), 0.001);
        Assert.AreEqual(180020, Clipper.Area(kHSource));
        
        // Gap removal test
        PathsD gR = GeoWrangler.gapRemoval(kH, 100);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, kH);
        SvgUtils.AddSolution(svgSrc, gR, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "combotest_gr.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(2, gR.Count);
        Assert.AreEqual(Math.Abs(Clipper.Area(kHSource)), Math.Abs(Clipper.Area(gR)));

        // Sliver removal test
        PathsD sR = GeoWrangler.gapRemoval(kH, -100);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, kH);
        SvgUtils.AddSolution(svgSrc, sR, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "combotest_sr.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(1, sR.Count);
        Assert.AreEqual(-179989.9975, Clipper.Area(sR), 0.001);
    }

    private static void simple_islandTest()
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
        Assert.AreEqual(2, kH.Count);
        Assert.AreEqual(-239979.995, Clipper.Area(kH), 0.001);
        Assert.AreEqual(240000, Clipper.Area(kHSource));

        // Gap removal test
        PathsD gR = GeoWrangler.gapRemoval(kH, 100);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, kH);
        SvgUtils.AddSolution(svgSrc, gR, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "simpleislandtest_gr.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(4, gR.Count);
        Assert.AreEqual(Math.Abs(Clipper.Area(kHSource)), Math.Abs(Clipper.Area(gR)));
        
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
        Assert.AreEqual(2, sR.Count);
        Assert.AreEqual(80000, Clipper.Area(sR));
    }

    private static void complex_islandTest()
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
        Assert.AreEqual(1, kH.Count);
        Assert.AreEqual(-180009.9975, area_kh, 0.001);
        Assert.AreEqual(180020, Clipper.Area(kHSource));

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
        Assert.AreEqual(1, kH2.Count);
        Assert.AreEqual(-119989.9975, area_kh2, 0.001);
        Assert.AreEqual(120000, Clipper.Area(kHSource));

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
        Assert.AreEqual(1, kH3.Count);
        Assert.AreEqual(-199979.995, area_kh3, 0.001);
        Assert.AreEqual(200000, Clipper.Area(kHSource));

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
        Assert.AreEqual(2, kH12.Count);
        Assert.AreEqual(area_kh + area_kh2, Clipper.Area(kH12));
        Assert.AreEqual(300020, Clipper.Area(kHSource));
        
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
        Assert.AreEqual(2, kH23h1.Count);
        double area_ = Clipper.Area(kH23h1);
        Assert.AreEqual(-339979.995, Clipper.Area(kH23h1), 0.001);
        Assert.AreEqual(340000, Clipper.Area(kHSource));

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
        Assert.AreEqual(2, kH23h2.Count);
        Assert.AreEqual(-339979.995, Clipper.Area(kH23h2), 0.001);
        Assert.AreEqual(340000, Clipper.Area(kHSource));

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
        Assert.AreEqual(2, kH13.Count);
        Assert.AreEqual(area_kh + area_kh3, Clipper.Area(kH13));
        Assert.AreEqual(380020, Clipper.Area(kHSource));

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
        Assert.AreEqual(kH123.Count, 3);
        Assert.AreEqual(Math.Abs(Clipper.Area(kH123)), Math.Abs(area_kh) + Math.Abs(area_kh2) + Math.Abs(area_kh3));
        Assert.AreEqual(500020, Clipper.Area(kHSource));
        
        // Gap removal test
        PathsD gR = GeoWrangler.gapRemoval(kH123, 100);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, kH123);
        SvgUtils.AddSolution(svgSrc, gR, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "complexislandtest_gr.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(gR.Count, 7);
        Assert.AreEqual(Clipper.Area(kHSource), -Clipper.Area(gR));
        
        // Sliver removal test
        PathsD sR = GeoWrangler.gapRemoval(kH123, -100);
        svgSrc.ClearAll();
        SvgUtils.AddClip(svgSrc, kH123);
        SvgUtils.AddSolution(svgSrc, sR, true);
        SvgUtils.SaveToFile(svgSrc, root_loc + "complexislandtest_sr.svg", FillRule.NonZero, 800, 800, 10);
        Assert.AreEqual(sR.Count, 3);
        Assert.AreEqual(-499959.989, Clipper.Area(sR), 0.001);
    }

    private static void multiHoleTest()
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
        Assert.AreEqual(1, kH.Count);
        Assert.AreEqual(-4987.99, area_kh, 0.001);
        Assert.AreEqual(-4990, Clipper.Area(kHSource));
    }
}