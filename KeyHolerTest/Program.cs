using Clipper2Lib;
using geoWrangler;
using System.Collections.Generic;

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

    private static void Main(string[] args)
    {
        singleTest();
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
        Path64 outer = new()
        {
            new(-200000, -200000),
            new(200000, -200000),
            new(200000, 200000),
            new(-200000, 200000),
            new(-200000, -200000)
        };

        Path64 inner1 = new()
        {
            new(-100000, -100000),
            new(-100000, 100000),
            new(100000, 100000),
            new(100000, -100000),
            new(-100000, -100000)
        };

        // Segment the paths to match real-world case.
        Fragmenter f = new(10000);
        Path64 outer_f = f.fragmentPath(outer);

        Path64 inner1_f = f.fragmentPath(inner1);

        Paths64 kHSource = new()
        {
            outer,
            inner1
        };

        // Generate keyholed geometry
        Paths64 kH = GeoWrangler.makeKeyHole(kHSource);

        // Generate sliver geometry.
        Paths64 sL = new();
        Clipper c = new();
        c.AddSubject(outer);
        c.AddClip(kH);
        c.Execute(ClipType.Difference, FillRule.EvenOdd, sL);

        // Gap removal test
        Paths64 gR = GeoWrangler.gapRemoval(kH, 100);

        // Sliver removal test
        Paths64 sR = GeoWrangler.gapRemoval(sL, -100);
    }

    private static void multiTest()
    {
        Path64 outer = new()
        {
            new(-300000, -200000),
            new(300000, -200000),
            new(300000, 200000),
            new(-300000, 200000),
            new(-300000, -200000)
        };

        Path64 inner1 = new()
        {
            new(-200000, -100000),
            new(-200000, 100000),
            new(-100000, 100000),
            new(-100000, -100000),
            new(-200000, -100000)
        };

        Path64 inner2 = new()
        {
            new(100000, -100000),
            new(100000, 100000),
            new(200000, 100000),
            new(200000, -100000),
            new(100000, -100000)
        };

        // Segment the paths to match real-world case.
        Fragmenter f = new(10000);
        Path64 outer_f = f.fragmentPath(outer);

        Path64 inner1_f = f.fragmentPath(inner1);
        Path64 inner2_f = f.fragmentPath(inner2);

        Paths64 kHSource = new()
        {
            outer,
            inner1,
            inner2
        };

        // Generate keyholed geometry
        Paths64 kH = GeoWrangler.makeKeyHole(new Paths64(kHSource));

        // Generate sliver geometry.
        Paths64 sL = new();
        Clipper c = new();
        c.AddSubject(outer);
        c.AddClip(kH);
        c.Execute(ClipType.Difference, FillRule.EvenOdd, sL);

        // Gap removal test
        Paths64 gR = GeoWrangler.gapRemoval(kH, 100);

        // Sliver removal test
        Paths64 sR = GeoWrangler.gapRemoval(sL, -100);
    }

    private static void multiCutTest()
    {
        Path64 outer = new()
        {
            new(0, 0),
            new(400000, 00000),
            new(400000, 400000),
            new(0, 400000),
            new(0, 0)
        };

        Path64 inner1 = new()
        {
            new(50000, 150000),
            new(50000, 250000),
            new(350000, 250000),
            new(350000, 150000),
            new(50000, 150000)
        };

        Path64 inner2 = new()
        {
            new(150000, 50000),
            new(150000, 350000),
            new(250000, 350000),
            new(250000, 50000),
            new(150000, 50000)
        };

        Paths64 kHSource = new()
        {
            outer,
            inner1
        };

        // Generate keyholed geometry
        Paths64 kH = GeoWrangler.makeKeyHole(new Paths64(kHSource));

        Paths64 kHSource2 = new();
        kHSource2.AddRange(kH);
        kHSource2.Add(inner2);

        // Generate keyholed geometry
        Paths64 kH2 = GeoWrangler.makeKeyHole(new Paths64(kHSource2));

        // Generate sliver geometry.
        Clipper c = new();
        Paths64 sL = new();
        c.AddSubject(outer);
        c.AddClip(kH);
        c.Execute(ClipType.Difference, FillRule.EvenOdd, sL);

        c.Clear();
        Paths64 sL2 = new();
        c.AddSubject(outer);
        c.AddClip(kH2);
        c.Execute(ClipType.Difference, FillRule.EvenOdd, sL2);

        // Gap removal test
        Paths64 gR = GeoWrangler.gapRemoval(kH, 100);
        Paths64 gR2 = GeoWrangler.gapRemoval(kH2, 100);

        Paths64 kHSource3 = new();
        kHSource3.AddRange(gR);

        Paths64 kHSource4 = new();
        kHSource4.AddRange(gR2);

        // Generate keyholed geometry
        Paths64 kH3 = GeoWrangler.makeKeyHole(new Paths64(kHSource3));
        Paths64 kH4 = GeoWrangler.makeKeyHole(new Paths64(kHSource4));

        // Sliver removal test
        Paths64 sR = GeoWrangler.gapRemoval(sL, -100);
        Paths64 sR2 = GeoWrangler.gapRemoval(sL2, -100);
    }

    private static void selfOverlapTest()
    {
        Path64 outer = new()
        {
            new(0, 0),
            new(110000, 0),
            new(110000, 50000),
            new(50000, 50000),
            new(50000, 150000),
            new(150000, 150000),
            new(150000, 50000),
            new(90000, 50000),
            new(90000, 0),
            new(200000, 0),
            new(200000, 200000),
            new(0, 200000),
            new(0, 0)
        };

        // decomposer test
        Paths64 dSource = new() {outer};
        Paths64 decomp = GeoWrangler.decompose(dSource);
        Paths64 kHD = GeoWrangler.makeKeyHole(dSource);

        // keyholer test
        Paths64 kHSource = new() {outer};
        Paths64 kH = GeoWrangler.makeKeyHole(kHSource);

        Clipper c = new();
        c.AddSubject(outer);

        // no good
        Paths64 unionRes = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, unionRes);
        Paths64 unionRes_kH = GeoWrangler.makeKeyHole(unionRes);
        Paths64 unionResc = GeoWrangler.close(unionRes);
        Paths64 unionResc_kH = GeoWrangler.makeKeyHole(unionResc);

        // seems good - get keyhole
        Paths64 unionResP = new();
        c.Execute(ClipType.Union, FillRule.Positive, unionResP);
        Paths64 unionResP_kH = GeoWrangler.makeKeyHole(unionResP);
        Paths64 unionResPc = GeoWrangler.close(unionResP);
        Paths64 unionResPc_kH = GeoWrangler.makeKeyHole(unionResPc);

        // seems good - get keyhole
        Paths64 unionResNZ = new();
        c.Execute(ClipType.Union, FillRule.NonZero, unionResNZ);
        Paths64 unionResNZ_kH = GeoWrangler.makeKeyHole(unionResNZ);
        Paths64 unionResNZc = GeoWrangler.close(unionResNZ);
        Paths64 unionResNZc_kH = GeoWrangler.makeKeyHole(unionResNZc);

        // no good - no result
        Paths64 simplifyRes = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, simplifyRes);
        simplifyRes = GeoWrangler.stripColinear(simplifyRes);
        Paths64 simplifyRes_kH = GeoWrangler.makeKeyHole(simplifyRes);
        Paths64 simplifyResc = GeoWrangler.close(simplifyRes);
        Paths64 simplifyResc_kH = GeoWrangler.makeKeyHole(simplifyResc);

        Paths64 simplifyRes2 = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, simplifyRes2);
        Paths64 simplifyRes2_kH = GeoWrangler.makeKeyHole(simplifyRes2);
        Paths64 simplifyRes2c = GeoWrangler.close(simplifyRes2);
        Paths64 simplifyRes2c_kH = GeoWrangler.makeKeyHole(simplifyRes2c);

        // no good - no result
        Paths64 intRes = new();
        c.Execute(ClipType.Intersection, FillRule.EvenOdd, intRes);
        Paths64 intRes_kH = GeoWrangler.makeKeyHole(intRes);
        Paths64 intResc = GeoWrangler.close(intRes);
        Paths64 intResc_kH = GeoWrangler.makeKeyHole(intResc);

        // no good - no result
        Paths64 intResP = new();
        c.Execute(ClipType.Intersection, FillRule.Positive, intResP);
        Paths64 intResP_kH = GeoWrangler.makeKeyHole(intResP);
        Paths64 intResPc = GeoWrangler.close(intResP);
        Paths64 intResPc_kH = GeoWrangler.makeKeyHole(intResPc);

        // no good - no result
        Paths64 intResNZ = new();
        c.Execute(ClipType.Intersection, FillRule.NonZero, intResNZ);
        Paths64 intResNZ_kH = GeoWrangler.makeKeyHole(intResNZ);
        Paths64 intResNZc = GeoWrangler.close(intResNZ);
        Paths64 intResNZc_kH = GeoWrangler.makeKeyHole(intResNZc);

        Rect64 bounds = ClipperFunc.GetBounds(new Paths64 { outer });
        Path64 bb = new()
        {
            new(bounds.left, bounds.bottom),
            new(bounds.left, bounds.top),
            new(bounds.right, bounds.top),
            new(bounds.right, bounds.bottom)
        };

        c.Clear();
        c.AddSubject(bb);
        c.AddClip(outer);

        // no good - overlap region is a gap.
        Paths64 intRes2 = new();
        c.Execute(ClipType.Intersection, FillRule.EvenOdd, intRes2);
        Paths64 intRes2_kH = GeoWrangler.makeKeyHole(intRes2);
        Paths64 intRes2c = GeoWrangler.close(intRes2);
        Paths64 intRes2c_kH = GeoWrangler.makeKeyHole(intRes2c);

        // seems good - get keyhole
        Paths64 intRes2P = new();
        c.Execute(ClipType.Intersection, FillRule.Positive, intRes2P);
        Paths64 intRes2P_kH = GeoWrangler.makeKeyHole(intRes2P);
        Paths64 intRes2Pc = GeoWrangler.close(intRes2P);
        Paths64 intRes2Pc_kH = GeoWrangler.makeKeyHole(intRes2Pc);

        // seems good - get keyhole
        Paths64 intRes2NZ = new();
        c.Execute(ClipType.Intersection, FillRule.NonZero, intRes2NZ);
        Paths64 intRes2NZ_kH = GeoWrangler.makeKeyHole(intRes2NZ);
        Paths64 intRes2NZc = GeoWrangler.close(intRes2NZ);
        Paths64 intRes2NZc_kH = GeoWrangler.makeKeyHole(intRes2NZc);
    }

    private static void selfOverlapTest_reversed()
    {
        Path64 outer = new()
        {
            new(0, 0),
            new(110000, 0),
            new(110000, 50000),
            new(50000, 50000),
            new(50000, 150000),
            new(150000, 150000),
            new(150000, 50000),
            new(90000, 50000),
            new(90000, 0),
            new(200000, 0),
            new(200000, 200000),
            new(0, 200000),
            new(0, 0)
        };

        outer.Reverse();

        // decomposer test
        Paths64 dSource = new() {outer};
        Paths64 decomp = GeoWrangler.decompose(dSource);
        Paths64 kHD = GeoWrangler.makeKeyHole(dSource);

        // keyholer test
        Paths64 kHSource = new() {outer};
        Paths64 kH = GeoWrangler.makeKeyHole(kHSource);

        Clipper c = new();
        c.AddSubject(outer);

        // no good - overlap region is a gap.
        Paths64 unionRes = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, unionRes);
        Paths64 unionRes_kH = GeoWrangler.makeKeyHole(unionRes);
        Paths64 unionResc = GeoWrangler.close(unionRes);
        Paths64 unionResc_kH = GeoWrangler.makeKeyHole(unionResc);

        // no good - no result
        Paths64 unionResP = new();
        c.Execute(ClipType.Union, FillRule.Positive, unionResP);
        Paths64 unionResP_kH = GeoWrangler.makeKeyHole(unionResP);
        Paths64 unionResPc = GeoWrangler.close(unionResP);
        Paths64 unionResPc_kH = GeoWrangler.makeKeyHole(unionResPc);

        // seems good - get keyhole
        Paths64 unionResNZ = new();
        c.Execute(ClipType.Union, FillRule.NonZero, unionResNZ);
        Paths64 unionResNZ_kH = GeoWrangler.makeKeyHole(unionResNZ);
        Paths64 unionResNZc = GeoWrangler.close(unionResNZ);
        Paths64 unionResNZc_kH = GeoWrangler.makeKeyHole(unionResNZc);

        // no good - overlap region is a gap.
        Paths64 simplifyRes = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, simplifyRes);
        simplifyRes = GeoWrangler.stripColinear(simplifyRes);
        Paths64 simplifyRes_kH = GeoWrangler.makeKeyHole(simplifyRes);
        Paths64 simplifyResc = GeoWrangler.close(simplifyRes);
        Paths64 simplifyResc_kH = GeoWrangler.makeKeyHole(simplifyResc);

        Paths64 simplifyRes2 = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, simplifyRes2);
        Paths64 simplifyRes2_kH = GeoWrangler.makeKeyHole(simplifyRes2);
        Paths64 simplifyRes2c = GeoWrangler.close(simplifyRes2);
        Paths64 simplifyRes2c_kH = GeoWrangler.makeKeyHole(simplifyRes2c);

        // no good - no result
        Paths64 intRes = new();
        c.Execute(ClipType.Intersection, FillRule.EvenOdd, intRes);
        Paths64 intRes_kH = GeoWrangler.makeKeyHole(intRes);
        Paths64 intResc = GeoWrangler.close(intRes);
        Paths64 intResc_kH = GeoWrangler.makeKeyHole(intResc);

        // no good - no result
        Paths64 intResP = new();
        c.Execute(ClipType.Intersection, FillRule.Positive, intResP);
        Paths64 intResP_kH = GeoWrangler.makeKeyHole(intResP);
        Paths64 intResPc = GeoWrangler.close(intResP);
        Paths64 intResPc_kH = GeoWrangler.makeKeyHole(intResPc);

        // no good - no result
        Paths64 intResNZ = new();
        c.Execute(ClipType.Intersection, FillRule.NonZero, intResNZ);
        Paths64 intResNZ_kH = GeoWrangler.makeKeyHole(intResNZ);
        Paths64 intResNZc = GeoWrangler.close(intResNZ);
        Paths64 intResNZc_kH = GeoWrangler.makeKeyHole(intResNZc);

        Rect64 bounds = ClipperFunc.GetBounds(new Paths64 { outer });
        Path64 bb = new()
        {
            new(bounds.left, bounds.bottom),
            new(bounds.left, bounds.top),
            new(bounds.right, bounds.top),
            new(bounds.right, bounds.bottom)
        };

        c.Clear();
        c.AddSubject(bb);
        c.AddClip(outer);

        // no good - overlap region is a gap.
        Paths64 intRes2 = new();
        c.Execute(ClipType.Intersection, FillRule.EvenOdd, intRes2);
        Paths64 intRes2_kH = GeoWrangler.makeKeyHole(intRes2);
        Paths64 intRes2c = GeoWrangler.close(intRes2);
        Paths64 intRes2c_kH = GeoWrangler.makeKeyHole(intRes2c);

        // no good - no result
        Paths64 intRes2P = new();
        c.Execute(ClipType.Intersection, FillRule.Positive, intRes2P);
        Paths64 intRes2P_kH = GeoWrangler.makeKeyHole(intRes2P);
        Paths64 intRes2Pc = GeoWrangler.close(intRes2P);
        Paths64 intRes2Pc_kH = GeoWrangler.makeKeyHole(intRes2Pc);

        // seems good - get keyhole
        Paths64 intRes2NZ = new();
        c.Execute(ClipType.Intersection, FillRule.NonZero, intRes2NZ);
        Paths64 intRes2NZ_kH = GeoWrangler.makeKeyHole(intRes2NZ);
        Paths64 intRes2NZc = GeoWrangler.close(intRes2NZ);
        Paths64 intRes2NZc_kH = GeoWrangler.makeKeyHole(intRes2NZc);
    }

    private static void comboTest()
    {
        // Manually create the sliver.
        Path64 outer = new()
        {
            new(-200000, -200000),
            new(300000, -200000),
            new(300000, 200000),
            new(-200000, 200000),
            new(-200000, -99900),
            new(-300000, -99900),
            new(-300000, -100100),
            new(-200000, -100100),
            new(-200000, -200000)
        };

        Path64 inner1 = new()
        {
            new(100000, -100000),
            new(100000, 100000),
            new(200000, 100000),
            new(200000, -100000),
            new(100000, -100000)
        };

        // Segment the paths to match real-world case.
        Fragmenter f = new(10000);
        Path64 outer_f = f.fragmentPath(outer);

        Path64 inner1_f = f.fragmentPath(inner1);

        Paths64 kHSource = new()
        {
            outer,
            inner1
        };

        Paths64 kH = GeoWrangler.makeKeyHole(kHSource);

        // Gap removal test
        Paths64 gR = GeoWrangler.gapRemoval(kH, 100);

        // Sliver removal test
        Paths64 sR = GeoWrangler.gapRemoval(kH, -100);
    }

    private static void simple_islandTest()
    {
        Path64 outer1 = new()
        {
            new(-200000, -200000),
            new(200000, -200000),
            new(200000, 200000),
            new(-200000, 200000),
            new(-200000, -200000)
        };

        Path64 inner1 = new()
        {
            new(-100000, -100000),
            new(-100000, 100000),
            new(100000, 100000),
            new(100000, -100000),
            new(-100000, -100000)
        };

        Path64 outer2 = new()
        {
            new(-200000, 400000),
            new(200000, 400000),
            new(200000, 800000),
            new(-200000, 800000),
            new(-200000, 400000)
        };

        Path64 inner2 = new()
        {
            new(-100000, 500000),
            new(-100000, 700000),
            new(100000, 700000),
            new(100000, 500000),
            new(-100000, 500000)
        };

        Paths64 kHSource = new()
        {
            outer1,
            inner1,
            outer2,
            inner2
        };

        // Generate keyholed geometry
        Paths64 kH = GeoWrangler.makeKeyHole(kHSource);

        // Generate sliver geometry.
        Paths64 sL = new();
        Clipper c = new();
        c.AddSubject(outer1);
        c.AddSubject(outer2);
        c.AddClip(kH);
        c.Execute(ClipType.Difference, FillRule.EvenOdd, sL);

        // Gap removal test
        Paths64 gR = GeoWrangler.gapRemoval(kH, 100);

        // Sliver removal test
        Paths64 sR = GeoWrangler.gapRemoval(sL, -100);
    }

    private static void complex_islandTest()
    {
        // Island 1 - mix of sliver and gap at the end.
        // Manually create the sliver.
        Path64 outer = new()
        {
            new(-200000, -200000),
            new(300000, -200000),
            new(300000, 200000),
            new(-200000, 200000),
            new(-200000, -99900),
            new(-300000, -99900),
            new(-300000, -100100),
            new(-200000, -100100),
            new(-200000, -200000)
        };

        Path64 inner1 = new()
        {
            new(100000, -100000),
            new(100000, 100000),
            new(200000, 100000),
            new(200000, -100000),
            new(100000, -100000)
        };

        Paths64 kHSource = new()
        {
            outer,
            inner1
        };

        // Island 2 - simple single hole
        Path64 outer2 = new()
        {
            new(-200000, 400000),
            new(200000, 400000),
            new(200000, 800000),
            new(-200000, 800000),
            new(-200000, 400000)
        };

        Path64 inner2 = new()
        {
            new(-100000, 500000),
            new(-100000, 700000),
            new(100000, 700000),
            new(100000, 500000),
            new(-100000, 500000)
        };

        kHSource.Add(outer2);
        kHSource.Add(inner2);

        // Island 3 - dual hole hole
        Path64 outer3 = new()
        {
            new(-300000, 1200000),
            new(300000, 1200000),
            new(300000, 1600000),
            new(-300000, 1600000),
            new(-300000, 1200000)
        };

        Path64 inner3_1 = new()
        {
            new(-200000, 1300000),
            new(-200000, 1500000),
            new(-100000, 1500000),
            new(-100000, 1300000),
            new(-200000, 1300000)
        };

        Path64 inner3_2 = new()
        {
            new(100000, 1300000),
            new(100000, 1500000),
            new(200000, 1500000),
            new(200000, 1300000),
            new(100000, 1300000)
        };

        kHSource.Add(outer3);
        kHSource.Add(inner3_1);
        kHSource.Add(inner3_2);

        Paths64 kH = GeoWrangler.makeKeyHole(kHSource);

        // Gap removal test
        Paths64 gR = GeoWrangler.gapRemoval(kH, 100);

        // Sliver removal test
        Paths64 sR = GeoWrangler.gapRemoval(kH, -100);
    }

    private static void multiHoleTest()
    {
        Paths64 paths = new();
        Path64 path_1 = new()
        {
            new(50000, 0),
            new(0, 0),
            new(0, 155000),
            new(50000, 155000),
            new(50000, 0)
        };
        paths.Add(path_1);

        Path64 path_2 = new()
        {
            new(35000, 135000),
            new(35000, 150000),
            new(5000, 150000),
            new(5000, 135000),
            new(35000, 135000)
        };
        paths.Add(path_2);

        Path64 path_3 = new()
        {
            new(22000, 95000),
            new(22000, 125000),
            new(5000, 125000),
            new(5000, 95000),
            new(22000, 95000)
        };
        paths.Add(path_3);

        Path64 path_4 = new()
        {
            new(35000, 45000),
            new(35000, 75000),
            new(5000, 75000),
            new(5000, 45000),
            new(35000, 45000)
        };
        paths.Add(path_4);

        Path64 path_5 = new()
        {
            new(35000, 5000),
            new(35000, 35000),
            new(5000, 35000),
            new(5000, 5000),
            new(35000, 5000)
        };
        paths.Add(path_5);

        // Generate keyholed geometry
        Paths64 kH = GeoWrangler.makeKeyHole(paths);
    }
}