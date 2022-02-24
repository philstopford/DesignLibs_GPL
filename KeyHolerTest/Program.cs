using ClipperLib1;
using geoWrangler;
using System.Collections.Generic;

namespace KeyHolerTest;

using Path = List<IntPoint>;
using Paths = List<List<IntPoint>>;

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
        Path outer = new()
        {
            new IntPoint(-200000, -200000),
            new IntPoint(200000, -200000),
            new IntPoint(200000, 200000),
            new IntPoint(-200000, 200000),
            new IntPoint(-200000, -200000)
        };

        Path inner1 = new()
        {
            new IntPoint(-100000, -100000),
            new IntPoint(-100000, 100000),
            new IntPoint(100000, 100000),
            new IntPoint(100000, -100000),
            new IntPoint(-100000, -100000)
        };

        // Segment the paths to match real-world case.
        Fragmenter f = new(10000);
        Path outer_f = f.fragmentPath(outer);

        Path inner1_f = f.fragmentPath(inner1);

        Paths kHSource = new()
        {
            outer,
            inner1
        };

        // Generate keyholed geometry
        Paths kH = GeoWrangler.makeKeyHole(kHSource);

        // Generate sliver geometry.
        Paths sL = new();
        Clipper c = new();
        c.AddPath(outer, PolyType.ptSubject, true);
        c.AddPaths(kH, PolyType.ptClip, true);
        c.Execute(ClipType.ctDifference, sL);

        // Gap removal test
        Paths gR = GeoWrangler.gapRemoval(kH, 100);

        // Sliver removal test
        Paths sR = GeoWrangler.gapRemoval(sL, -100);
    }

    private static void multiTest()
    {
        Path outer = new()
        {
            new IntPoint(-300000, -200000),
            new IntPoint(300000, -200000),
            new IntPoint(300000, 200000),
            new IntPoint(-300000, 200000),
            new IntPoint(-300000, -200000)
        };

        Path inner1 = new()
        {
            new IntPoint(-200000, -100000),
            new IntPoint(-200000, 100000),
            new IntPoint(-100000, 100000),
            new IntPoint(-100000, -100000),
            new IntPoint(-200000, -100000)
        };

        Path inner2 = new()
        {
            new IntPoint(100000, -100000),
            new IntPoint(100000, 100000),
            new IntPoint(200000, 100000),
            new IntPoint(200000, -100000),
            new IntPoint(100000, -100000)
        };

        // Segment the paths to match real-world case.
        Fragmenter f = new(10000);
        Path outer_f = f.fragmentPath(outer);

        Path inner1_f = f.fragmentPath(inner1);
        Path inner2_f = f.fragmentPath(inner2);

        Paths kHSource = new()
        {
            outer,
            inner1,
            inner2
        };

        // Generate keyholed geometry
        Paths kH = GeoWrangler.makeKeyHole(new Paths(kHSource));

        // Generate sliver geometry.
        Paths sL = new();
        Clipper c = new();
        c.AddPath(outer, PolyType.ptSubject, true);
        c.AddPaths(kH, PolyType.ptClip, true);
        c.Execute(ClipType.ctDifference, sL);

        // Gap removal test
        Paths gR = GeoWrangler.gapRemoval(kH, 100);

        // Sliver removal test
        Paths sR = GeoWrangler.gapRemoval(sL, -100);
    }

    private static void multiCutTest()
    {
        Path outer = new()
        {
            new IntPoint(0, 0),
            new IntPoint(400000, 00000),
            new IntPoint(400000, 400000),
            new IntPoint(0, 400000),
            new IntPoint(0, 0)
        };

        Path inner1 = new()
        {
            new IntPoint(50000, 150000),
            new IntPoint(50000, 250000),
            new IntPoint(350000, 250000),
            new IntPoint(350000, 150000),
            new IntPoint(50000, 150000)
        };

        Path inner2 = new()
        {
            new IntPoint(150000, 50000),
            new IntPoint(150000, 350000),
            new IntPoint(250000, 350000),
            new IntPoint(250000, 50000),
            new IntPoint(150000, 50000)
        };

        Paths kHSource = new()
        {
            outer,
            inner1
        };

        // Generate keyholed geometry
        Paths kH = GeoWrangler.makeKeyHole(new Paths(kHSource));

        Paths kHSource2 = new();
        kHSource2.AddRange(kH);
        kHSource2.Add(inner2);

        // Generate keyholed geometry
        Paths kH2 = GeoWrangler.makeKeyHole(new Paths(kHSource2));

        // Generate sliver geometry.
        Clipper c = new();
        Paths sL = new();
        c.AddPath(outer, PolyType.ptSubject, true);
        c.AddPaths(kH, PolyType.ptClip, true);
        c.Execute(ClipType.ctDifference, sL);

        c.Clear();
        Paths sL2 = new();
        c.AddPath(outer, PolyType.ptSubject, true);
        c.AddPaths(kH2, PolyType.ptClip, true);
        c.Execute(ClipType.ctDifference, sL2);

        // Gap removal test
        Paths gR = GeoWrangler.gapRemoval(kH, 100);
        Paths gR2 = GeoWrangler.gapRemoval(kH2, 100);

        Paths kHSource3 = new();
        kHSource3.AddRange(gR);

        Paths kHSource4 = new();
        kHSource4.AddRange(gR2);

        // Generate keyholed geometry
        Paths kH3 = GeoWrangler.makeKeyHole(new Paths(kHSource3));
        Paths kH4 = GeoWrangler.makeKeyHole(new Paths(kHSource4));

        // Sliver removal test
        Paths sR = GeoWrangler.gapRemoval(sL, -100);
        Paths sR2 = GeoWrangler.gapRemoval(sL2, -100);
    }

    private static void selfOverlapTest()
    {
        Path outer = new()
        {
            new IntPoint(0, 0),
            new IntPoint(110000, 0),
            new IntPoint(110000, 50000),
            new IntPoint(50000, 50000),
            new IntPoint(50000, 150000),
            new IntPoint(150000, 150000),
            new IntPoint(150000, 50000),
            new IntPoint(90000, 50000),
            new IntPoint(90000, 0),
            new IntPoint(200000, 0),
            new IntPoint(200000, 200000),
            new IntPoint(0, 200000),
            new IntPoint(0, 0)
        };

        // decomposer test
        Paths dSource = new() {outer};
        Paths decomp = GeoWrangler.decompose(dSource);
        Paths kHD = GeoWrangler.makeKeyHole(dSource);

        // keyholer test
        Paths kHSource = new() {outer};
        Paths kH = GeoWrangler.makeKeyHole(kHSource);

        Clipper c = new()
        {
            PreserveCollinear = true
        };
        c.AddPath(outer, PolyType.ptSubject, true);

        // no good
        Paths unionRes = new();
        c.Execute(ClipType.ctUnion, unionRes);
        Paths unionRes_kH = GeoWrangler.makeKeyHole(unionRes);
        Paths unionResc = GeoWrangler.close(unionRes);
        Paths unionResc_kH = GeoWrangler.makeKeyHole(unionResc);

        // seems good - get keyhole
        Paths unionResP = new();
        c.Execute(ClipType.ctUnion, unionResP, PolyFillType.pftPositive);
        Paths unionResP_kH = GeoWrangler.makeKeyHole(unionResP);
        Paths unionResPc = GeoWrangler.close(unionResP);
        Paths unionResPc_kH = GeoWrangler.makeKeyHole(unionResPc);

        // seems good - get keyhole
        Paths unionResNZ = new();
        c.Execute(ClipType.ctUnion, unionResNZ, PolyFillType.pftNonZero);
        Paths unionResNZ_kH = GeoWrangler.makeKeyHole(unionResNZ);
        Paths unionResNZc = GeoWrangler.close(unionResNZ);
        Paths unionResNZc_kH = GeoWrangler.makeKeyHole(unionResNZc);

        // no good - no result
        Paths simplifyRes = new();
        simplifyRes = Clipper.SimplifyPolygon(outer, preserveColinear: false);
        Paths simplifyRes_kH = GeoWrangler.makeKeyHole(simplifyRes);
        Paths simplifyResc = GeoWrangler.close(simplifyRes);
        Paths simplifyResc_kH = GeoWrangler.makeKeyHole(simplifyResc);

        Paths simplifyRes2 = new();
        simplifyRes2 = Clipper.SimplifyPolygon(outer, preserveColinear: true);
        Paths simplifyRes2_kH = GeoWrangler.makeKeyHole(simplifyRes2);
        Paths simplifyRes2c = GeoWrangler.close(simplifyRes2);
        Paths simplifyRes2c_kH = GeoWrangler.makeKeyHole(simplifyRes2c);

        // no good - no result
        Paths intRes = new();
        c.Execute(ClipType.ctIntersection, intRes);
        Paths intRes_kH = GeoWrangler.makeKeyHole(intRes);
        Paths intResc = GeoWrangler.close(intRes);
        Paths intResc_kH = GeoWrangler.makeKeyHole(intResc);

        // no good - no result
        Paths intResP = new();
        c.Execute(ClipType.ctIntersection, intResP, PolyFillType.pftPositive);
        Paths intResP_kH = GeoWrangler.makeKeyHole(intResP);
        Paths intResPc = GeoWrangler.close(intResP);
        Paths intResPc_kH = GeoWrangler.makeKeyHole(intResPc);

        // no good - no result
        Paths intResNZ = new();
        c.Execute(ClipType.ctIntersection, intResNZ, PolyFillType.pftNonZero);
        Paths intResNZ_kH = GeoWrangler.makeKeyHole(intResNZ);
        Paths intResNZc = GeoWrangler.close(intResNZ);
        Paths intResNZc_kH = GeoWrangler.makeKeyHole(intResNZc);

        IntRect bounds = ClipperBase.GetBounds(new Paths { outer });
        Path bb = new()
        {
            new IntPoint(bounds.left, bounds.bottom),
            new IntPoint(bounds.left, bounds.top),
            new IntPoint(bounds.right, bounds.top),
            new IntPoint(bounds.right, bounds.bottom)
        };

        c.Clear();
        c.AddPath(bb, PolyType.ptSubject, true);
        c.AddPath(outer, PolyType.ptClip, true);

        // no good - overlap region is a gap.
        Paths intRes2 = new();
        c.Execute(ClipType.ctIntersection, intRes2);
        Paths intRes2_kH = GeoWrangler.makeKeyHole(intRes2);
        Paths intRes2c = GeoWrangler.close(intRes2);
        Paths intRes2c_kH = GeoWrangler.makeKeyHole(intRes2c);

        // seems good - get keyhole
        Paths intRes2P = new();
        c.Execute(ClipType.ctIntersection, intRes2P, PolyFillType.pftPositive);
        Paths intRes2P_kH = GeoWrangler.makeKeyHole(intRes2P);
        Paths intRes2Pc = GeoWrangler.close(intRes2P);
        Paths intRes2Pc_kH = GeoWrangler.makeKeyHole(intRes2Pc);

        // seems good - get keyhole
        Paths intRes2NZ = new();
        c.Execute(ClipType.ctIntersection, intRes2NZ, PolyFillType.pftNonZero);
        Paths intRes2NZ_kH = GeoWrangler.makeKeyHole(intRes2NZ);
        Paths intRes2NZc = GeoWrangler.close(intRes2NZ);
        Paths intRes2NZc_kH = GeoWrangler.makeKeyHole(intRes2NZc);
    }

    private static void selfOverlapTest_reversed()
    {
        Path outer = new()
        {
            new IntPoint(0, 0),
            new IntPoint(110000, 0),
            new IntPoint(110000, 50000),
            new IntPoint(50000, 50000),
            new IntPoint(50000, 150000),
            new IntPoint(150000, 150000),
            new IntPoint(150000, 50000),
            new IntPoint(90000, 50000),
            new IntPoint(90000, 0),
            new IntPoint(200000, 0),
            new IntPoint(200000, 200000),
            new IntPoint(0, 200000),
            new IntPoint(0, 0)
        };

        outer.Reverse();

        // decomposer test
        Paths dSource = new() {outer};
        Paths decomp = GeoWrangler.decompose(dSource);
        Paths kHD = GeoWrangler.makeKeyHole(dSource);

        // keyholer test
        Paths kHSource = new() {outer};
        Paths kH = GeoWrangler.makeKeyHole(kHSource);

        Clipper c = new()
        {
            PreserveCollinear = true
        };
        c.AddPath(outer, PolyType.ptSubject, true);

        // no good - overlap region is a gap.
        Paths unionRes = new();
        c.Execute(ClipType.ctUnion, unionRes);
        Paths unionRes_kH = GeoWrangler.makeKeyHole(unionRes);
        Paths unionResc = GeoWrangler.close(unionRes);
        Paths unionResc_kH = GeoWrangler.makeKeyHole(unionResc);

        // no good - no result
        Paths unionResP = new();
        c.Execute(ClipType.ctUnion, unionResP, PolyFillType.pftPositive);
        Paths unionResP_kH = GeoWrangler.makeKeyHole(unionResP);
        Paths unionResPc = GeoWrangler.close(unionResP);
        Paths unionResPc_kH = GeoWrangler.makeKeyHole(unionResPc);

        // seems good - get keyhole
        Paths unionResNZ = new();
        c.Execute(ClipType.ctUnion, unionResNZ, PolyFillType.pftNonZero);
        Paths unionResNZ_kH = GeoWrangler.makeKeyHole(unionResNZ);
        Paths unionResNZc = GeoWrangler.close(unionResNZ);
        Paths unionResNZc_kH = GeoWrangler.makeKeyHole(unionResNZc);

        // no good - overlap region is a gap.
        Paths simplifyRes = new();
        simplifyRes = Clipper.SimplifyPolygon(outer, preserveColinear: false);
        Paths simplifyRes_kH = GeoWrangler.makeKeyHole(simplifyRes);
        Paths simplifyResc = GeoWrangler.close(simplifyRes);
        Paths simplifyResc_kH = GeoWrangler.makeKeyHole(simplifyResc);

        Paths simplifyRes2 = new();
        simplifyRes2 = Clipper.SimplifyPolygon(outer, preserveColinear: true);
        Paths simplifyRes2_kH = GeoWrangler.makeKeyHole(simplifyRes2);
        Paths simplifyRes2c = GeoWrangler.close(simplifyRes2);
        Paths simplifyRes2c_kH = GeoWrangler.makeKeyHole(simplifyRes2c);

        // no good - no result
        Paths intRes = new();
        c.Execute(ClipType.ctIntersection, intRes);
        Paths intRes_kH = GeoWrangler.makeKeyHole(intRes);
        Paths intResc = GeoWrangler.close(intRes);
        Paths intResc_kH = GeoWrangler.makeKeyHole(intResc);

        // no good - no result
        Paths intResP = new();
        c.Execute(ClipType.ctIntersection, intResP, PolyFillType.pftPositive);
        Paths intResP_kH = GeoWrangler.makeKeyHole(intResP);
        Paths intResPc = GeoWrangler.close(intResP);
        Paths intResPc_kH = GeoWrangler.makeKeyHole(intResPc);

        // no good - no result
        Paths intResNZ = new();
        c.Execute(ClipType.ctIntersection, intResNZ, PolyFillType.pftNonZero);
        Paths intResNZ_kH = GeoWrangler.makeKeyHole(intResNZ);
        Paths intResNZc = GeoWrangler.close(intResNZ);
        Paths intResNZc_kH = GeoWrangler.makeKeyHole(intResNZc);

        IntRect bounds = ClipperBase.GetBounds(new Paths { outer });
        Path bb = new()
        {
            new IntPoint(bounds.left, bounds.bottom),
            new IntPoint(bounds.left, bounds.top),
            new IntPoint(bounds.right, bounds.top),
            new IntPoint(bounds.right, bounds.bottom)
        };

        c.Clear();
        c.AddPath(bb, PolyType.ptSubject, true);
        c.AddPath(outer, PolyType.ptClip, true);

        // no good - overlap region is a gap.
        Paths intRes2 = new();
        c.Execute(ClipType.ctIntersection, intRes2);
        Paths intRes2_kH = GeoWrangler.makeKeyHole(intRes2);
        Paths intRes2c = GeoWrangler.close(intRes2);
        Paths intRes2c_kH = GeoWrangler.makeKeyHole(intRes2c);

        // no good - no result
        Paths intRes2P = new();
        c.Execute(ClipType.ctIntersection, intRes2P, PolyFillType.pftPositive);
        Paths intRes2P_kH = GeoWrangler.makeKeyHole(intRes2P);
        Paths intRes2Pc = GeoWrangler.close(intRes2P);
        Paths intRes2Pc_kH = GeoWrangler.makeKeyHole(intRes2Pc);

        // seems good - get keyhole
        Paths intRes2NZ = new();
        c.Execute(ClipType.ctIntersection, intRes2NZ, PolyFillType.pftNonZero);
        Paths intRes2NZ_kH = GeoWrangler.makeKeyHole(intRes2NZ);
        Paths intRes2NZc = GeoWrangler.close(intRes2NZ);
        Paths intRes2NZc_kH = GeoWrangler.makeKeyHole(intRes2NZc);
    }

    private static void comboTest()
    {
        // Manually create the sliver.
        Path outer = new()
        {
            new IntPoint(-200000, -200000),
            new IntPoint(300000, -200000),
            new IntPoint(300000, 200000),
            new IntPoint(-200000, 200000),
            new IntPoint(-200000, -99900),
            new IntPoint(-300000, -99900),
            new IntPoint(-300000, -100100),
            new IntPoint(-200000, -100100),
            new IntPoint(-200000, -200000)
        };

        Path inner1 = new()
        {
            new IntPoint(100000, -100000),
            new IntPoint(100000, 100000),
            new IntPoint(200000, 100000),
            new IntPoint(200000, -100000),
            new IntPoint(100000, -100000)
        };

        // Segment the paths to match real-world case.
        Fragmenter f = new(10000);
        Path outer_f = f.fragmentPath(outer);

        Path inner1_f = f.fragmentPath(inner1);

        Paths kHSource = new()
        {
            outer,
            inner1
        };

        Paths kH = GeoWrangler.makeKeyHole(kHSource);

        // Gap removal test
        Paths gR = GeoWrangler.gapRemoval(kH, 100);

        // Sliver removal test
        Paths sR = GeoWrangler.gapRemoval(kH, -100);
    }

    private static void simple_islandTest()
    {
        Path outer1 = new()
        {
            new IntPoint(-200000, -200000),
            new IntPoint(200000, -200000),
            new IntPoint(200000, 200000),
            new IntPoint(-200000, 200000),
            new IntPoint(-200000, -200000)
        };

        Path inner1 = new()
        {
            new IntPoint(-100000, -100000),
            new IntPoint(-100000, 100000),
            new IntPoint(100000, 100000),
            new IntPoint(100000, -100000),
            new IntPoint(-100000, -100000)
        };

        Path outer2 = new()
        {
            new IntPoint(-200000, 400000),
            new IntPoint(200000, 400000),
            new IntPoint(200000, 800000),
            new IntPoint(-200000, 800000),
            new IntPoint(-200000, 400000)
        };

        Path inner2 = new()
        {
            new IntPoint(-100000, 500000),
            new IntPoint(-100000, 700000),
            new IntPoint(100000, 700000),
            new IntPoint(100000, 500000),
            new IntPoint(-100000, 500000)
        };

        Paths kHSource = new()
        {
            outer1,
            inner1,
            outer2,
            inner2
        };

        // Generate keyholed geometry
        Paths kH = GeoWrangler.makeKeyHole(kHSource);

        // Generate sliver geometry.
        Paths sL = new();
        Clipper c = new();
        c.AddPath(outer1, PolyType.ptSubject, true);
        c.AddPath(outer2, PolyType.ptSubject, true);
        c.AddPaths(kH, PolyType.ptClip, true);
        c.Execute(ClipType.ctDifference, sL);

        // Gap removal test
        Paths gR = GeoWrangler.gapRemoval(kH, 100);

        // Sliver removal test
        Paths sR = GeoWrangler.gapRemoval(sL, -100);
    }

    private static void complex_islandTest()
    {
        // Island 1 - mix of sliver and gap at the end.
        // Manually create the sliver.
        Path outer = new()
        {
            new IntPoint(-200000, -200000),
            new IntPoint(300000, -200000),
            new IntPoint(300000, 200000),
            new IntPoint(-200000, 200000),
            new IntPoint(-200000, -99900),
            new IntPoint(-300000, -99900),
            new IntPoint(-300000, -100100),
            new IntPoint(-200000, -100100),
            new IntPoint(-200000, -200000)
        };

        Path inner1 = new()
        {
            new IntPoint(100000, -100000),
            new IntPoint(100000, 100000),
            new IntPoint(200000, 100000),
            new IntPoint(200000, -100000),
            new IntPoint(100000, -100000)
        };

        Paths kHSource = new()
        {
            outer,
            inner1
        };

        // Island 2 - simple single hole
        Path outer2 = new()
        {
            new IntPoint(-200000, 400000),
            new IntPoint(200000, 400000),
            new IntPoint(200000, 800000),
            new IntPoint(-200000, 800000),
            new IntPoint(-200000, 400000)
        };

        Path inner2 = new()
        {
            new IntPoint(-100000, 500000),
            new IntPoint(-100000, 700000),
            new IntPoint(100000, 700000),
            new IntPoint(100000, 500000),
            new IntPoint(-100000, 500000)
        };

        kHSource.Add(outer2);
        kHSource.Add(inner2);

        // Island 3 - dual hole hole
        Path outer3 = new()
        {
            new IntPoint(-300000, 1200000),
            new IntPoint(300000, 1200000),
            new IntPoint(300000, 1600000),
            new IntPoint(-300000, 1600000),
            new IntPoint(-300000, 1200000)
        };

        Path inner3_1 = new()
        {
            new IntPoint(-200000, 1300000),
            new IntPoint(-200000, 1500000),
            new IntPoint(-100000, 1500000),
            new IntPoint(-100000, 1300000),
            new IntPoint(-200000, 1300000)
        };

        Path inner3_2 = new()
        {
            new IntPoint(100000, 1300000),
            new IntPoint(100000, 1500000),
            new IntPoint(200000, 1500000),
            new IntPoint(200000, 1300000),
            new IntPoint(100000, 1300000)
        };

        kHSource.Add(outer3);
        kHSource.Add(inner3_1);
        kHSource.Add(inner3_2);

        Paths kH = GeoWrangler.makeKeyHole(kHSource);

        // Gap removal test
        Paths gR = GeoWrangler.gapRemoval(kH, 100);

        // Sliver removal test
        Paths sR = GeoWrangler.gapRemoval(kH, -100);
    }

    private static void multiHoleTest()
    {
        Paths paths = new();
        Path path_1 = new()
        {
            new IntPoint(50000, 0),
            new IntPoint(0, 0),
            new IntPoint(0, 155000),
            new IntPoint(50000, 155000),
            new IntPoint(50000, 0)
        };
        paths.Add(path_1);

        Path path_2 = new()
        {
            new IntPoint(35000, 135000),
            new IntPoint(35000, 150000),
            new IntPoint(5000, 150000),
            new IntPoint(5000, 135000),
            new IntPoint(35000, 135000)
        };
        paths.Add(path_2);

        Path path_3 = new()
        {
            new IntPoint(22000, 95000),
            new IntPoint(22000, 125000),
            new IntPoint(5000, 125000),
            new IntPoint(5000, 95000),
            new IntPoint(22000, 95000)
        };
        paths.Add(path_3);

        Path path_4 = new()
        {
            new IntPoint(35000, 45000),
            new IntPoint(35000, 75000),
            new IntPoint(5000, 75000),
            new IntPoint(5000, 45000),
            new IntPoint(35000, 45000)
        };
        paths.Add(path_4);

        Path path_5 = new()
        {
            new IntPoint(35000, 5000),
            new IntPoint(35000, 35000),
            new IntPoint(5000, 35000),
            new IntPoint(5000, 5000),
            new IntPoint(35000, 5000)
        };
        paths.Add(path_5);

        // Generate keyholed geometry
        Paths kH = GeoWrangler.makeKeyHole(paths);
    }
}