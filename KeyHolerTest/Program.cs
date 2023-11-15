using System;
using Clipper2Lib;
using geoWrangler;
using NUnit.Framework;

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
        PathsD kH = GeoWrangler.makeKeyHole(kHSource, reverseEval:false, biDirectionalEval:true);
        Assert.AreEqual(kH.Count, 1);
        // Expected area should be 120000
        double expectedArea = Math.Abs(Clipper.Area(outer)) - Math.Abs(Clipper.Area(inner1));
        // There is a small delta due to the keyhole, but it should be negligible.
        Assert.LessOrEqual(expectedArea - Math.Abs(Clipper.Area(kH)), 10.01);
        // Let's make sure our keyhole didn't move too much as well, that could be annoying.
        Assert.LessOrEqual(Math.Abs(kH[0][2].x - -100.05), 0.001);
        Assert.LessOrEqual(Math.Abs(kH[0][3].x - -100.05), 0.001);
        Assert.LessOrEqual(Math.Abs(kH[0][8].x - -99.95), 0.001);
        Assert.LessOrEqual(Math.Abs(kH[0][9].x - -99.95), 0.001);
        
        // Generate sliver geometry.
        PathsD sL = new();
        ClipperD c = new(Constants.roundingDecimalPrecision);
        c.AddSubject(outer);
        c.AddClip(kH);
        c.Execute(ClipType.Difference, FillRule.EvenOdd, sL);
        Assert.AreEqual(sL.Count, 1);
        Assert.LessOrEqual(Clipper.Area(sL) - 40010.0025, 0.0001);

        // Gap removal test
        PathsD gR = GeoWrangler.gapRemoval(kH, 100);
        Assert.AreEqual(gR.Count, 2);
        Assert.AreEqual(Math.Abs(Clipper.Area(outer)), Math.Abs(Clipper.Area(gR[0])));
        Assert.AreEqual(Math.Abs(Clipper.Area(inner1)), Math.Abs(Clipper.Area(gR[1])));
        Assert.False(Clipper.IsPositive(gR[0]));
        Assert.True(Clipper.IsPositive(gR[1]));
        
        // Sliver removal test
        PathsD sR = GeoWrangler.gapRemoval(sL, -100);
        Assert.AreEqual(Clipper.Area(sR), 40000);
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
        PathsD kH = GeoWrangler.makeKeyHole(new PathsD(kHSource), reverseEval:false, biDirectionalEval:true);
        Assert.AreEqual(kH.Count, 1);
        Assert.LessOrEqual(Clipper.Area(kH) - -199979.995d, 0.001);
        
        // Generate sliver geometry.
        PathsD sL = new();
        ClipperD c = new(Constants.roundingDecimalPrecision);
        c.AddSubject(outer);
        c.AddClip(kH);
        c.Execute(ClipType.Difference, FillRule.EvenOdd, sL);
        
        /* Expected output
        sL = {List<List<Point64>>} Count = 2
         [0] = {List<Point64>} Count = 8
          [0] = {Point64} -300,100.5,0 
          [1] = {Point64} -199.5,100.5,0 
          [2] = {Point64} -199.5,100,0 
          [3] = {Point64} -100,100,0 
          [4] = {Point64} -100,-100,0 
          [5] = {Point64} -200,-100,0 
          [6] = {Point64} -200,99.5,0 
          [7] = {Point64} -300,99.5,0 
         [1] = {List<Point64>} Count = 8
          [0] = {Point64} 199.5,-100.5,0 
          [1] = {Point64} 199.5,-100,0 
          [2] = {Point64} 100,-100,0 
          [3] = {Point64} 100,100,0 
          [4] = {Point64} 200,100,0 
          [5] = {Point64} 200,-99.5,0 
          [6] = {Point64} 300,-99.5,0 
          [7] = {Point64} 300,-100.5,0 
           */

        // Gap removal test
        PathsD gR = GeoWrangler.gapRemoval(kH, 100);

        /* Expected output
        gR = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 21
          [0] = {Point64} 300,-200,0 
          [1] = {Point64} 300,-100.5,0 
          [2] = {Point64} 199.5,-100.5,0 
          [3] = {Point64} 199.5,-100,0 
          [4] = {Point64} 100,-100,0 
          [5] = {Point64} 100,100,0 
          [6] = {Point64} 200,100,0 
          [7] = {Point64} 200,-99.5,0 
          [8] = {Point64} 300,-99.5,0 
          [9] = {Point64} 300,200,0 
          [10] = {Point64} -300,200,0 
          [11] = {Point64} -300,100.5,0 
          [12] = {Point64} -199.5,100.5,0 
          [13] = {Point64} -199.5,100,0 
          [14] = {Point64} -100,100,0 
          [15] = {Point64} -100,-100,0 
          [16] = {Point64} -200,-100,0 
          [17] = {Point64} -200,99.5,0 
          [18] = {Point64} -300,99.5,0 
          [19] = {Point64} -300,-200,0 
          [20] = {Point64} 300,-200,0 
           */
        
        // Sliver removal test
        PathsD sR = GeoWrangler.gapRemoval(sL, -100);
        
        /* Expected output
        sR = {List<List<Point64>>} Count = 2
         [0] = {List<Point64>} Count = 9
          [0] = {Point64} -300,99.5,0 
          [1] = {Point64} -300,100.5,0 
          [2] = {Point64} -199.5,100.5,0 
          [3] = {Point64} -199.5,100,0 
          [4] = {Point64} -100,100,0 
          [5] = {Point64} -100,-100,0 
          [6] = {Point64} -200,-100,0 
          [7] = {Point64} -200,99.5,0 
          [8] = {Point64} -300,99.5,0 
         [1] = {List<Point64>} Count = 9
          [0] = {Point64} 100,-100,0 
          [1] = {Point64} 100,100,0 
          [2] = {Point64} 200,100,0 
          [3] = {Point64} 200,-99.5,0 
          [4] = {Point64} 300,-99.5,0 
          [5] = {Point64} 300,-100.5,0 
          [6] = {Point64} 199.5,-100.5,0 
          [7] = {Point64} 199.5,-100,0 
          [8] = {Point64} 100,-100,0 
           */
    }

    private static void multiCutTest()
    {
        PathD outer = new()
        {
            new(0, 0),
            new(400, 0),
            new(400, 400),
            new(0, 400),
            new(0, 0)
        };

        PathD inner1 = new()
        {
            new(50, 150),
            new(50, 250),
            new(350, 250),
            new(350, 150),
            new(50, 150)
        };

        PathD inner2 = new()
        {
            new(150, 50),
            new(150, 350),
            new(250, 350),
            new(250, 50),
            new(150, 50)
        };

        PathsD kHSource = new()
        {
            outer,
            inner1
        };
        
        /* Expected
        kHSource = {List<List<Point64>>} Count = 2
         [0] = {List<Point64>} Count = 5
          [0] = {Point64} 0,0,0 
          [1] = {Point64} 400,0,0 
          [2] = {Point64} 400,400,0 
          [3] = {Point64} 0,400,0 
          [4] = {Point64} 0,0,0 
         [1] = {List<Point64>} Count = 5
          [0] = {Point64} 50,150,0 
          [1] = {Point64} 50,250,0 
          [2] = {Point64} 350,250,0 
          [3] = {Point64} 350,150,0 
          [4] = {Point64} 50,150,0 
           */

        // Generate keyholed geometry
        PathsD kH = GeoWrangler.makeKeyHole(new PathsD(kHSource), reverseEval:false, biDirectionalEval:true);

        /* Expected output
        kH = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 13
          [0] = {Point64} 400,0,0 
          [1] = {Point64} 400,149.5,0 
          [2] = {Point64} 349.5,149.5,0 
          [3] = {Point64} 349.5,150,0 
          [4] = {Point64} 50,150,0 
          [5] = {Point64} 50,250,0 
          [6] = {Point64} 350,250,0 
          [7] = {Point64} 350,150.5,0 
          [8] = {Point64} 400,150.5,0 
          [9] = {Point64} 400,400,0 
          [10] = {Point64} 0,400,0 
          [11] = {Point64} 0,0,0 
          [12] = {Point64} 400,0,0 
           */
        
        PathsD kHSource2 = new();
        kHSource2.AddRange(kH);
        kHSource2.Add(inner2);

        /* Expected
        kHSource2 = {List<List<Point64>>} Count = 2
         [0] = {List<Point64>} Count = 13
          [0] = {Point64} 400,0,0 
          [1] = {Point64} 400,149.5,0 
          [2] = {Point64} 349.5,149.5,0 
          [3] = {Point64} 349.5,150,0 
          [4] = {Point64} 50,150,0 
          [5] = {Point64} 50,250,0 
          [6] = {Point64} 350,250,0 
          [7] = {Point64} 350,150.5,0 
          [8] = {Point64} 400,150.5,0 
          [9] = {Point64} 400,400,0 
          [10] = {Point64} 0,400,0 
          [11] = {Point64} 0,0,0 
          [12] = {Point64} 400,0,0 
         [1] = {List<Point64>} Count = 5
          [0] = {Point64} 150,50,0 
          [1] = {Point64} 150,350,0 
          [2] = {Point64} 250,350,0 
          [3] = {Point64} 250,50,0 
          [4] = {Point64} 150,50,0 
           */
        
        // Generate keyholed geometry
        PathsD kH2 = GeoWrangler.makeKeyHole(new PathsD(kHSource2), reverseEval:false, biDirectionalEval:true);

        /* Expected output
        kH2 = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 21
          [0] = {Point64} 400,0,0 
          [1] = {Point64} 400,149.5,0 
          [2] = {Point64} 349.5,149.5,0 
          [3] = {Point64} 349.5,150,0 
          [4] = {Point64} 250,150,0 
          [5] = {Point64} 250,50,0 
          [6] = {Point64} 150,50,0 
          [7] = {Point64} 150,150,0 
          [8] = {Point64} 50,150,0 
          [9] = {Point64} 50,250,0 
          [10] = {Point64} 150,250,0 
          [11] = {Point64} 150,350,0 
          [12] = {Point64} 250,350,0 
          [13] = {Point64} 250,250,0 
          [14] = {Point64} 350,250,0 
          [15] = {Point64} 350,150.5,0 
          [16] = {Point64} 400,150.5,0 
          [17] = {Point64} 400,400,0 
          [18] = {Point64} 0,400,0 
          [19] = {Point64} 0,0,0 
          [20] = {Point64} 400,0,0 
           */
        
        // Generate sliver geometry.
        ClipperD c = new(Constants.roundingDecimalPrecision);
        PathsD sL = new();
        c.AddSubject(outer);
        c.AddClip(kH);
        c.Execute(ClipType.Difference, FillRule.EvenOdd, sL);
        
        /* Expected output
        sL = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 8
          [0] = {Point64} 349.5,149.5,0 
          [1] = {Point64} 349.5,150,0 
          [2] = {Point64} 50,150,0 
          [3] = {Point64} 50,250,0 
          [4] = {Point64} 350,250,0 
          [5] = {Point64} 350,150.5,0 
          [6] = {Point64} 400,150.5,0 
          [7] = {Point64} 400,149.5,0 
           */

        c.Clear();
        PathsD sL2 = new();
        c.AddSubject(outer);
        c.AddClip(kH2);
        c.Execute(ClipType.Difference, FillRule.EvenOdd, sL2);
        
        /* Expected output
        sL2 = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 16
          [0] = {Point64} 349.5,149.5,0 
          [1] = {Point64} 349.5,150,0 
          [2] = {Point64} 250,150,0 
          [3] = {Point64} 250,50,0 
          [4] = {Point64} 150,50,0 
          [5] = {Point64} 150,150,0 
          [6] = {Point64} 50,100,0 
          [7] = {Point64} 50,200,0 
          [8] = {Point64} 150,250,0 
          [9] = {Point64} 150,350,0 
          [10] = {Point64} 250,350,0 
          [11] = {Point64} 250,250,0 
          [12] = {Point64} 350,250,0 
          [13] = {Point64} 350,150.5,0 
          [14] = {Point64} 400,150.5,0 
          [15] = {Point64} 400,149.5,0 
           */

        // Gap removal test
        PathsD gR = GeoWrangler.gapRemoval(kH, 100);
        
        /* Expected output
        gR = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 13
          [0] = {Point64} 400,0,0 
          [1] = {Point64} 400,149.5,0 
          [2] = {Point64} 349.5,149.5,0 
          [3] = {Point64} 349.5,150,0 
          [4] = {Point64} 50,150,0 
          [5] = {Point64} 50,250,0 
          [6] = {Point64} 350,250,0 
          [7] = {Point64} 350,150.5,0 
          [8] = {Point64} 400,150.5,0 
          [9] = {Point64} 400,400,0 
          [10] = {Point64} 0,400,0 
          [11] = {Point64} 0,0,0 
          [12] = {Point64} 400,0,0 
           */
        
        PathsD gR2 = GeoWrangler.gapRemoval(kH2, 100);
        
        /* Expected output
        gR2 = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 21
          [0] = {Point64} 400,0,0 
          [1] = {Point64} 400,149.5,0 
          [2] = {Point64} 349.5,149.5,0 
          [3] = {Point64} 349.5,150,0 
          [4] = {Point64} 250,150,0 
          [5] = {Point64} 250,50,0 
          [6] = {Point64} 150,50,0 
          [7] = {Point64} 150,150,0 
          [8] = {Point64} 50,150,0 
          [9] = {Point64} 50,250,0 
          [10] = {Point64} 150,250,0 
          [11] = {Point64} 150,350,0 
          [12] = {Point64} 250,350,0 
          [13] = {Point64} 250,250,0 
          [14] = {Point64} 350,250,0 
          [15] = {Point64} 350,150.5,0 
          [16] = {Point64} 400,150.5,0 
          [17] = {Point64} 400,400,0 
          [18] = {Point64} 0,400,0 
          [19] = {Point64} 0,0,0 
          [20] = {Point64} 400,0,0 
           */

        PathsD kHSource3 = new();
        kHSource3.AddRange(gR);

        PathsD kHSource4 = new();
        kHSource4.AddRange(gR2);

        // Generate keyholed geometry
        PathsD kH3 = GeoWrangler.makeKeyHole(new PathsD(kHSource3), reverseEval:false, biDirectionalEval:true);
        
        /* Expected output
        kH3 = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 13
          [0] = {Point64} 400,0,0 
          [1] = {Point64} 400,149.5,0 
          [2] = {Point64} 349.5,149.5,0 
          [3] = {Point64} 349.5,150,0 
          [4] = {Point64} 50,150,0 
          [5] = {Point64} 50,250,0 
          [6] = {Point64} 350,250,0 
          [7] = {Point64} 350,150.5,0 
          [8] = {Point64} 400,150.5,0 
          [9] = {Point64} 400,400,0 
          [10] = {Point64} 0,400,0 
          [11] = {Point64} 0,0,0 
          [12] = {Point64} 400,0,0 
           */
        
        PathsD kH4 = GeoWrangler.makeKeyHole(new PathsD(kHSource4), reverseEval:false, biDirectionalEval:true);

        /* Expected output
        kH4 = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 21
          [0] = {Point64} 400,0,0 
          [1] = {Point64} 400,149.5,0 
          [2] = {Point64} 349.5,149.5,0 
          [3] = {Point64} 349.5,150,0 
          [4] = {Point64} 250,150,0 
          [5] = {Point64} 250,50,0 
          [6] = {Point64} 150,50,0 
          [7] = {Point64} 150,150,0 
          [8] = {Point64} 50,150,0 
          [9] = {Point64} 50,250,0 
          [10] = {Point64} 150,250,0 
          [11] = {Point64} 150,350,0 
          [12] = {Point64} 250,350,0 
          [13] = {Point64} 250,250,0 
          [14] = {Point64} 350,250,0 
          [15] = {Point64} 350,150.5,0 
          [16] = {Point64} 400,150.5,0 
          [17] = {Point64} 400,400,0 
          [18] = {Point64} 0,400,0 
          [19] = {Point64} 0,0,0 
          [20] = {Point64} 400,0,0 
           */
        
        // Sliver removal test
        PathsD sR = GeoWrangler.gapRemoval(sL, -100);
        
        /* Expected output
        sR = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 9
          [0] = {Point64} 50,150,0 
          [1] = {Point64} 50,250,0 
          [2] = {Point64} 350,250,0 
          [3] = {Point64} 350,150.5,0 
          [4] = {Point64} 400,150.5,0 
          [5] = {Point64} 400,149.5,0 
          [6] = {Point64} 349.5,149.5,0 
          [7] = {Point64} 349.5,150,0 
          [8] = {Point64} 50,150,0 
           */
        
        PathsD sR2 = GeoWrangler.gapRemoval(sL2, -100);
        
        /* Expected output
        sR2 = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 17
          [0] = {Point64} 50,150,0 
          [1] = {Point64} 50,250,0 
          [2] = {Point64} 150,250,0 
          [3] = {Point64} 150,350,0 
          [4] = {Point64} 250,350,0 
          [5] = {Point64} 250,250,0 
          [6] = {Point64} 350,250,0 
          [7] = {Point64} 350,150.5,0 
          [8] = {Point64} 400,150.5,0 
          [9] = {Point64} 400,149.5,0 
          [10] = {Point64} 349.5,149.5,0 
          [11] = {Point64} 349.5,150,0 
          [12] = {Point64} 250,150,0 
          [13] = {Point64} 250,50,0 
          [14] = {Point64} 150,50,0 
          [15] = {Point64} 150,150,0 
          [16] = {Point64} 50,150,0 
           */
    }

    private static void selfOverlapTest()
    {
        PathD outer = new()
        {
            new(0, 0),
            new(110, 0),
            new(110, 50),
            new(50, 50),
            new(50, 150),
            new(150, 150),
            new(150, 50),
            new(90, 50),
            new(90, 0),
            new(200, 0),
            new(200, 200),
            new(0, 200),
            new(0, 0)
        };

        // decomposer test
        PathsD dSource = new() {outer};
        PathsD decomp = GeoWrangler.decompose(dSource);
        
        /* Expected output
        decomp = {List<List<Point64>>} Count = 2
         [0] = {List<Point64>} Count = 6
          [0] = {Point64} 0,0,0 
          [1] = {Point64} 0,200,0 
          [2] = {Point64} 200,200,0 
          [3] = {Point64} 200,0,0 
          [4] = {Point64} 110,0,0 
          [5] = {Point64} 0,0,0 
         [1] = {List<Point64>} Count = 6
          [0] = {Point64} 50,50,0 
          [1] = {Point64} 90,50,0 
          [2] = {Point64} 150,50,0 
          [3] = {Point64} 150,150,0 
          [4] = {Point64} 50,150,0 
          [5] = {Point64} 50,50,0 
           */
        
        PathsD kHD = GeoWrangler.makeKeyHole(dSource, reverseEval:false, biDirectionalEval:true);

        /* Expected output
        kHD = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 13
          [0] = {Point64} 200,0,0 
          [1] = {Point64} 200,49.5,0 
          [2] = {Point64} 149.5,49.5,0 
          [3] = {Point64} 149.5,50,0 
          [4] = {Point64} 50,50,0 
          [5] = {Point64} 50,150,0 
          [6] = {Point64} 150,150,0 
          [7] = {Point64} 150,50.5,0 
          [8] = {Point64} 200,50.5,0 
          [9] = {Point64} 200,200,0 
          [10] = {Point64} 0,200,0 
          [11] = {Point64} 0,0,0 
          [12] = {Point64} 200,0,0 
           */
        
        // keyholer test
        PathsD kHSource = new() {outer};
        PathsD kH = GeoWrangler.makeKeyHole(kHSource, reverseEval:false, biDirectionalEval:true);

        /* Expected output
        kH = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 13
          [0] = {Point64} 200,0,0 
          [1] = {Point64} 200,49.5,0 
          [2] = {Point64} 149.5,49.5,0 
          [3] = {Point64} 149.5,50,0 
          [4] = {Point64} 50,50,0 
          [5] = {Point64} 50,15,0 
          [6] = {Point64} 15,15,0 
          [7] = {Point64} 15,50.5,0 
          [8] = {Point64} 200,50.5,0 
          [9] = {Point64} 200,200,0 
          [10] = {Point64} 0,200,0 
          [11] = {Point64} 0,0,0 
          [12] = {Point64} 200,0,0 
           */
        
        ClipperD c = new(Constants.roundingDecimalPrecision);
        c.AddSubject(outer);

        PathsD unionRes = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, unionRes);
        
        /* Expected output
        unionRes = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 13
          [0] = {Point64} 110,50,0 
          [1] = {Point64} 150,50,0 
          [2] = {Point64} 150,150,0 
          [3] = {Point64} 50,150,0 
          [4] = {Point64} 50,50,0 
          [5] = {Point64} 90,50,0 
          [6] = {Point64} 90,0,0 
          [7] = {Point64} 0,0,0 
          [8] = {Point64} 0,200,0 
          [9] = {Point64} 200,200,0 
          [10] = {Point64} 200,0,0 
          [11] = {Point64} 110,0,0 
          [12] = {Point64} 110,50,0 
           */
        
        PathsD unionRes_kH = GeoWrangler.makeKeyHole(unionRes, reverseEval:false, biDirectionalEval:true);
        
        /* Expected output
        unionRes_kH = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 13
          [0] = {Point64} 0,0,0 
          [1] = {Point64} 0,200,0 
          [2] = {Point64} 200,200,0 
          [3] = {Point64} 200,0,0 
          [4] = {Point64} 110,0,0 
          [5] = {Point64} 110,50,0 
          [6] = {Point64} 150,50,0 
          [7] = {Point64} 150,150,0 
          [8] = {Point64} 50,150,0 
          [9] = {Point64} 50,50,0 
          [10] = {Point64} 90,50,0 
          [11] = {Point64} 90,0,0 
          [12] = {Point64} 0,0,0 
           */
        
        PathsD unionResc = GeoWrangler.close(unionRes);
        PathsD unionResc_kH = GeoWrangler.makeKeyHole(unionResc, reverseEval:false, biDirectionalEval:true);

        /* Expected output
        unionResc_kH = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 13
          [0] = {Point64} 0,0,0 
          [1] = {Point64} 0,200,0 
          [2] = {Point64} 200,200,0 
          [3] = {Point64} 200,0,0 
          [4] = {Point64} 110,0,0 
          [5] = {Point64} 110,50,0 
          [6] = {Point64} 150,50,0 
          [7] = {Point64} 150,150,0 
          [8] = {Point64} 50,150,0 
          [9] = {Point64} 50,50,0 
          [10] = {Point64} 90,50,0 
          [11] = {Point64} 90,0,0 
          [12] = {Point64} 0,0,0 
           */
        
        PathsD unionResP = new();
        c.Execute(ClipType.Union, FillRule.Positive, unionResP);
        
        /* Expected output
         No geometry in unionResP
         */
        
        // no keyhole for any of the below
        PathsD unionResP_kH = GeoWrangler.makeKeyHole(unionResP, reverseEval:false, biDirectionalEval:true);
        PathsD unionResPc = GeoWrangler.close(unionResP);
        PathsD unionResPc_kH = GeoWrangler.makeKeyHole(unionResPc, reverseEval:false, biDirectionalEval:true);

        // seems good - get keyhole
        PathsD unionResNZ = new();
        c.Execute(ClipType.Union, FillRule.NonZero, unionResNZ);
        
        /* Expected output
        unionResNZ = {List<List<Point64>>} Count = 2
         [0] = {List<Point64>} Count = 6
          [0] = {Point64} 0,0,0 
          [1] = {Point64} 0,200,0 
          [2] = {Point64} 200,200,0 
          [3] = {Point64} 200,0,0 
          [4] = {Point64} 110,0,0 
          [5] = {Point64} 0,0,0 
         [1] = {List<Point64>} Count = 6
          [0] = {Point64} 90,50,0 
          [1] = {Point64} 150,50,0 
          [2] = {Point64} 150,150,0 
          [3] = {Point64} 50,150,0 
          [4] = {Point64} 50,50,0 
          [5] = {Point64} 90,50,0 
           */
        
        PathsD unionResNZ_kH = GeoWrangler.makeKeyHole(unionResNZ, reverseEval:false, biDirectionalEval:true);
        
        /* Expected output
        unionResNZ_kH = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 13
          [0] = {Point64} 200,0,0 
          [1] = {Point64} 200,49.5,0 
          [2] = {Point64} 149.5,49.5,0 
          [3] = {Point64} 149.5,50,0 
          [4] = {Point64} 50,50,0 
          [5] = {Point64} 50,150,0 
          [6] = {Point64} 150,150,0 
          [7] = {Point64} 150,50.5,0 
          [8] = {Point64} 200,50.5,0 
          [9] = {Point64} 200,200,0 
          [10] = {Point64} 0,200,0 
          [11] = {Point64} 0,0,0 
          [12] = {Point64} 200,0,0 
           */
        
        PathsD unionResNZc = GeoWrangler.close(unionResNZ);
        
        /* Expected output
        unionResNZc = {List<List<Point64>>} Count = 2
         [0] = {List<Point64>} Count = 6
          [0] = {Point64} 0,0,0 
          [1] = {Point64} 0,200,0 
          [2] = {Point64} 200,200,0 
          [3] = {Point64} 200,0,0 
          [4] = {Point64} 110,0,0 
          [5] = {Point64} 0,0,0 
         [1] = {List<Point64>} Count = 6
          [0] = {Point64} 90,50,0 
          [1] = {Point64} 150,50,0 
          [2] = {Point64} 150,150,0 
          [3] = {Point64} 50,150,0 
          [4] = {Point64} 50,50,0 
          [5] = {Point64} 90,50,0 
           */
        
        PathsD unionResNZc_kH = GeoWrangler.makeKeyHole(unionResNZc, reverseEval:false, biDirectionalEval:true);
        
        /* Expected result
        unionResNZc_kH = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 13
          [0] = {Point64} 200,0,0 
          [1] = {Point64} 200,49.5,0 
          [2] = {Point64} 149.5,49.5,0 
          [3] = {Point64} 149.5,50,0 
          [4] = {Point64} 50,50,0 
          [5] = {Point64} 50,150,0 
          [6] = {Point64} 150,150,0 
          [7] = {Point64} 150,50.5,0 
          [8] = {Point64} 200,50.5,0 
          [9] = {Point64} 200,200,0 
          [10] = {Point64} 0,200,0 
          [11] = {Point64} 0,0,0 
          [12] = {Point64} 200,0,0 
           */

        PathsD simplifyRes = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, simplifyRes);
        simplifyRes = GeoWrangler.stripCollinear(simplifyRes);

        /* Expected output
        simplifyRes = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 13
          [0] = {Point64} 110,50,0 
          [1] = {Point64} 150,50,0 
          [2] = {Point64} 150,150,0 
          [3] = {Point64} 50,150,0 
          [4] = {Point64} 50,50,0 
          [5] = {Point64} 90,50,0 
          [6] = {Point64} 90,0,0 
          [7] = {Point64} 0,0,0 
          [8] = {Point64} 0,200,0 
          [9] = {Point64} 200,200,0 
          [10] = {Point64} 200,0,0 
          [11] = {Point64} 110,0,0 
          [12] = {Point64} 110,50,0 
           */

        PathsD simplifyRes_kH = GeoWrangler.makeKeyHole(simplifyRes, reverseEval:false, biDirectionalEval:true);
        
        /* Expected output
        simplifyRes_kH = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 13
          [0] = {Point64} 0,0,0 
          [1] = {Point64} 0,200,0 
          [2] = {Point64} 200,200,0 
          [3] = {Point64} 200,0,0 
          [4] = {Point64} 110,0,0 
          [5] = {Point64} 110,50,0 
          [6] = {Point64} 150,50,0 
          [7] = {Point64} 150,150,0 
          [8] = {Point64} 50,150,0 
          [9] = {Point64} 50,50,0 
          [10] = {Point64} 90,50,0 
          [11] = {Point64} 90,0,0 
          [12] = {Point64} 0,0,0 
           */
        
        PathsD simplifyResc = GeoWrangler.close(simplifyRes);
        
        /* Expected output
        simplifyResc = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 13
          [0] = {Point64} 110,50,0 
          [1] = {Point64} 150,50,0 
          [2] = {Point64} 150,150,0 
          [3] = {Point64} 50,150,0 
          [4] = {Point64} 50,50,0 
          [5] = {Point64} 90,50,0 
          [6] = {Point64} 90,0,0 
          [7] = {Point64} 0,0,0 
          [8] = {Point64} 0,200,0 
          [9] = {Point64} 200,200,0 
          [10] = {Point64} 200,0,0 
          [11] = {Point64} 110,0,0 
          [12] = {Point64} 110,50,0 
           */
        
        PathsD simplifyResc_kH = GeoWrangler.makeKeyHole(simplifyResc, reverseEval:false, biDirectionalEval:true);

        /* Expected output
        simplifyResc_kH = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 13
          [0] = {Point64} 0,0,0 
          [1] = {Point64} 0,200,0 
          [2] = {Point64} 200,200,0 
          [3] = {Point64} 200,0,0 
          [4] = {Point64} 110,0,0 
          [5] = {Point64} 110,50,0 
          [6] = {Point64} 150,50,0 
          [7] = {Point64} 150,150,0 
          [8] = {Point64} 50,150,0 
          [9] = {Point64} 50,50,0 
          [10] = {Point64} 90,50,0 
          [11] = {Point64} 90,0,0 
          [12] = {Point64} 0,0,0 
           */
        
        PathsD simplifyRes2 = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, simplifyRes2);
        PathsD simplifyRes2_kH = GeoWrangler.makeKeyHole(simplifyRes2, reverseEval:false, biDirectionalEval:true);
        PathsD simplifyRes2c = GeoWrangler.close(simplifyRes2);
        PathsD simplifyRes2c_kH = GeoWrangler.makeKeyHole(simplifyRes2c, reverseEval:false, biDirectionalEval:true);

        // no good - no result
        PathsD intRes = new();
        c.Execute(ClipType.Intersection, FillRule.EvenOdd, intRes);
        PathsD intRes_kH = GeoWrangler.makeKeyHole(intRes, reverseEval:false, biDirectionalEval:true);
        PathsD intResc = GeoWrangler.close(intRes);
        PathsD intResc_kH = GeoWrangler.makeKeyHole(intResc, reverseEval:false, biDirectionalEval:true);

        // no good - no result
        PathsD intResP = new();
        c.Execute(ClipType.Intersection, FillRule.Positive, intResP);
        PathsD intResP_kH = GeoWrangler.makeKeyHole(intResP, reverseEval:false, biDirectionalEval:true);
        PathsD intResPc = GeoWrangler.close(intResP);
        PathsD intResPc_kH = GeoWrangler.makeKeyHole(intResPc, reverseEval:false, biDirectionalEval:true);

        // no good - no result
        PathsD intResNZ = new();
        c.Execute(ClipType.Intersection, FillRule.NonZero, intResNZ);
        PathsD intResNZ_kH = GeoWrangler.makeKeyHole(intResNZ, reverseEval:false, biDirectionalEval:true);
        PathsD intResNZc = GeoWrangler.close(intResNZ);
        PathsD intResNZc_kH = GeoWrangler.makeKeyHole(intResNZc, reverseEval:false, biDirectionalEval:true);

        RectD bounds = Clipper.GetBounds(new PathsD { outer });
        PathD bb = new()
        {
            new(bounds.left, bounds.bottom),
            new(bounds.left, bounds.top),
            new(bounds.right, bounds.top),
            new(bounds.right, bounds.bottom)
        };

        c.Clear();
        c.AddSubject(bb);
        c.AddClip(outer);

        PathsD intRes2 = new();
        c.Execute(ClipType.Intersection, FillRule.EvenOdd, intRes2);
        
        /* Expected output
         intRes2 = {List<List<Point64>>} Count = 2
          [0] = {List<Point64>} Count = 6
           [0] = {Point64} 0,0,0 
           [1] = {Point64} 0,200,0 
           [2] = {Point64} 200,200,0 
           [3] = {Point64} 200,0,0 
           [4] = {Point64} 110,0,0 
           [5] = {Point64} 0,0,0 
          [1] = {List<Point64>} Count = 9
           [0] = {Point64} 90,0,0 
           [1] = {Point64} 110,0,0 
           [2] = {Point64} 110,50,0 
           [3] = {Point64} 150,50,0 
           [4] = {Point64} 150,150,0 
           [5] = {Point64} 50,150,0 
           [6] = {Point64} 50,50,0 
           [7] = {Point64} 90,50,0 
           [8] = {Point64} 90,0,0 
           */
        
        PathsD intRes2_kH = GeoWrangler.makeKeyHole(intRes2, reverseEval:false, biDirectionalEval:true);
        
        /* Expected output
         intRes2_kH = {List<List<Point64>>} Count = 1
          [0] = {List<Point64>} Count = 13
           [0] = {Point64} 0,0,0 
           [1] = {Point64} 0,200,0 
           [2] = {Point64} 200,200,0 
           [3] = {Point64} 200,0,0 
           [4] = {Point64} 110,0,0 
           [5] = {Point64} 110,50,0 
           [6] = {Point64} 150,50,0 
           [7] = {Point64} 150,150,0 
           [8] = {Point64} 50,150,0 
           [9] = {Point64} 50,50,0 
           [10] = {Point64} 90,50,0 
           [11] = {Point64} 90,0,0 
           [12] = {Point64} 0,0,0 
           */
        
        PathsD intRes2c = GeoWrangler.close(intRes2);
        
        /* Expected output
         intRes2c = {List<List<Point64>>} Count = 2
          [0] = {List<Point64>} Count = 6
           [0] = {Point64} 0,0,0 
           [1] = {Point64} 0,200,0 
           [2] = {Point64} 200,200,0 
           [3] = {Point64} 200,0,0 
           [4] = {Point64} 110,0,0 
           [5] = {Point64} 0,0,0 
          [1] = {List<Point64>} Count = 9
           [0] = {Point64} 90,0,0 
           [1] = {Point64} 110,0,0 
           [2] = {Point64} 110,50,0 
           [3] = {Point64} 150,50,0 
           [4] = {Point64} 150,150,0 
           [5] = {Point64} 50,150,0 
           [6] = {Point64} 50,50,0 
           [7] = {Point64} 90,50,0 
           [8] = {Point64} 90,0,0 
           */
        
        PathsD intRes2c_kH = GeoWrangler.makeKeyHole(intRes2c, reverseEval:false, biDirectionalEval:true);

        /* Expected output
         intRes2c_kH = {List<List<Point64>>} Count = 1
          [0] = {List<Point64>} Count = 13
           [0] = {Point64} 0,0,0 
           [1] = {Point64} 0,200,0 
           [2] = {Point64} 200,200,0 
           [3] = {Point64} 200,0,0 
           [4] = {Point64} 110,0,0 
           [5] = {Point64} 110,50,0 
           [6] = {Point64} 150,50,0 
           [7] = {Point64} 150,150,0 
           [8] = {Point64} 50,150,0 
           [9] = {Point64} 50,50,0 
           [10] = {Point64} 90,50,0 
           [11] = {Point64} 90,0,0 
           [12] = {Point64} 0,0,0 
           */
        
        PathsD intRes2P = new();
        c.Execute(ClipType.Intersection, FillRule.Positive, intRes2P);
        
        /* Expected output
         No geometry in intRes2P
         */
        
        // No keyholes as no geometry.
        PathsD intRes2P_kH = GeoWrangler.makeKeyHole(intRes2P, reverseEval:false, biDirectionalEval:true);
        PathsD intRes2Pc = GeoWrangler.close(intRes2P);
        PathsD intRes2Pc_kH = GeoWrangler.makeKeyHole(intRes2Pc, reverseEval:false, biDirectionalEval:true);

        // seems good - get keyhole
        PathsD intRes2NZ = new();
        c.Execute(ClipType.Intersection, FillRule.NonZero, intRes2NZ);
        
        /* Expected output
         intRes2NZ = {List<List<Point64>>} Count = 2
          [0] = {List<Point64>} Count = 6
           [0] = {Point64} 0,0,0 
           [1] = {Point64} 0,200,0 
           [2] = {Point64} 200,200,0 
           [3] = {Point64} 200,0,0 
           [4] = {Point64} 110,0,0 
           [5] = {Point64} 0,0,0 
          [1] = {List<Point64>} Count = 6
           [0] = {Point64} 90,50,0 
           [1] = {Point64} 150,50,0 
           [2] = {Point64} 150,150,0 
           [3] = {Point64} 50,150,0 
           [4] = {Point64} 50,50,0 
           [5] = {Point64} 90,50,0 
           */
        
        PathsD intRes2NZ_kH = GeoWrangler.makeKeyHole(intRes2NZ, reverseEval:false, biDirectionalEval:true);
        
        /* Expected output
         intRes2NZ_kH = {List<List<Point64>>} Count = 1
          [0] = {List<Point64>} Count = 13
           [0] = {Point64} 200,0,0 
           [1] = {Point64} 200,49.5,0 
           [2] = {Point64} 149.5,49.5,0 
           [3] = {Point64} 149.5,50,0 
           [4] = {Point64} 50,50,0 
           [5] = {Point64} 50,150,0 
           [6] = {Point64} 150,150,0 
           [7] = {Point64} 150,50.5,0 
           [8] = {Point64} 200,50.5,0 
           [9] = {Point64} 200,200,0 
           [10] = {Point64} 0,200,0 
           [11] = {Point64} 0,0,0 
           [12] = {Point64} 200,0,0 
           */
        
        PathsD intRes2NZc = GeoWrangler.close(intRes2NZ);
        
        /* Expected output
         intRes2NZc = {List<List<Point64>>} Count = 2
          [0] = {List<Point64>} Count = 6
           [0] = {Point64} 0,0,0 
           [1] = {Point64} 0,200,0 
           [2] = {Point64} 200,200,0 
           [3] = {Point64} 200,0,0 
           [4] = {Point64} 110,0,0 
           [5] = {Point64} 0,0,0 
          [1] = {List<Point64>} Count = 6
           [0] = {Point64} 90,50,0 
           [1] = {Point64} 150,50,0 
           [2] = {Point64} 150,150,0 
           [3] = {Point64} 50,150,0 
           [4] = {Point64} 50,50,0 
           [5] = {Point64} 90,50,0 
           */
        
        PathsD intRes2NZc_kH = GeoWrangler.makeKeyHole(intRes2NZc, reverseEval:false, biDirectionalEval:true);
        
         /* Expected output
          intRes2NZc_kH = {List<List<Point64>>} Count = 1
           [0] = {List<Point64>} Count = 13
            [0] = {Point64} 200,0,0 
            [1] = {Point64} 200,49.5,0 
            [2] = {Point64} 149.5,49.5,0 
            [3] = {Point64} 149.5,50,0 
            [4] = {Point64} 50,50,0 
            [5] = {Point64} 50,150,0 
            [6] = {Point64} 150,150,0 
            [7] = {Point64} 150,50.5,0 
            [8] = {Point64} 200,50.5,0 
            [9] = {Point64} 200,200,0 
            [10] = {Point64} 0,200,0 
            [11] = {Point64} 0,0,0 
            [12] = {Point64} 200,0,0 
          */
        
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

        // decomposer test
        PathsD dSource = new() {outer};

        /* Expected
         dSource = {List<List<Point64>>} Count = 1
          [0] = {List<Point64>} Count = 13
           [0] = {Point64} 0,0,0 
           [1] = {Point64} 0,200,0 
           [2] = {Point64} 200,200,0 
           [3] = {Point64} 200,0,0 
           [4] = {Point64} 90,0,0 
           [5] = {Point64} 90,50,0 
           [6] = {Point64} 150,50,0 
           [7] = {Point64} 150,150,0 
           [8] = {Point64} 50,150,0 
           [9] = {Point64} 50,50,0 
           [10] = {Point64} 110,50,0 
           [11] = {Point64} 110,0,0 
           [12] = {Point64} 0,0,0 
           */
        
        PathsD decomp = GeoWrangler.decompose(dSource);
        
        /* Expected output
         decomp = {List<List<Point64>>} Count = 2
          [0] = {List<Point64>} Count = 6
           [0] = {Point64} 0,0,0 
           [1] = {Point64} 0,200,0 
           [2] = {Point64} 200,200,0 
           [3] = {Point64} 200,0,0 
           [4] = {Point64} 110,0,0 
           [5] = {Point64} 0,0,0 
          [1] = {List<Point64>} Count = 6
           [0] = {Point64} 50,50,0 
           [1] = {Point64} 90,50,0 
           [2] = {Point64} 150,50,0 
           [3] = {Point64} 150,150,0 
           [4] = {Point64} 50,150,0 
           [5] = {Point64} 50,50,0 
           */
        
        PathsD kHD = GeoWrangler.makeKeyHole(dSource, reverseEval:false, biDirectionalEval:true);
        
        /* Expected output
         kHD = {List<List<Point64>>} Count = 1
          [0] = {List<Point64>} Count = 13
           [0] = {Point64} 200,0,0 
           [1] = {Point64} 200,49.5,0 
           [2] = {Point64} 149.5,49.5,0 
           [3] = {Point64} 149.5,50,0 
           [4] = {Point64} 50,50.,0 
           [5] = {Point64} 50,150,0 
           [6] = {Point64} 150,150,0 
           [7] = {Point64} 150,50.5,0 
           [8] = {Point64} 200,50.5,0 
           [9] = {Point64} 200,200,0 
           [10] = {Point64} 0,200,0 
           [11] = {Point64} 0,0,0 
           [12] = {Point64} 200,0,0 
           */

        // keyholer test
        PathsD kHSource = new() {outer};
        
        /* Expected
         kHSource = {List<List<Point64>>} Count = 1
          [0] = {List<Point64>} Count = 13
           [0] = {Point64} 0,0,0 
           [1] = {Point64} 0,200,0 
           [2] = {Point64} 200,200,0 
           [3] = {Point64} 200,0,0 
           [4] = {Point64} 90,0,0 
           [5] = {Point64} 90,50,0 
           [6] = {Point64} 150,50,0 
           [7] = {Point64} 150,150,0 
           [8] = {Point64} 50,150,0 
           [9] = {Point64} 50,50,0 
           [10] = {Point64} 110,50,0 
           [11] = {Point64} 110,0,0 
           [12] = {Point64} 0,0,0 
           */
        
        PathsD kH = GeoWrangler.makeKeyHole(kHSource, reverseEval:false, biDirectionalEval:true);

        /* Expected output
         kH = {List<List<Point64>>} Count = 1
          [0] = {List<Point64>} Count = 13
           [0] = {Point64} 200,0,0 
           [1] = {Point64} 200,49.5,0 
           [2] = {Point64} 149.5,49.5,0 
           [3] = {Point64} 149.5,50,0 
           [4] = {Point64} 50,5,0 
           [5] = {Point64} 50,150,0 
           [6] = {Point64} 150,150,0 
           [7] = {Point64} 150,50.5,0 
           [8] = {Point64} 200,50.5,0 
           [9] = {Point64} 200,200,0 
           [10] = {Point64} 0,200,0 
           [11] = {Point64} 0,0,0 
           [12] = {Point64} 200,0,0 
           */
        
        ClipperD c = new(Constants.roundingDecimalPrecision);
        c.AddSubject(outer);

        PathsD unionRes = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, unionRes);

        /* Expected output
         unionRes = {List<List<Point64>>} Count = 1
          [0] = {List<Point64>} Count = 13
           [0] = {Point64} 90000,0,0 
           [1] = {Point64} 0,0,0 
           [2] = {Point64} 0,200,0 
           [3] = {Point64} 200,200,0 
           [4] = {Point64} 200,0,0 
           [5] = {Point64} 110,0,0 
           [6] = {Point64} 110,50,0 
           [7] = {Point64} 150,50,0 
           [8] = {Point64} 150,150,0 
           [9] = {Point64} 50,150,0 
           [10] = {Point64} 50,50,0 
           [11] = {Point64} 90,50,0 
           [12] = {Point64} 90,0,0 
           */
        
        PathsD unionRes_kH = GeoWrangler.makeKeyHole(unionRes, reverseEval:false, biDirectionalEval:true);
        
        /* Expected output
         unionRes_kH = {List<List<Point64>>} Count = 1
          [0] = {List<Point64>} Count = 13
           [0] = {Point64} 0,0,0 
           [1] = {Point64} 0,200,0 
           [2] = {Point64} 200,200,0 
           [3] = {Point64} 200,0,0 
           [4] = {Point64} 110,0,0 
           [5] = {Point64} 110,50,0 
           [6] = {Point64} 150,50,0 
           [7] = {Point64} 150,150,0 
           [8] = {Point64} 50,150,0 
           [9] = {Point64} 50,50,0 
           [10] = {Point64} 90,50,0 
           [11] = {Point64} 90,0,0 
           [12] = {Point64} 0,0,0 
           */
        
        PathsD unionResc = GeoWrangler.close(unionRes);
        
        /* Expected output
         unionResc = {List<List<Point64>>} Count = 1
          [0] = {List<Point64>} Count = 13
           [0] = {Point64} 90,0,0 
           [1] = {Point64} 0,0,0 
           [2] = {Point64} 0,200,0 
           [3] = {Point64} 200,200,0 
           [4] = {Point64} 200,0,0 
           [5] = {Point64} 110,0,0 
           [6] = {Point64} 110,50,0 
           [7] = {Point64} 150,50,0 
           [8] = {Point64} 150,150,0 
           [9] = {Point64} 50,150,0 
           [10] = {Point64} 50,50,0 
           [11] = {Point64} 90,50,0 
           [12] = {Point64} 90,0,0 
           */
        
        PathsD unionResc_kH = GeoWrangler.makeKeyHole(unionResc, reverseEval:false, biDirectionalEval:true);

        /* Expected output
         unionResc_kH = {List<List<Point64>>} Count = 1
          [0] = {List<Point64>} Count = 13
           [0] = {Point64} 0,0,0 
           [1] = {Point64} 0,200,0 
           [2] = {Point64} 200,200,0 
           [3] = {Point64} 200,0,0 
           [4] = {Point64} 110,0,0 
           [5] = {Point64} 110,50,0 
           [6] = {Point64} 150,50,0 
           [7] = {Point64} 150,150,0 
           [8] = {Point64} 50,150,0 
           [9] = {Point64} 50,50,0 
           [10] = {Point64} 90,50,0 
           [11] = {Point64} 90,0,0 
           [12] = {Point64} 0,0,0 
           */
        
        PathsD unionResP = new();
        c.Execute(ClipType.Union, FillRule.Positive, unionResP);

        /* Expected output
           unionResP = {List<List<Point64>>} Count = 2
            [0] = {List<Point64>} Count = 6
             [0] = {Point64} 0,0,0 
             [1] = {Point64} 0,200,0 
             [2] = {Point64} 200,200,0 
             [3] = {Point64} 200,0,0 
             [4] = {Point64} 110,0,0 
             [5] = {Point64} 0,0,0 
            [1] = {List<Point64>} Count = 6
             [0] = {Point64} 90,50,0 
             [1] = {Point64} 150,50,0 
             [2] = {Point64} 150,150,0 
             [3] = {Point64} 50,150,0 
             [4] = {Point64} 50,50,0 
             [5] = {Point64} 90,50,0 
           */
        
        PathsD unionResP_kH = GeoWrangler.makeKeyHole(unionResP, reverseEval:false, biDirectionalEval:true);
        
        /* Expected output
         unionResP_kH = {List<List<Point64>>} Count = 1
          [0] = {List<Point64>} Count = 13
           [0] = {Point64} 200,0,0 
           [1] = {Point64} 200,49.5,0 
           [2] = {Point64} 149.5,49.5,0 
           [3] = {Point64} 149.5,50,0 
           [4] = {Point64} 50,50,0 
           [5] = {Point64} 50,150,0 
           [6] = {Point64} 150,150,0 
           [7] = {Point64} 150,50.5,0 
           [8] = {Point64} 200,50.5,0 
           [9] = {Point64} 200,200,0 
           [10] = {Point64} 0,200,0 
           [11] = {Point64} 0,0,0 
           [12] = {Point64} 200,0,0 
         */
        
        PathsD unionResPc = GeoWrangler.close(unionResP);
        
        /* Expected output
          unionResPc = {List<List<Point64>>} Count = 2
           [0] = {List<Point64>} Count = 6
            [0] = {Point64} 0,0,0 
            [1] = {Point64} 0,200,0 
            [2] = {Point64} 200,200,0 
            [3] = {Point64} 200,0,0 
            [4] = {Point64} 110,0,0 
            [5] = {Point64} 0,0,0 
           [1] = {List<Point64>} Count = 6
            [0] = {Point64} 90,50,0 
            [1] = {Point64} 150,50,0 
            [2] = {Point64} 150,150,0 
            [3] = {Point64} 50,150,0 
            [4] = {Point64} 50,50,0 
            [5] = {Point64} 90,50,0 
           */
        
        PathsD unionResPc_kH = GeoWrangler.makeKeyHole(unionResPc, reverseEval:false, biDirectionalEval:true);

        /* Expected output
         unionResPc_kH = {List<List<Point64>>} Count = 1
          [0] = {List<Point64>} Count = 13
           [0] = {Point64} 200,0,0 
           [1] = {Point64} 200,49.5,0 
           [2] = {Point64} 149.5,49.500,0 
           [3] = {Point64} 149.5,50,0 
           [4] = {Point64} 50,50,0 
           [5] = {Point64} 50,150,0 
           [6] = {Point64} 150,150,0 
           [7] = {Point64} 150,50.5,0 
           [8] = {Point64} 200,50.5,0 
           [9] = {Point64} 200,200,0 
           [10] = {Point64} 0,200,0 
           [11] = {Point64} 0,0,0 
           [12] = {Point64} 200,0,0 
           */
        
        // seems good - get keyhole
        PathsD unionResNZ = new();
        c.Execute(ClipType.Union, FillRule.NonZero, unionResNZ);
        
        /* Expected output
           unionResNZ = {List<List<Point64>>} Count = 2
            [0] = {List<Point64>} Count = 6
             [0] = {Point64} 0,0,0 
             [1] = {Point64} 0,200,0 
             [2] = {Point64} 200,200,0 
             [3] = {Point64} 200,0,0 
             [4] = {Point64} 110,0,0 
             [5] = {Point64} 0,0,0 
            [1] = {List<Point64>} Count = 6
             [0] = {Point64} 90,50,0 
             [1] = {Point64} 150,50,0 
             [2] = {Point64} 150,150,0 
             [3] = {Point64} 50,150,0 
             [4] = {Point64} 50,50,0 
             [5] = {Point64} 90,50,0 
           */
        
        PathsD unionResNZ_kH = GeoWrangler.makeKeyHole(unionResNZ, reverseEval:false, biDirectionalEval:true);
        
        /* Expected output
           unionResNZ_kH = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 13
             [0] = {Point64} 200,0,0 
             [1] = {Point64} 200,49.5,0 
             [2] = {Point64} 149.5,49.5,0 
             [3] = {Point64} 149.5,50,0 
             [4] = {Point64} 50,50,0 
             [5] = {Point64} 50,150,0 
             [6] = {Point64} 150,150,0 
             [7] = {Point64} 150,50.5,0 
             [8] = {Point64} 200,50.5,0 
             [9] = {Point64} 200,200,0 
             [10] = {Point64} 0,200,0 
             [11] = {Point64} 0,0,0 
             [12] = {Point64} 200,0,0 
           */
        
        PathsD unionResNZc = GeoWrangler.close(unionResNZ);
        
        /* Expected output
           unionResNZc = {List<List<Point64>>} Count = 2
            [0] = {List<Point64>} Count = 6
             [0] = {Point64} 0,0,0 
             [1] = {Point64} 0,200,0 
             [2] = {Point64} 200,200,0 
             [3] = {Point64} 200,0,0 
             [4] = {Point64} 110,0,0 
             [5] = {Point64} 0,0,0 
            [1] = {List<Point64>} Count = 6
             [0] = {Point64} 90,50,0 
             [1] = {Point64} 150,50,0 
             [2] = {Point64} 150,150,0 
             [3] = {Point64} 50,150,0 
             [4] = {Point64} 50,50,0 
             [5] = {Point64} 90,50,0 
           */
        
        PathsD unionResNZc_kH = GeoWrangler.makeKeyHole(unionResNZc, reverseEval:false, biDirectionalEval:true);

        /* Expected output
           unionResNZc_kH = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 13
             [0] = {Point64} 200,0,0 
             [1] = {Point64} 200,49.5,0 
             [2] = {Point64} 149.5,49.5,0 
             [3] = {Point64} 149.5,50,0 
             [4] = {Point64} 50,50,0 
             [5] = {Point64} 50,150,0 
             [6] = {Point64} 150,150,0 
             [7] = {Point64} 150,50.5,0 
             [8] = {Point64} 200,50.5,0 
             [9] = {Point64} 200,200,0 
             [10] = {Point64} 0,200,0 
             [11] = {Point64} 0,0,0 
             [12] = {Point64} 200,0,0 
           */
        
        // no good - overlap region is a gap.
        PathsD simplifyRes = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, simplifyRes);
        simplifyRes = GeoWrangler.stripCollinear(simplifyRes);
        
        /* Expected output
           simplifyRes = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 13
             [0] = {Point64} 90,0,0 
             [1] = {Point64} 0,0,0 
             [2] = {Point64} 0,200,0 
             [3] = {Point64} 200,200,0 
             [4] = {Point64} 200,0,0 
             [5] = {Point64} 110,0,0 
             [6] = {Point64} 110,50,0 
             [7] = {Point64} 150,50,0 
             [8] = {Point64} 150,150,0 
             [9] = {Point64} 50,150,0 
             [10] = {Point64} 50,50,0 
             [11] = {Point64} 90,50,0 
             [12] = {Point64} 90,0,0 
           */
        
        PathsD simplifyRes_kH = GeoWrangler.makeKeyHole(simplifyRes, reverseEval:false, biDirectionalEval:true);
        
        /* Expected output
           simplifyRes_kH = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 13
             [0] = {Point64} 0,0,0 
             [1] = {Point64} 0,200,0 
             [2] = {Point64} 200,200,0 
             [3] = {Point64} 200,0,0 
             [4] = {Point64} 110,0,0 
             [5] = {Point64} 110,50,0 
             [6] = {Point64} 150,50,0 
             [7] = {Point64} 150,150,0 
             [8] = {Point64} 50,150,0 
             [9] = {Point64} 50,50,0 
             [10] = {Point64} 90,50,0 
             [11] = {Point64} 90,0,0 
             [12] = {Point64} 0,0,0 
           */
        
        PathsD simplifyResc = GeoWrangler.close(simplifyRes);
        
        /* Expected output
           simplifyResc = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 13
             [0] = {Point64} 90,0,0 
             [1] = {Point64} 0,0,0 
             [2] = {Point64} 0,200,0 
             [3] = {Point64} 200,200,0 
             [4] = {Point64} 200,0,0 
             [5] = {Point64} 110,0,0 
             [6] = {Point64} 110,50,0 
             [7] = {Point64} 150,50,0 
             [8] = {Point64} 150,150,0 
             [9] = {Point64} 50,150,0 
             [10] = {Point64} 50,50,0 
             [11] = {Point64} 90,50,0 
             [12] = {Point64} 90,0,0 
           */
        
        PathsD simplifyResc_kH = GeoWrangler.makeKeyHole(simplifyResc, reverseEval:false, biDirectionalEval:true);
        
        /* Expected output
           simplifyResc_kH = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 13
             [0] = {Point64} 0,0,0 
             [1] = {Point64} 0,200,0 
             [2] = {Point64} 200,200,0 
             [3] = {Point64} 200,0,0 
             [4] = {Point64} 110,0,0 
             [5] = {Point64} 110,50,0 
             [6] = {Point64} 150,50,0 
             [7] = {Point64} 150,150,0 
             [8] = {Point64} 50,150,0 
             [9] = {Point64} 50,50,0 
             [10] = {Point64} 90,50,0 
             [11] = {Point64} 90,0,0 
             [12] = {Point64} 0,0,0 
           */

        PathsD simplifyRes2 = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, simplifyRes2);
        
        /* Expected output
           simplifyRes2 = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 13
             [0] = {Point64} 90,0,0 
             [1] = {Point64} 0,0,0 
             [2] = {Point64} 0,200,0 
             [3] = {Point64} 200,200,0 
             [4] = {Point64} 200,0,0 
             [5] = {Point64} 110,0,0 
             [6] = {Point64} 110,50,0 
             [7] = {Point64} 150,50,0 
             [8] = {Point64} 150,150,0 
             [9] = {Point64} 50,150,0 
             [10] = {Point64} 50,50,0 
             [11] = {Point64} 90,50,0 
             [12] = {Point64} 90,0,0 
           */
        
        PathsD simplifyRes2_kH = GeoWrangler.makeKeyHole(simplifyRes2, reverseEval:false, biDirectionalEval:true);
        
        /* Expected output
           simplifyRes2_kH = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 13
             [0] = {Point64} 0,0,0 
             [1] = {Point64} 0,200,0 
             [2] = {Point64} 200,200,0 
             [3] = {Point64} 200,0,0 
             [4] = {Point64} 110,0,0 
             [5] = {Point64} 110,50,0 
             [6] = {Point64} 150,50,0 
             [7] = {Point64} 150,150,0 
             [8] = {Point64} 50,150,0 
             [9] = {Point64} 50,50,0 
             [10] = {Point64} 90,50,0 
             [11] = {Point64} 90,0,0 
             [12] = {Point64} 0,0,0 
           */
        
        PathsD simplifyRes2c = GeoWrangler.close(simplifyRes2);
        
        /* Expected output
           simplifyRes2c = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 13
             [0] = {Point64} 90,0,0 
             [1] = {Point64} 0,0,0 
             [2] = {Point64} 0,200,0 
             [3] = {Point64} 200,200,0 
             [4] = {Point64} 200,0,0 
             [5] = {Point64} 110,0,0 
             [6] = {Point64} 110,50,0 
             [7] = {Point64} 150,50,0 
             [8] = {Point64} 150,150,0 
             [9] = {Point64} 50,150,0 
             [10] = {Point64} 50,50,0 
             [11] = {Point64} 90,50,0 
             [12] = {Point64} 90,0,0 
           */
        
        PathsD simplifyRes2c_kH = GeoWrangler.makeKeyHole(simplifyRes2c, reverseEval:false, biDirectionalEval:true);
        
        /* Expected output
           simplifyRes2c_kH = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 13
             [0] = {Point64} 0,0,0 
             [1] = {Point64} 0,200,0 
             [2] = {Point64} 200,200,0 
             [3] = {Point64} 200,0,0 
             [4] = {Point64} 110,0,0 
             [5] = {Point64} 110,50,0 
             [6] = {Point64} 150,50,0 
             [7] = {Point64} 150,150,0 
             [8] = {Point64} 50,150,0 
             [9] = {Point64} 50,50,0 
             [10] = {Point64} 90,50,0 
             [11] = {Point64} 90,0,0 
             [12] = {Point64} 0,0,0 
           */

        // no good - no result
        PathsD intRes = new();
        c.Execute(ClipType.Intersection, FillRule.EvenOdd, intRes);
        PathsD intRes_kH = GeoWrangler.makeKeyHole(intRes, reverseEval:false, biDirectionalEval:true);
        PathsD intResc = GeoWrangler.close(intRes);
        PathsD intResc_kH = GeoWrangler.makeKeyHole(intResc, reverseEval:false, biDirectionalEval:true);

        // no good - no result
        PathsD intResP = new();
        c.Execute(ClipType.Intersection, FillRule.Positive, intResP);
        PathsD intResP_kH = GeoWrangler.makeKeyHole(intResP, reverseEval:false, biDirectionalEval:true);
        PathsD intResPc = GeoWrangler.close(intResP);
        PathsD intResPc_kH = GeoWrangler.makeKeyHole(intResPc, reverseEval:false, biDirectionalEval:true);

        // no good - no result
        PathsD intResNZ = new();
        c.Execute(ClipType.Intersection, FillRule.NonZero, intResNZ);
        PathsD intResNZ_kH = GeoWrangler.makeKeyHole(intResNZ, reverseEval:false, biDirectionalEval:true);
        PathsD intResNZc = GeoWrangler.close(intResNZ);
        PathsD intResNZc_kH = GeoWrangler.makeKeyHole(intResNZc, reverseEval:false, biDirectionalEval:true);

        RectD bounds = Clipper.GetBounds(new PathsD { outer });
        PathD bb = new()
        {
            new(bounds.left, bounds.bottom),
            new(bounds.left, bounds.top),
            new(bounds.right, bounds.top),
            new(bounds.right, bounds.bottom)
        };

        c.Clear();
        c.AddSubject(bb);
        c.AddClip(outer);

        PathsD intRes2 = new();
        c.Execute(ClipType.Intersection, FillRule.EvenOdd, intRes2);
        
        /* Expected output
           intRes2 = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 13
             [0] = {Point64} 90,0,0 
             [1] = {Point64} 0,0,0 
             [2] = {Point64} 0,200,0 
             [3] = {Point64} 200,200,0 
             [4] = {Point64} 200,0,0 
             [5] = {Point64} 110,0,0 
             [6] = {Point64} 110,50,0 
             [7] = {Point64} 150,50,0 
             [8] = {Point64} 150,150,0 
             [9] = {Point64} 50,150,0 
             [10] = {Point64} 50,50,0 
             [11] = {Point64} 90,50,0 
             [12] = {Point64} 90,0,0 
           */
        
        PathsD intRes2_kH = GeoWrangler.makeKeyHole(intRes2, reverseEval:false, biDirectionalEval:true);
        
        /* Expected output
           intRes2_kH = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 13
             [0] = {Point64} 0,0,0 
             [1] = {Point64} 0,200,0 
             [2] = {Point64} 200,200,0 
             [3] = {Point64} 200,0,0 
             [4] = {Point64} 110,0,0 
             [5] = {Point64} 110,50,0 
             [6] = {Point64} 150,50,0 
             [7] = {Point64} 150,150,0 
             [8] = {Point64} 50,150,0 
             [9] = {Point64} 50,50,0 
             [10] = {Point64} 90,50,0 
             [11] = {Point64} 90,0,0 
             [12] = {Point64} 0,0,0 
           */
        
        PathsD intRes2c = GeoWrangler.close(intRes2);
        
        /* Expected output
           intRes2c = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 13
             [0] = {Point64} 90,0,0 
             [1] = {Point64} 0,0,0 
             [2] = {Point64} 0,200,0 
             [3] = {Point64} 200,200,0 
             [4] = {Point64} 200,0,0 
             [5] = {Point64} 110,0,0 
             [6] = {Point64} 110,50,0 
             [7] = {Point64} 150,50,0 
             [8] = {Point64} 150,150,0 
             [9] = {Point64} 50,150,0 
             [10] = {Point64} 50,50,0 
             [11] = {Point64} 90,50,0 
             [12] = {Point64} 90,0,0 
           */
        
        PathsD intRes2c_kH = GeoWrangler.makeKeyHole(intRes2c, reverseEval:false, biDirectionalEval:true);
        
        /* Expected output
           intRes2c_kH = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 13
             [0] = {Point64} 0,0,0 
             [1] = {Point64} 0,200,0 
             [2] = {Point64} 200,200,0 
             [3] = {Point64} 200,0,0 
             [4] = {Point64} 110,0,0 
             [5] = {Point64} 110,50,0 
             [6] = {Point64} 150,50,0 
             [7] = {Point64} 150,150,0 
             [8] = {Point64} 50,150,0 
             [9] = {Point64} 50,50,0 
             [10] = {Point64} 90,50,0 
             [11] = {Point64} 90,0,0 
             [12] = {Point64} 0,0,0 
           */

        // no good - no result
        PathsD intRes2P = new();
        c.Execute(ClipType.Intersection, FillRule.Positive, intRes2P);
        
        /* Expected output
         No geometry in intRes2P
         */
        
        // No results - no geometry
        PathsD intRes2P_kH = GeoWrangler.makeKeyHole(intRes2P, reverseEval:false, biDirectionalEval:true);
        PathsD intRes2Pc = GeoWrangler.close(intRes2P);
        PathsD intRes2Pc_kH = GeoWrangler.makeKeyHole(intRes2Pc, reverseEval:false, biDirectionalEval:true);

        // seems good - get keyhole
        PathsD intRes2NZ = new();
        c.Execute(ClipType.Intersection, FillRule.NonZero, intRes2NZ);
        
        /* Expected output
           intRes2NZ = {List<List<Point64>>} Count = 2
            [0] = {List<Point64>} Count = 6
             [0] = {Point64} 0,0,0 
             [1] = {Point64} 0,200,0 
             [2] = {Point64} 200,200,0 
             [3] = {Point64} 200,0,0 
             [4] = {Point64} 110,0,0 
             [5] = {Point64} 0,0,0 
            [1] = {List<Point64>} Count = 6
             [0] = {Point64} 90,50,0 
             [1] = {Point64} 150,50,0 
             [2] = {Point64} 150,150,0 
             [3] = {Point64} 50,150,0 
             [4] = {Point64} 50,50,0 
             [5] = {Point64} 90,50,0 
           */
        
        PathsD intRes2NZ_kH = GeoWrangler.makeKeyHole(intRes2NZ, reverseEval:false, biDirectionalEval:true);
        
        /* Expected output
           intRes2NZ_kH = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 13
             [0] = {Point64} 200,0,0 
             [1] = {Point64} 200,49.5,0 
             [2] = {Point64} 149.5,49.5,0 
             [3] = {Point64} 149.5,50,0 
             [4] = {Point64} 50,50,0 
             [5] = {Point64} 50,150,0 
             [6] = {Point64} 150,150,0 
             [7] = {Point64} 150,50.5,0 
             [8] = {Point64} 200,50.5,0 
             [9] = {Point64} 200,200,0 
             [10] = {Point64} 0,200,0 
             [11] = {Point64} 0,0,0 
             [12] = {Point64} 200,0,0 
           */
        
        PathsD intRes2NZc = GeoWrangler.close(intRes2NZ);
        
        /* Expected output
           intRes2NZc = {List<List<Point64>>} Count = 2
            [0] = {List<Point64>} Count = 6
             [0] = {Point64} 0,0,0 
             [1] = {Point64} 0,200,0 
             [2] = {Point64} 200,200,0 
             [3] = {Point64} 200,0,0 
             [4] = {Point64} 110,0,0 
             [5] = {Point64} 0,0,0 
            [1] = {List<Point64>} Count = 6
             [0] = {Point64} 90,50,0 
             [1] = {Point64} 150,50,0 
             [2] = {Point64} 150,150,0 
             [3] = {Point64} 50,150,0 
             [4] = {Point64} 50,50,0 
             [5] = {Point64} 90,50,0 
           */
        
        PathsD intRes2NZc_kH = GeoWrangler.makeKeyHole(intRes2NZc, reverseEval:false, biDirectionalEval:true);
        
        /* Expected output
           intRes2NZc_kH = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 13
             [0] = {Point64} 200,0,0 
             [1] = {Point64} 200,49.5,0 
             [2] = {Point64} 149.5,49.5,0 
             [3] = {Point64} 149.5,50,0 
             [4] = {Point64} 50,50,0 
             [5] = {Point64} 50,150,0 
             [6] = {Point64} 150,150,0 
             [7] = {Point64} 150,50.5,0 
             [8] = {Point64} 200,50.5,0 
             [9] = {Point64} 200,200,0 
             [10] = {Point64} 0,200,0 
             [11] = {Point64} 0,0,0 
             [12] = {Point64} 200,0,0 
           */
        
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

        // Segment the paths to match real-world case.
        /*
        Fragmenter f = new(10000);
        PathD outer_f = f.fragmentPath(outer);

        PathD inner1_f = f.fragmentPath(inner1);
        */
        PathsD kHSource = new()
        {
            outer,
            inner1
        };

        /* Expected
           kHSource = {List<List<Point64>>} Count = 2
            [0] = {List<Point64>} Count = 9
             [0] = {Point64} -200,-200,0 
             [1] = {Point64} 300,-200,0 
             [2] = {Point64} 300,200,0 
             [3] = {Point64} -200,200,0 
             [4] = {Point64} -200,-99.9,0 
             [5] = {Point64} -300,-99.9,0 
             [6] = {Point64} -300,-100.1,0 
             [7] = {Point64} -200,-100.1,0 
             [8] = {Point64} -200,-200,0 
            [1] = {List<Point64>} Count = 5
             [0] = {Point64} 100,-100,0 
             [1] = {Point64} 100,100,0 
             [2] = {Point64} 200,100,0 
             [3] = {Point64} 200,-100,0 
             [4] = {Point64} 100,-100,0 
           */
        
        PathsD kH = GeoWrangler.makeKeyHole(kHSource, reverseEval:false, biDirectionalEval:true);

        /* Expected output
           kH = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 17
             [0] = {Point64} -300,-100.1,0 
             [1] = {Point64} -200,-100.1,0 
             [2] = {Point64} -200,-200,0 
             [3] = {Point64} 300,-200,0 
             [4] = {Point64} 300,-100.5,0 
             [5] = {Point64} 199.5,-100.5,0 
             [6] = {Point64} 199.5,-100,0 
             [7] = {Point64} 100,-100,0 
             [8] = {Point64} 100,100,0 
             [9] = {Point64} 200,100,0 
             [10] = {Point64} 200,-99.5,0 
             [11] = {Point64} 300,-99.5,0 
             [12] = {Point64} 300,200,0 
             [13] = {Point64} -200,200,0 
             [14] = {Point64} -200,-99.9,0 
             [15] = {Point64} -300,-99.9,0 
             [16] = {Point64} -300,-100.1,0 
           */
        
        // Gap removal test
        PathsD gR = GeoWrangler.gapRemoval(kH, 100);

        /* Expected output
           gR = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 17
             [0] = {Point64} -300,-100.1,0 
             [1] = {Point64} -200,-100.1,0 
             [2] = {Point64} -200,-200,0 
             [3] = {Point64} 300,-200,0 
             [4] = {Point64} 300,-100.5,0 
             [5] = {Point64} 199.5,-100.5,0 
             [6] = {Point64} 199.5,-100,0 
             [7] = {Point64} 100,-100,0 
             [8] = {Point64} 100,100,0 
             [9] = {Point64} 200,100,0 
             [10] = {Point64} 200,-99.5,0 
             [11] = {Point64} 300,-99.5,0 
             [12] = {Point64} 300,200,0 
             [13] = {Point64} -200,200,0 
             [14] = {Point64} -200,-99.9,0 
             [15] = {Point64} -300,-99.9,0 
             [16] = {Point64} -300,-100.1,0 
           */
        
        // Sliver removal test
        PathsD sR = GeoWrangler.gapRemoval(kH, -100);
        
        /* Expected output
           sR = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 13
             [0] = {Point64} 300,-200,0 
             [1] = {Point64} 300,-100.5,0 
             [2] = {Point64} 199.5,-100.5,0 
             [3] = {Point64} 199.5,-100,0 
             [4] = {Point64} 100,-100,0 
             [5] = {Point64} 100,100,0 
             [6] = {Point64} 200,100,0 
             [7] = {Point64} 200,-99.5,0 
             [8] = {Point64} 300,-99.5,0 
             [9] = {Point64} 300,200,0 
             [10] = {Point64} -200,200,0 
             [11] = {Point64} -200,-200,0 
             [12] = {Point64} 300,-200,0 
           */
        
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
        
        /* Expected
           kHSource = {List<List<Point64>>} Count = 4
            [0] = {List<Point64>} Count = 5
             [0] = {Point64} -200,-200,0 
             [1] = {Point64} 200,-200,0 
             [2] = {Point64} 200,200,0 
             [3] = {Point64} -200,200,0 
             [4] = {Point64} -200,-200,0 
            [1] = {List<Point64>} Count = 5
             [0] = {Point64} -100,-100,0 
             [1] = {Point64} -100,100,0 
             [2] = {Point64} 100,100,0 
             [3] = {Point64} 100,-100,0 
             [4] = {Point64} -100,-100,0 
            [2] = {List<Point64>} Count = 5
             [0] = {Point64} -200,400,0 
             [1] = {Point64} 200,400,0 
             [2] = {Point64} 200,800,0 
             [3] = {Point64} -200,800,0 
             [4] = {Point64} -200,400,0 
            [3] = {List<Point64>} Count = 5
             [0] = {Point64} -100,500,0 
             [1] = {Point64} -100,700,0 
             [2] = {Point64} 100,700,0 
             [3] = {Point64} 100,500,0 
             [4] = {Point64} -100,500,0 
           */

        // Generate keyholed geometry
        PathsD kH = GeoWrangler.makeKeyHole(kHSource, reverseEval:false, biDirectionalEval:true);

        /* Expected output
           kH = {List<List<Point64>>} Count = 2
            [0] = {List<Point64>} Count = 13
             [0] = {Point64} 200,400,0 
             [1] = {Point64} 200,499.5,0 
             [2] = {Point64} 99.5,499.5,0 
             [3] = {Point64} 99.5,500,0 
             [4] = {Point64} -100,500,0 
             [5] = {Point64} -100,700,0 
             [6] = {Point64} 100,700,0 
             [7] = {Point64} 100,500.5,0 
             [8] = {Point64} 200,500.5,0 
             [9] = {Point64} 200,800,0 
             [10] = {Point64} -200,800,0 
             [11] = {Point64} -200,400,0 
             [12] = {Point64} 200,400,0 
            [1] = {List<Point64>} Count = 13
             [0] = {Point64} 200,-200,0 
             [1] = {Point64} 200,-100.5,0 
             [2] = {Point64} 99.5,-100.5,0 
             [3] = {Point64} 99.5,-100,0 
             [4] = {Point64} -100,-100,0 
             [5] = {Point64} -100,100,0 
             [6] = {Point64} 100,100,0 
             [7] = {Point64} 100,-99.5,0 
             [8] = {Point64} 200,-99.5,0 
             [9] = {Point64} 200,200,0 
             [10] = {Point64} -200,200,0 
             [11] = {Point64} -200,-200,0 
             [12] = {Point64} 200,-200,0 
           */

        // Gap removal test
        PathsD gR = GeoWrangler.gapRemoval(kH, 100);

        /* Expected output
           gR = {List<List<Point64>>} Count = 2
            [0] = {List<Point64>} Count = 13
             [0] = {Point64} 200,400,0 
             [1] = {Point64} 200,499.5,0 
             [2] = {Point64} 99.5,499.5,0 
             [3] = {Point64} 99.5,500,0 
             [4] = {Point64} -100,500,0 
             [5] = {Point64} -100,700,0 
             [6] = {Point64} 100,700,0 
             [7] = {Point64} 100,500.5,0 
             [8] = {Point64} 200,500.5,0 
             [9] = {Point64} 200,800,0 
             [10] = {Point64} -200,800,0 
             [11] = {Point64} -200,400,0 
             [12] = {Point64} 200,400,0 
            [1] = {List<Point64>} Count = 13
             [0] = {Point64} 200,-200,0 
             [1] = {Point64} 200,-100.5,0 
             [2] = {Point64} 99.5,-100.5,0 
             [3] = {Point64} 99.5,-100,0 
             [4] = {Point64} -100,-100,0 
             [5] = {Point64} -100,100,0 
             [6] = {Point64} 100,100,0 
             [7] = {Point64} 100,-99.5,0 
             [8] = {Point64} 200,-99.5,0 
             [9] = {Point64} 200,200,0 
             [10] = {Point64} -200,200,0 
             [11] = {Point64} -200,-200,0 
             [12] = {Point64} 200,-200,0 
           */
        
        // Generate sliver geometry.
        PathsD sL = new();
        ClipperD c = new(Constants.roundingDecimalPrecision);
        c.AddSubject(outer1);
        c.AddSubject(outer2);
        c.AddClip(kH);
        c.Execute(ClipType.Difference, FillRule.EvenOdd, sL);
        
        /* Expected output
           sL = {List<List<Point64>>} Count = 2
            [0] = {List<Point64>} Count = 8
             [0] = {Point64} 99.5,499.5,0 
             [1] = {Point64} 99.5,500,0 
             [2] = {Point64} -100,500,0 
             [3] = {Point64} -100,700,0 
             [4] = {Point64} 100,700,0 
             [5] = {Point64} 100,500.5,0 
             [6] = {Point64} 200,500.5,0 
             [7] = {Point64} 200,499.5,0 
            [1] = {List<Point64>} Count = 8
             [0] = {Point64} 99.5,-100.5,0 
             [1] = {Point64} 99.5,-100,0 
             [2] = {Point64} -100,-100,0 
             [3] = {Point64} -100,100,0 
             [4] = {Point64} 100,100,0 
             [5] = {Point64} 100,-99.5,0 
             [6] = {Point64} 200,-99.5,0 
             [7] = {Point64} 200,-100.5,0 
           */
        
        // Sliver removal test
        PathsD sR = GeoWrangler.gapRemoval(sL, -100);
        
        /* Expected output
           sR = {List<List<Point64>>} Count = 2
            [0] = {List<Point64>} Count = 9
             [0] = {Point64} -100,500,0 
             [1] = {Point64} -100,700,0 
             [2] = {Point64} 100,700,0 
             [3] = {Point64} 100,500.5,0 
             [4] = {Point64} 200,500.5,0 
             [5] = {Point64} 200,499.5,0 
             [6] = {Point64} 99.5,499.5,0 
             [7] = {Point64} 99.5,500,0 
             [8] = {Point64} -100,500,0 
            [1] = {List<Point64>} Count = 9
             [0] = {Point64} -100,-100,0 
             [1] = {Point64} -100,100,0 
             [2] = {Point64} 100,100,0 
             [3] = {Point64} 100,-99.5,0 
             [4] = {Point64} 200,-99.5,0 
             [5] = {Point64} 200,-100.5,0 
             [6] = {Point64} 99.5,-100.5,0 
             [7] = {Point64} 99.5,-100,0 
             [8] = {Point64} -100,-100,0 
           */
    }

    private static void complex_islandTest()
    {
        // Island 1 - mix of sliver and gap at the end.
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

        kHSource.Add(outer2);
        kHSource.Add(inner2);

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

        kHSource.Add(outer3);
        kHSource.Add(inner3_1);
        kHSource.Add(inner3_2);

        /* Expected
           kHSource = {List<List<Point64>>} Count = 7
            [0] = {List<Point64>} Count = 9
             [0] = {Point64} -200,-200,0 
             [1] = {Point64} 300,-200,0 
             [2] = {Point64} 300,200,0 
             [3] = {Point64} -200,200,0 
             [4] = {Point64} -200,-99.9,0 
             [5] = {Point64} -300,-99.9,0 
             [6] = {Point64} -300,-100.1,0 
             [7] = {Point64} -200,-100.1,0 
             [8] = {Point64} -200,-200,0 
            [1] = {List<Point64>} Count = 5
             [0] = {Point64} 100,-100,0 
             [1] = {Point64} 100,100,0 
             [2] = {Point64} 200,100,0 
             [3] = {Point64} 200,-100,0 
             [4] = {Point64} 100,-100,0 
            [2] = {List<Point64>} Count = 5
             [0] = {Point64} -200,400,0 
             [1] = {Point64} 200,400,0 
             [2] = {Point64} 200,800,0 
             [3] = {Point64} -200,800,0 
             [4] = {Point64} -200,400,0 
            [3] = {List<Point64>} Count = 5
             [0] = {Point64} -100,500,0 
             [1] = {Point64} -100,700,0 
             [2] = {Point64} 100,700,0 
             [3] = {Point64} 100,500,0 
             [4] = {Point64} -100,500,0 
            [4] = {List<Point64>} Count = 5
             [0] = {Point64} -300,1200,0 
             [1] = {Point64} 300,1200,0 
             [2] = {Point64} 300,1600,0 
             [3] = {Point64} -300,1600,0 
             [4] = {Point64} -300,1200,0 
            [5] = {List<Point64>} Count = 5
             [0] = {Point64} -200,1300,0 
             [1] = {Point64} -200,1500,0 
             [2] = {Point64} -100,1500,0 
             [3] = {Point64} -100,1300,0 
             [4] = {Point64} -200,1300,0 
            [6] = {List<Point64>} Count = 5
             [0] = {Point64} 100,1300,0 
             [1] = {Point64} 100,1500,0 
             [2] = {Point64} 200,1500,0 
             [3] = {Point64} 200,1300,0 
             [4] = {Point64} 100,1300,0 
           */
        
        PathsD kH = GeoWrangler.makeKeyHole(kHSource, reverseEval:false, biDirectionalEval:true);

        /* Expected output
           kH = {List<List<Point64>>} Count = 3
            [0] = {List<Point64>} Count = 21
             [0] = {Point64} 300,1200,0 
             [1] = {Point64} 300,1299.5,0 
             [2] = {Point64} 199.5,1299.5,0 
             [3] = {Point64} 199.5,1300,0 
             [4] = {Point64} 100,1300,0 
             [5] = {Point64} 100,1500,0 
             [6] = {Point64} 200,1500,0 
             [7] = {Point64} 200,1300.5,0 
             [8] = {Point64} 300,1300.5,0 
             [9] = {Point64} 300,1600,0 
             [10] = {Point64} -300,1600,0 
             [11] = {Point64} -300,1500.5,0 
             [12] = {Point64} -199.5,1500.5,0 
             [13] = {Point64} -199.5,1500,0 
             [14] = {Point64} -100,1500,0 
             [15] = {Point64} -100,1300,0 
             [16] = {Point64} -200,1300,0 
             [17] = {Point64} -200,1499.5,0 
             [18] = {Point64} -300,1499.5,0 
             [19] = {Point64} -300,1200,0 
             [20] = {Point64} 300,1200,0 
            [1] = {List<Point64>} Count = 13
             [0] = {Point64} 200,400,0 
             [1] = {Point64} 200,499.5,0 
             [2] = {Point64} 99.5,499.5,0 
             [3] = {Point64} 99.5,500,0 
             [4] = {Point64} -100,500,0 
             [5] = {Point64} -100,700,0 
             [6] = {Point64} 100,700,0 
             [7] = {Point64} 100,500.5,0 
             [8] = {Point64} 200,500.5,0 
             [9] = {Point64} 200,800,0 
             [10] = {Point64} -200,800,0 
             [11] = {Point64} -200,400,0 
             [12] = {Point64} 200,400,0 
            [2] = {List<Point64>} Count = 17
             [0] = {Point64} -300,-100.1,0 
             [1] = {Point64} -200,-100.1,0 
             [2] = {Point64} -200,-200,0 
             [3] = {Point64} 300,-200,0 
             [4] = {Point64} 300,-100.5,0 
             [5] = {Point64} 199.5,-100.5,0 
             [6] = {Point64} 199.5,-100,0 
             [7] = {Point64} 100,-100,0 
             [8] = {Point64} 100,100,0 
             [9] = {Point64} 200,100,0 
             [10] = {Point64} 200,-99.5,0 
             [11] = {Point64} 300,-99.5,0 
             [12] = {Point64} 300,200,0 
             [13] = {Point64} -200,200,0 
             [14] = {Point64} -200,-99.9,0 
             [15] = {Point64} -300,-99.9,0 
             [16] = {Point64} -300,-100.1,0 
           */
        
        // Gap removal test
        PathsD gR = GeoWrangler.gapRemoval(kH, 100);

        /* Expected output
           gR = {List<List<Point64>>} Count = 3
            [0] = {List<Point64>} Count = 21
             [0] = {Point64} 300,1200,0 
             [1] = {Point64} 300,1299.5,0 
             [2] = {Point64} 199.5,1299.5,0 
             [3] = {Point64} 199.5,1300,0 
             [4] = {Point64} 100,1300,0 
             [5] = {Point64} 100,1500,0 
             [6] = {Point64} 200,1500,0 
             [7] = {Point64} 200,1300.5,0 
             [8] = {Point64} 300,1300.5,0 
             [9] = {Point64} 300,1600,0 
             [10] = {Point64} -300,1600,0 
             [11] = {Point64} -300,1500.5,0 
             [12] = {Point64} -199.5,1500.5,0 
             [13] = {Point64} -199.5,1500,0 
             [14] = {Point64} -100,1500,0 
             [15] = {Point64} -100,1300,0 
             [16] = {Point64} -200,1300,0 
             [17] = {Point64} -200,1499.5,0 
             [18] = {Point64} -300,1499.5,0 
             [19] = {Point64} -300,1200,0 
             [20] = {Point64} 300,1200,0 
            [1] = {List<Point64>} Count = 13
             [0] = {Point64} 200,400,0 
             [1] = {Point64} 200,499.5,0 
             [2] = {Point64} 99.5,499.5,0 
             [3] = {Point64} 99.5,500,0 
             [4] = {Point64} -100,500,0 
             [5] = {Point64} -100,700,0 
             [6] = {Point64} 100,700,0 
             [7] = {Point64} 100,500.5,0 
             [8] = {Point64} 200,500.5,0 
             [9] = {Point64} 200,800,0 
             [10] = {Point64} -200,800,0 
             [11] = {Point64} -200,400,0 
             [12] = {Point64} 200,400,0 
            [2] = {List<Point64>} Count = 17
             [0] = {Point64} -300,-100.1,0 
             [1] = {Point64} -200,-100.1,0 
             [2] = {Point64} -200,-200,0 
             [3] = {Point64} 300,-200,0 
             [4] = {Point64} 300,-100.5,0 
             [5] = {Point64} 199.5,-100.5,0 
             [6] = {Point64} 199.5,-100,0 
             [7] = {Point64} 100,-100,0 
             [8] = {Point64} 100,100,0 
             [9] = {Point64} 200,100,0 
             [10] = {Point64} 200,-99.5,0 
             [11] = {Point64} 300,-99.5,0 
             [12] = {Point64} 300,200,0 
             [13] = {Point64} -200,200,0 
             [14] = {Point64} -200,-99.9,0 
             [15] = {Point64} -300,-99.9,0 
             [16] = {Point64} -300,-100.1,0 
           */
        
        // Sliver removal test
        PathsD sR = GeoWrangler.gapRemoval(kH, -100);
        
        /* Expected output
           sR = {List<List<Point64>>} Count = 3
            [0] = {List<Point64>} Count = 21
             [0] = {Point64} 300,1200,0 
             [1] = {Point64} 300,1299.5,0 
             [2] = {Point64} 199.5,1299.5,0 
             [3] = {Point64} 199.5,1300,0 
             [4] = {Point64} 100,1300,0 
             [5] = {Point64} 100,1500,0 
             [6] = {Point64} 200,1500,0 
             [7] = {Point64} 200,1300.5,0 
             [8] = {Point64} 300,1300.5,0 
             [9] = {Point64} 300,1600,0 
             [10] = {Point64} -300,1600,0 
             [11] = {Point64} -300,1500.5,0 
             [12] = {Point64} -199.5,1500.5,0 
             [13] = {Point64} -199.5,1500,0 
             [14] = {Point64} -100,1500,0 
             [15] = {Point64} -100,1300,0 
             [16] = {Point64} -200,1300,0 
             [17] = {Point64} -200,1499.5,0 
             [18] = {Point64} -300,1499.5,0 
             [19] = {Point64} -300,1200,0 
             [20] = {Point64} 300,1200,0 
            [1] = {List<Point64>} Count = 13
             [0] = {Point64} 200,400,0 
             [1] = {Point64} 200,499.5,0 
             [2] = {Point64} 99.5,499.5,0 
             [3] = {Point64} 99.5,500,0 
             [4] = {Point64} -100,500,0 
             [5] = {Point64} -100,700,0 
             [6] = {Point64} 100,700,0 
             [7] = {Point64} 100,500.5,0 
             [8] = {Point64} 200,500.5,0 
             [9] = {Point64} 200,800,0 
             [10] = {Point64} -200,800,0 
             [11] = {Point64} -200,400,0 
             [12] = {Point64} 200,400,0 
            [2] = {List<Point64>} Count = 13
             [0] = {Point64} 300,-200,0 
             [1] = {Point64} 300,-100.5,0 
             [2] = {Point64} 199.5,-100.5,0 
             [3] = {Point64} 199.5,-100,0 
             [4] = {Point64} 100,-100,0 
             [5] = {Point64} 100,100,0 
             [6] = {Point64} 200,100,0 
             [7] = {Point64} 200,-99.5,0 
             [8] = {Point64} 300,-99.5,0 
             [9] = {Point64} 300,200,0 
             [10] = {Point64} -200,200,0 
             [11] = {Point64} -200,-200,0 
             [12] = {Point64} 300,-200,0 
           */
    }

    private static void multiHoleTest()
    {
        PathsD paths = new();
        PathD path_1 = Clipper.MakePath(new double[]
        {
            50, 0,
            0, 0,
            0, 155,
            50, 155,
            50, 0
        });
        paths.Add(path_1);

        PathD path_2 = Clipper.MakePath(new double[]
        {
            35, 135,
            35, 150,
            5, 150,
            5, 135,
            35, 135
        });
        paths.Add(path_2);

        PathD path_3 = Clipper.MakePath(new double[]
        {
            22, 95,
            22, 125,
            5, 125,
            5, 95,
            22, 95
        });
        paths.Add(path_3);

        PathD path_4 = Clipper.MakePath(new double[]
        {
            35, 45,
            35, 75,
            5, 75,
            5, 45,
            35, 45
        });
        paths.Add(path_4);

        PathD path_5 = Clipper.MakePath(new double[]
        {
            35, 5,
            35, 35,
            5, 35,
            5, 50,
            35, 5
        });
        paths.Add(path_5);

        /* Expected
           paths = {List<List<Point64>>} Count = 5
            [0] = {List<Point64>} Count = 5
             [0] = {Point64} 5000,0,0 
             [1] = {Point64} 0,0,0 
             [2] = {Point64} 0,155000,0 
             [3] = {Point64} 5000,155000,0 
             [4] = {Point64} 5000,0,0 
            [1] = {List<Point64>} Count = 5
             [0] = {Point64} 3500,13500,0 
             [1] = {Point64} 3500,1500,0 
             [2] = {Point64} 5000,1500,0 
             [3] = {Point64} 5000,13500,0 
             [4] = {Point64} 3500,13500,0 
            [2] = {List<Point64>} Count = 5
             [0] = {Point64} 22000,95000,0 
             [1] = {Point64} 22000,12500,0 
             [2] = {Point64} 5000,12500,0 
             [3] = {Point64} 5000,95000,0 
             [4] = {Point64} 22000,95000,0 
            [3] = {List<Point64>} Count = 5
             [0] = {Point64} 3500,45000,0 
             [1] = {Point64} 3500,75000,0 
             [2] = {Point64} 5000,75000,0 
             [3] = {Point64} 5000,45000,0 
             [4] = {Point64} 3500,45000,0 
            [4] = {List<Point64>} Count = 5
             [0] = {Point64} 3500,5000,0 
             [1] = {Point64} 3500,3500,0 
             [2] = {Point64} 5000,3500,0 
             [3] = {Point64} 5000,5000,0 
             [4] = {Point64} 3500,5000,0 
           */
        
        // Generate keyholed geometry
        PathsD kH = GeoWrangler.makeKeyHole(paths, reverseEval:false, biDirectionalEval:true);
        
        /* Expected output
           kH = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 37
             [0] = {Point64} 5000,0,0 
             [1] = {Point64} 5000,155000,0 
             [2] = {Point64} 0,155000,0 
             [3] = {Point64} 0,150500,0 
             [4] = {Point64} 5500,150500,0 
             [5] = {Point64} 5500,1500,0 
             [6] = {Point64} 3500,1500,0 
             [7] = {Point64} 3500,13500,0 
             [8] = {Point64} 5000,13500,0 
             [9] = {Point64} 5000,149500,0 
             [10] = {Point64} 0,149500,0 
             [11] = {Point64} 0,125500,0 
             [12] = {Point64} 5500,125500,0 
             [13] = {Point64} 5500,12500,0 
             [14] = {Point64} 22000,12500,0 
             [15] = {Point64} 22000,95000,0 
             [16] = {Point64} 5000,95000,0 
             [17] = {Point64} 5000,124500,0 
             [18] = {Point64} 0,124500,0 
             [19] = {Point64} 0,75500,0 
             [20] = {Point64} 5500,75500,0 
             [21] = {Point64} 5500,75000,0 
             [22] = {Point64} 3500,75000,0 
             [23] = {Point64} 3500,45000,0 
             [24] = {Point64} 5000,45000,0 
             [25] = {Point64} 5000,74500,0 
             [26] = {Point64} 0,74500,0 
             [27] = {Point64} 0,35500,0 
             [28] = {Point64} 5500,35500,0 
             [29] = {Point64} 5500,3500,0 
             [30] = {Point64} 3500,3500,0 
             [31] = {Point64} 3500,5000,0 
             [32] = {Point64} 5000,5000,0 
             [33] = {Point64} 5000,34500,0 
             [34] = {Point64} 0,34500,0 
             [35] = {Point64} 0,0,0 
             [36] = {Point64} 5000,0,0 
           */
    }
}