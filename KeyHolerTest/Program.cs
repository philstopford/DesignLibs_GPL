using Clipper2Lib;
using geoWrangler;
using System.Collections.Generic;

namespace KeyHolerTest;

using Path = Path64;
using Paths = Paths64;

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
            new(-200000, -200000),
            new(200000, -200000),
            new(200000, 200000),
            new(-200000, 200000),
            new(-200000, -200000)
        };

        Path inner1 = new()
        {
            new(-100000, -100000),
            new(-100000, 100000),
            new(100000, 100000),
            new(100000, -100000),
            new(-100000, -100000)
        };

        // Segment the paths to match real-world case.
        /*
        Fragmenter f = new(10000);
        Path outer_f = f.fragmentPath(outer);

        Path inner1_f = f.fragmentPath(inner1);
        */

        Paths kHSource = new()
        {
            outer,
            inner1
        };

        /* Expect to have

        kHSource = {List<List<Point64>>} Count = 2
         [0] = {List<Point64>} Count = 5
          [0] = {Point64} -200000,-200000,0 
          [1] = {Point64} 200000,-200000,0 
          [2] = {Point64} 200000,200000,0 
          [3] = {Point64} -200000,200000,0 
          [4] = {Point64} -200000,-200000,0 
         [1] = {List<Point64>} Count = 5
          [0] = {Point64} -100000,-100000,0 
          [1] = {Point64} -100000,100000,0 
          [2] = {Point64} 100000,100000,0 
          [3] = {Point64} 100000,-100000,0 
          [4] = {Point64} -100000,-100000,0 
         */

        // Generate keyholed geometry
        Paths kH = GeoWrangler.makeKeyHole(kHSource, reverseEval:false, biDirectionalEval:true);

        /* Expected output

        kH = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 13
          [0] = {Point64} 200000,-200000,0 
          [1] = {Point64} 200000,-100500,0 
          [2] = {Point64} 99500,-100500,0 
          [3] = {Point64} 99500,-100000,0 
          [4] = {Point64} -100000,-100000,0 
          [5] = {Point64} -100000,100000,0 
          [6] = {Point64} 100000,100000,0 
          [7] = {Point64} 100000,-99500,0 
          [8] = {Point64} 200000,-99500,0 
          [9] = {Point64} 200000,200000,0 
          [10] = {Point64} -200000,200000,0 
          [11] = {Point64} -200000,-200000,0 
          [12] = {Point64} 200000,-200000,0         
        */
        
        // Generate sliver geometry.
        Paths sL = new();
        Clipper64 c = new();
        c.AddSubject(outer);
        c.AddClip(kH);
        c.Execute(ClipType.Difference, FillRule.EvenOdd, sL);

        /* Expected output
        sL = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 8
          [0] = {Point64} 99500,-100500,0 
          [1] = {Point64} 99500,-100000,0 
          [2] = {Point64} -100000,-100000,0 
          [3] = {Point64} -100000,100000,0 
          [4] = {Point64} 100000,100000,0 
          [5] = {Point64} 100000,-99500,0 
          [6] = {Point64} 200000,-99500,0 
          [7] = {Point64} 200000,-100500,0 
           */
        
        // Gap removal test
        Paths gR = GeoWrangler.gapRemoval(kH, 100);

        /* Expected output
        gR = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 13
          [0] = {Point64} 200000,-200000,0 
          [1] = {Point64} 200000,-100500,0 
          [2] = {Point64} 99500,-100500,0 
          [3] = {Point64} 99500,-100000,0 
          [4] = {Point64} -100000,-100000,0 
          [5] = {Point64} -100000,100000,0 
          [6] = {Point64} 100000,100000,0 
          [7] = {Point64} 100000,-99500,0 
          [8] = {Point64} 200000,-99500,0 
          [9] = {Point64} 200000,200000,0 
          [10] = {Point64} -200000,200000,0 
          [11] = {Point64} -200000,-200000,0 
          [12] = {Point64} 200000,-200000,0 
         */
    
        // Sliver removal test
        Paths sR = GeoWrangler.gapRemoval(sL, -100);
        
        /* Expected output
        sR = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 9
          [0] = {Point64} -100000,-100000,0 
          [1] = {Point64} -100000,100000,0 
          [2] = {Point64} 100000,100000,0 
          [3] = {Point64} 100000,-99500,0 
          [4] = {Point64} 200000,-99500,0 
          [5] = {Point64} 200000,-100500,0 
          [6] = {Point64} 99500,-100500,0 
          [7] = {Point64} 99500,-100000,0 
          [8] = {Point64} -100000,-100000,0 
           */
    }

    private static void multiTest()
    {
        Path outer = new()
        {
            new(-300000, -200000),
            new(300000, -200000),
            new(300000, 200000),
            new(-300000, 200000),
            new(-300000, -200000)
        };

        Path inner1 = new()
        {
            new(-200000, -100000),
            new(-200000, 100000),
            new(-100000, 100000),
            new(-100000, -100000),
            new(-200000, -100000)
        };

        Path inner2 = new()
        {
            new(100000, -100000),
            new(100000, 100000),
            new(200000, 100000),
            new(200000, -100000),
            new(100000, -100000)
        };

        // Segment the paths to match real-world case.
        /*
        Fragmenter f = new(10000);
        Path outer_f = f.fragmentPath(outer);

        Path inner1_f = f.fragmentPath(inner1);
        Path inner2_f = f.fragmentPath(inner2);
        */
        Paths kHSource = new()
        {
            outer,
            inner1,
            inner2
        };
        
        /* Expected
        kHSource = {List<List<Point64>>} Count = 3
         [0] = {List<Point64>} Count = 5
          [0] = {Point64} -300000,-200000,0 
          [1] = {Point64} 300000,-200000,0 
          [2] = {Point64} 300000,200000,0 
          [3] = {Point64} -300000,200000,0 
          [4] = {Point64} -300000,-200000,0 
         [1] = {List<Point64>} Count = 5
          [0] = {Point64} -200000,-100000,0 
          [1] = {Point64} -200000,100000,0 
          [2] = {Point64} -100000,100000,0 
          [3] = {Point64} -100000,-100000,0 
          [4] = {Point64} -200000,-100000,0 
         [2] = {List<Point64>} Count = 5
          [0] = {Point64} 100000,-100000,0 
          [1] = {Point64} 100000,100000,0 
          [2] = {Point64} 200000,100000,0 
          [3] = {Point64} 200000,-100000,0 
          [4] = {Point64} 100000,-100000,0 
           */

        // Generate keyholed geometry
        Paths kH = GeoWrangler.makeKeyHole(new Paths(kHSource), reverseEval:false, biDirectionalEval:true);

        /* Expected output
        kH = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 21
          [0] = {Point64} 300000,-200000,0 
          [1] = {Point64} 300000,-100500,0 
          [2] = {Point64} 199500,-100500,0 
          [3] = {Point64} 199500,-100000,0 
          [4] = {Point64} 100000,-100000,0 
          [5] = {Point64} 100000,100000,0 
          [6] = {Point64} 200000,100000,0 
          [7] = {Point64} 200000,-99500,0 
          [8] = {Point64} 300000,-99500,0 
          [9] = {Point64} 300000,200000,0 
          [10] = {Point64} -300000,200000,0 
          [11] = {Point64} -300000,100500,0 
          [12] = {Point64} -199500,100500,0 
          [13] = {Point64} -199500,100000,0 
          [14] = {Point64} -100000,100000,0 
          [15] = {Point64} -100000,-100000,0 
          [16] = {Point64} -200000,-100000,0 
          [17] = {Point64} -200000,99500,0 
          [18] = {Point64} -300000,99500,0 
          [19] = {Point64} -300000,-200000,0 
          [20] = {Point64} 300000,-200000,0 
           */
        
        // Generate sliver geometry.
        Paths sL = new();
        Clipper64 c = new();
        c.AddSubject(outer);
        c.AddClip(kH);
        c.Execute(ClipType.Difference, FillRule.EvenOdd, sL);
        
        /* Expected output
        sL = {List<List<Point64>>} Count = 2
         [0] = {List<Point64>} Count = 8
          [0] = {Point64} -300000,100500,0 
          [1] = {Point64} -199500,100500,0 
          [2] = {Point64} -199500,100000,0 
          [3] = {Point64} -100000,100000,0 
          [4] = {Point64} -100000,-100000,0 
          [5] = {Point64} -200000,-100000,0 
          [6] = {Point64} -200000,99500,0 
          [7] = {Point64} -300000,99500,0 
         [1] = {List<Point64>} Count = 8
          [0] = {Point64} 199500,-100500,0 
          [1] = {Point64} 199500,-100000,0 
          [2] = {Point64} 100000,-100000,0 
          [3] = {Point64} 100000,100000,0 
          [4] = {Point64} 200000,100000,0 
          [5] = {Point64} 200000,-99500,0 
          [6] = {Point64} 300000,-99500,0 
          [7] = {Point64} 300000,-100500,0 
           */

        // Gap removal test
        Paths gR = GeoWrangler.gapRemoval(kH, 100);

        /* Expected output
        gR = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 21
          [0] = {Point64} 300000,-200000,0 
          [1] = {Point64} 300000,-100500,0 
          [2] = {Point64} 199500,-100500,0 
          [3] = {Point64} 199500,-100000,0 
          [4] = {Point64} 100000,-100000,0 
          [5] = {Point64} 100000,100000,0 
          [6] = {Point64} 200000,100000,0 
          [7] = {Point64} 200000,-99500,0 
          [8] = {Point64} 300000,-99500,0 
          [9] = {Point64} 300000,200000,0 
          [10] = {Point64} -300000,200000,0 
          [11] = {Point64} -300000,100500,0 
          [12] = {Point64} -199500,100500,0 
          [13] = {Point64} -199500,100000,0 
          [14] = {Point64} -100000,100000,0 
          [15] = {Point64} -100000,-100000,0 
          [16] = {Point64} -200000,-100000,0 
          [17] = {Point64} -200000,99500,0 
          [18] = {Point64} -300000,99500,0 
          [19] = {Point64} -300000,-200000,0 
          [20] = {Point64} 300000,-200000,0 
           */
        
        // Sliver removal test
        Paths sR = GeoWrangler.gapRemoval(sL, -100);
        
        /* Expected output
        sR = {List<List<Point64>>} Count = 2
         [0] = {List<Point64>} Count = 9
          [0] = {Point64} -300000,99500,0 
          [1] = {Point64} -300000,100500,0 
          [2] = {Point64} -199500,100500,0 
          [3] = {Point64} -199500,100000,0 
          [4] = {Point64} -100000,100000,0 
          [5] = {Point64} -100000,-100000,0 
          [6] = {Point64} -200000,-100000,0 
          [7] = {Point64} -200000,99500,0 
          [8] = {Point64} -300000,99500,0 
         [1] = {List<Point64>} Count = 9
          [0] = {Point64} 100000,-100000,0 
          [1] = {Point64} 100000,100000,0 
          [2] = {Point64} 200000,100000,0 
          [3] = {Point64} 200000,-99500,0 
          [4] = {Point64} 300000,-99500,0 
          [5] = {Point64} 300000,-100500,0 
          [6] = {Point64} 199500,-100500,0 
          [7] = {Point64} 199500,-100000,0 
          [8] = {Point64} 100000,-100000,0 
           */
    }

    private static void multiCutTest()
    {
        Path outer = new()
        {
            new(0, 0),
            new(400000, 00000),
            new(400000, 400000),
            new(0, 400000),
            new(0, 0)
        };

        Path inner1 = new()
        {
            new(50000, 150000),
            new(50000, 250000),
            new(350000, 250000),
            new(350000, 150000),
            new(50000, 150000)
        };

        Path inner2 = new()
        {
            new(150000, 50000),
            new(150000, 350000),
            new(250000, 350000),
            new(250000, 50000),
            new(150000, 50000)
        };

        Paths kHSource = new()
        {
            outer,
            inner1
        };
        
        /* Expected
        kHSource = {List<List<Point64>>} Count = 2
         [0] = {List<Point64>} Count = 5
          [0] = {Point64} 0,0,0 
          [1] = {Point64} 400000,0,0 
          [2] = {Point64} 400000,400000,0 
          [3] = {Point64} 0,400000,0 
          [4] = {Point64} 0,0,0 
         [1] = {List<Point64>} Count = 5
          [0] = {Point64} 50000,150000,0 
          [1] = {Point64} 50000,250000,0 
          [2] = {Point64} 350000,250000,0 
          [3] = {Point64} 350000,150000,0 
          [4] = {Point64} 50000,150000,0 
           */

        // Generate keyholed geometry
        Paths kH = GeoWrangler.makeKeyHole(new Paths(kHSource), reverseEval:false, biDirectionalEval:true);

        /* Expected output
        kH = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 13
          [0] = {Point64} 400000,0,0 
          [1] = {Point64} 400000,149500,0 
          [2] = {Point64} 349500,149500,0 
          [3] = {Point64} 349500,150000,0 
          [4] = {Point64} 50000,150000,0 
          [5] = {Point64} 50000,250000,0 
          [6] = {Point64} 350000,250000,0 
          [7] = {Point64} 350000,150500,0 
          [8] = {Point64} 400000,150500,0 
          [9] = {Point64} 400000,400000,0 
          [10] = {Point64} 0,400000,0 
          [11] = {Point64} 0,0,0 
          [12] = {Point64} 400000,0,0 
           */
        
        Paths kHSource2 = new();
        kHSource2.AddRange(kH);
        kHSource2.Add(inner2);

        /* Expected
        kHSource2 = {List<List<Point64>>} Count = 2
         [0] = {List<Point64>} Count = 13
          [0] = {Point64} 400000,0,0 
          [1] = {Point64} 400000,149500,0 
          [2] = {Point64} 349500,149500,0 
          [3] = {Point64} 349500,150000,0 
          [4] = {Point64} 50000,150000,0 
          [5] = {Point64} 50000,250000,0 
          [6] = {Point64} 350000,250000,0 
          [7] = {Point64} 350000,150500,0 
          [8] = {Point64} 400000,150500,0 
          [9] = {Point64} 400000,400000,0 
          [10] = {Point64} 0,400000,0 
          [11] = {Point64} 0,0,0 
          [12] = {Point64} 400000,0,0 
         [1] = {List<Point64>} Count = 5
          [0] = {Point64} 150000,50000,0 
          [1] = {Point64} 150000,350000,0 
          [2] = {Point64} 250000,350000,0 
          [3] = {Point64} 250000,50000,0 
          [4] = {Point64} 150000,50000,0 
           */
        
        // Generate keyholed geometry
        Paths kH2 = GeoWrangler.makeKeyHole(new Paths(kHSource2), reverseEval:false, biDirectionalEval:true);

        /* Expected output
        kH2 = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 21
          [0] = {Point64} 400000,0,0 
          [1] = {Point64} 400000,149500,0 
          [2] = {Point64} 349500,149500,0 
          [3] = {Point64} 349500,150000,0 
          [4] = {Point64} 250000,150000,0 
          [5] = {Point64} 250000,50000,0 
          [6] = {Point64} 150000,50000,0 
          [7] = {Point64} 150000,150000,0 
          [8] = {Point64} 50000,150000,0 
          [9] = {Point64} 50000,250000,0 
          [10] = {Point64} 150000,250000,0 
          [11] = {Point64} 150000,350000,0 
          [12] = {Point64} 250000,350000,0 
          [13] = {Point64} 250000,250000,0 
          [14] = {Point64} 350000,250000,0 
          [15] = {Point64} 350000,150500,0 
          [16] = {Point64} 400000,150500,0 
          [17] = {Point64} 400000,400000,0 
          [18] = {Point64} 0,400000,0 
          [19] = {Point64} 0,0,0 
          [20] = {Point64} 400000,0,0 
           */
        
        // Generate sliver geometry.
        Clipper64 c = new();
        Paths sL = new();
        c.AddSubject(outer);
        c.AddClip(kH);
        c.Execute(ClipType.Difference, FillRule.EvenOdd, sL);
        
        /* Expected output
        sL = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 8
          [0] = {Point64} 349500,149500,0 
          [1] = {Point64} 349500,150000,0 
          [2] = {Point64} 50000,150000,0 
          [3] = {Point64} 50000,250000,0 
          [4] = {Point64} 350000,250000,0 
          [5] = {Point64} 350000,150500,0 
          [6] = {Point64} 400000,150500,0 
          [7] = {Point64} 400000,149500,0 
           */

        c.Clear();
        Paths sL2 = new();
        c.AddSubject(outer);
        c.AddClip(kH2);
        c.Execute(ClipType.Difference, FillRule.EvenOdd, sL2);
        
        /* Expected output
        sL2 = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 16
          [0] = {Point64} 349500,149500,0 
          [1] = {Point64} 349500,150000,0 
          [2] = {Point64} 250000,150000,0 
          [3] = {Point64} 250000,50000,0 
          [4] = {Point64} 150000,50000,0 
          [5] = {Point64} 150000,150000,0 
          [6] = {Point64} 50000,150000,0 
          [7] = {Point64} 50000,250000,0 
          [8] = {Point64} 150000,250000,0 
          [9] = {Point64} 150000,350000,0 
          [10] = {Point64} 250000,350000,0 
          [11] = {Point64} 250000,250000,0 
          [12] = {Point64} 350000,250000,0 
          [13] = {Point64} 350000,150500,0 
          [14] = {Point64} 400000,150500,0 
          [15] = {Point64} 400000,149500,0 
           */

        // Gap removal test
        Paths gR = GeoWrangler.gapRemoval(kH, 100);
        
        /* Expected output
        gR = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 13
          [0] = {Point64} 400000,0,0 
          [1] = {Point64} 400000,149500,0 
          [2] = {Point64} 349500,149500,0 
          [3] = {Point64} 349500,150000,0 
          [4] = {Point64} 50000,150000,0 
          [5] = {Point64} 50000,250000,0 
          [6] = {Point64} 350000,250000,0 
          [7] = {Point64} 350000,150500,0 
          [8] = {Point64} 400000,150500,0 
          [9] = {Point64} 400000,400000,0 
          [10] = {Point64} 0,400000,0 
          [11] = {Point64} 0,0,0 
          [12] = {Point64} 400000,0,0 
           */
        
        Paths gR2 = GeoWrangler.gapRemoval(kH2, 100);
        
        /* Expected output
        gR2 = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 21
          [0] = {Point64} 400000,0,0 
          [1] = {Point64} 400000,149500,0 
          [2] = {Point64} 349500,149500,0 
          [3] = {Point64} 349500,150000,0 
          [4] = {Point64} 250000,150000,0 
          [5] = {Point64} 250000,50000,0 
          [6] = {Point64} 150000,50000,0 
          [7] = {Point64} 150000,150000,0 
          [8] = {Point64} 50000,150000,0 
          [9] = {Point64} 50000,250000,0 
          [10] = {Point64} 150000,250000,0 
          [11] = {Point64} 150000,350000,0 
          [12] = {Point64} 250000,350000,0 
          [13] = {Point64} 250000,250000,0 
          [14] = {Point64} 350000,250000,0 
          [15] = {Point64} 350000,150500,0 
          [16] = {Point64} 400000,150500,0 
          [17] = {Point64} 400000,400000,0 
          [18] = {Point64} 0,400000,0 
          [19] = {Point64} 0,0,0 
          [20] = {Point64} 400000,0,0 
           */

        Paths kHSource3 = new();
        kHSource3.AddRange(gR);

        Paths kHSource4 = new();
        kHSource4.AddRange(gR2);

        // Generate keyholed geometry
        Paths kH3 = GeoWrangler.makeKeyHole(new Paths(kHSource3), reverseEval:false, biDirectionalEval:true);
        
        /* Expected output
        kH3 = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 13
          [0] = {Point64} 400000,0,0 
          [1] = {Point64} 400000,149500,0 
          [2] = {Point64} 349500,149500,0 
          [3] = {Point64} 349500,150000,0 
          [4] = {Point64} 50000,150000,0 
          [5] = {Point64} 50000,250000,0 
          [6] = {Point64} 350000,250000,0 
          [7] = {Point64} 350000,150500,0 
          [8] = {Point64} 400000,150500,0 
          [9] = {Point64} 400000,400000,0 
          [10] = {Point64} 0,400000,0 
          [11] = {Point64} 0,0,0 
          [12] = {Point64} 400000,0,0 
           */
        
        Paths kH4 = GeoWrangler.makeKeyHole(new Paths(kHSource4), reverseEval:false, biDirectionalEval:true);

        /* Expected output
        kH4 = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 21
          [0] = {Point64} 400000,0,0 
          [1] = {Point64} 400000,149500,0 
          [2] = {Point64} 349500,149500,0 
          [3] = {Point64} 349500,150000,0 
          [4] = {Point64} 250000,150000,0 
          [5] = {Point64} 250000,50000,0 
          [6] = {Point64} 150000,50000,0 
          [7] = {Point64} 150000,150000,0 
          [8] = {Point64} 50000,150000,0 
          [9] = {Point64} 50000,250000,0 
          [10] = {Point64} 150000,250000,0 
          [11] = {Point64} 150000,350000,0 
          [12] = {Point64} 250000,350000,0 
          [13] = {Point64} 250000,250000,0 
          [14] = {Point64} 350000,250000,0 
          [15] = {Point64} 350000,150500,0 
          [16] = {Point64} 400000,150500,0 
          [17] = {Point64} 400000,400000,0 
          [18] = {Point64} 0,400000,0 
          [19] = {Point64} 0,0,0 
          [20] = {Point64} 400000,0,0 
           */
        
        // Sliver removal test
        Paths sR = GeoWrangler.gapRemoval(sL, -100);
        
        /* Expected output
        sR = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 9
          [0] = {Point64} 50000,150000,0 
          [1] = {Point64} 50000,250000,0 
          [2] = {Point64} 350000,250000,0 
          [3] = {Point64} 350000,150500,0 
          [4] = {Point64} 400000,150500,0 
          [5] = {Point64} 400000,149500,0 
          [6] = {Point64} 349500,149500,0 
          [7] = {Point64} 349500,150000,0 
          [8] = {Point64} 50000,150000,0 
           */
        
        Paths sR2 = GeoWrangler.gapRemoval(sL2, -100);
        
        /* Expected output
        sR2 = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 17
          [0] = {Point64} 50000,150000,0 
          [1] = {Point64} 50000,250000,0 
          [2] = {Point64} 150000,250000,0 
          [3] = {Point64} 150000,350000,0 
          [4] = {Point64} 250000,350000,0 
          [5] = {Point64} 250000,250000,0 
          [6] = {Point64} 350000,250000,0 
          [7] = {Point64} 350000,150500,0 
          [8] = {Point64} 400000,150500,0 
          [9] = {Point64} 400000,149500,0 
          [10] = {Point64} 349500,149500,0 
          [11] = {Point64} 349500,150000,0 
          [12] = {Point64} 250000,150000,0 
          [13] = {Point64} 250000,50000,0 
          [14] = {Point64} 150000,50000,0 
          [15] = {Point64} 150000,150000,0 
          [16] = {Point64} 50000,150000,0 
           */
    }

    private static void selfOverlapTest()
    {
        Path outer = new()
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
        Paths dSource = new() {outer};
        Paths decomp = GeoWrangler.decompose(dSource);
        
        /* Expected output
        decomp = {List<List<Point64>>} Count = 2
         [0] = {List<Point64>} Count = 6
          [0] = {Point64} 0,0,0 
          [1] = {Point64} 0,200000,0 
          [2] = {Point64} 200000,200000,0 
          [3] = {Point64} 200000,0,0 
          [4] = {Point64} 110000,0,0 
          [5] = {Point64} 0,0,0 
         [1] = {List<Point64>} Count = 6
          [0] = {Point64} 50000,50000,0 
          [1] = {Point64} 90000,50000,0 
          [2] = {Point64} 150000,50000,0 
          [3] = {Point64} 150000,150000,0 
          [4] = {Point64} 50000,150000,0 
          [5] = {Point64} 50000,50000,0 
           */
        
        Paths kHD = GeoWrangler.makeKeyHole(dSource, reverseEval:false, biDirectionalEval:true);

        /* Expected output
        kHD = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 13
          [0] = {Point64} 200000,0,0 
          [1] = {Point64} 200000,49500,0 
          [2] = {Point64} 149500,49500,0 
          [3] = {Point64} 149500,50000,0 
          [4] = {Point64} 50000,50000,0 
          [5] = {Point64} 50000,150000,0 
          [6] = {Point64} 150000,150000,0 
          [7] = {Point64} 150000,50500,0 
          [8] = {Point64} 200000,50500,0 
          [9] = {Point64} 200000,200000,0 
          [10] = {Point64} 0,200000,0 
          [11] = {Point64} 0,0,0 
          [12] = {Point64} 200000,0,0 
           */
        
        // keyholer test
        Paths kHSource = new() {outer};
        Paths kH = GeoWrangler.makeKeyHole(kHSource, reverseEval:false, biDirectionalEval:true);

        /* Expected output
        kH = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 13
          [0] = {Point64} 200000,0,0 
          [1] = {Point64} 200000,49500,0 
          [2] = {Point64} 149500,49500,0 
          [3] = {Point64} 149500,50000,0 
          [4] = {Point64} 50000,50000,0 
          [5] = {Point64} 50000,150000,0 
          [6] = {Point64} 150000,150000,0 
          [7] = {Point64} 150000,50500,0 
          [8] = {Point64} 200000,50500,0 
          [9] = {Point64} 200000,200000,0 
          [10] = {Point64} 0,200000,0 
          [11] = {Point64} 0,0,0 
          [12] = {Point64} 200000,0,0 
           */
        
        Clipper64 c = new();
        c.AddSubject(outer);

        Paths unionRes = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, unionRes);
        
        /* Expected output
        unionRes = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 13
          [0] = {Point64} 110000,50000,0 
          [1] = {Point64} 150000,50000,0 
          [2] = {Point64} 150000,150000,0 
          [3] = {Point64} 50000,150000,0 
          [4] = {Point64} 50000,50000,0 
          [5] = {Point64} 90000,50000,0 
          [6] = {Point64} 90000,0,0 
          [7] = {Point64} 0,0,0 
          [8] = {Point64} 0,200000,0 
          [9] = {Point64} 200000,200000,0 
          [10] = {Point64} 200000,0,0 
          [11] = {Point64} 110000,0,0 
          [12] = {Point64} 110000,50000,0 
           */
        
        Paths unionRes_kH = GeoWrangler.makeKeyHole(unionRes, reverseEval:false, biDirectionalEval:true);
        
        /* Expected output
        unionRes_kH = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 13
          [0] = {Point64} 0,0,0 
          [1] = {Point64} 0,200000,0 
          [2] = {Point64} 200000,200000,0 
          [3] = {Point64} 200000,0,0 
          [4] = {Point64} 110000,0,0 
          [5] = {Point64} 110000,50000,0 
          [6] = {Point64} 150000,50000,0 
          [7] = {Point64} 150000,150000,0 
          [8] = {Point64} 50000,150000,0 
          [9] = {Point64} 50000,50000,0 
          [10] = {Point64} 90000,50000,0 
          [11] = {Point64} 90000,0,0 
          [12] = {Point64} 0,0,0 
           */
        
        Paths unionResc = GeoWrangler.close(unionRes);
        Paths unionResc_kH = GeoWrangler.makeKeyHole(unionResc, reverseEval:false, biDirectionalEval:true);

        /* Expected output
        unionResc_kH = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 13
          [0] = {Point64} 0,0,0 
          [1] = {Point64} 0,200000,0 
          [2] = {Point64} 200000,200000,0 
          [3] = {Point64} 200000,0,0 
          [4] = {Point64} 110000,0,0 
          [5] = {Point64} 110000,50000,0 
          [6] = {Point64} 150000,50000,0 
          [7] = {Point64} 150000,150000,0 
          [8] = {Point64} 50000,150000,0 
          [9] = {Point64} 50000,50000,0 
          [10] = {Point64} 90000,50000,0 
          [11] = {Point64} 90000,0,0 
          [12] = {Point64} 0,0,0 
           */
        
        Paths unionResP = new();
        c.Execute(ClipType.Union, FillRule.Positive, unionResP);
        
        /* Expected output
         No geometry in unionResP
         */
        
        // no keyhole for any of the below
        Paths unionResP_kH = GeoWrangler.makeKeyHole(unionResP, reverseEval:false, biDirectionalEval:true);
        Paths unionResPc = GeoWrangler.close(unionResP);
        Paths unionResPc_kH = GeoWrangler.makeKeyHole(unionResPc, reverseEval:false, biDirectionalEval:true);

        // seems good - get keyhole
        Paths unionResNZ = new();
        c.Execute(ClipType.Union, FillRule.NonZero, unionResNZ);
        
        /* Expected output
        unionResNZ = {List<List<Point64>>} Count = 2
         [0] = {List<Point64>} Count = 6
          [0] = {Point64} 0,0,0 
          [1] = {Point64} 0,200000,0 
          [2] = {Point64} 200000,200000,0 
          [3] = {Point64} 200000,0,0 
          [4] = {Point64} 110000,0,0 
          [5] = {Point64} 0,0,0 
         [1] = {List<Point64>} Count = 6
          [0] = {Point64} 90000,50000,0 
          [1] = {Point64} 150000,50000,0 
          [2] = {Point64} 150000,150000,0 
          [3] = {Point64} 50000,150000,0 
          [4] = {Point64} 50000,50000,0 
          [5] = {Point64} 90000,50000,0 
           */
        
        Paths unionResNZ_kH = GeoWrangler.makeKeyHole(unionResNZ, reverseEval:false, biDirectionalEval:true);
        
        /* Expected output
        unionResNZ_kH = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 13
          [0] = {Point64} 200000,0,0 
          [1] = {Point64} 200000,49500,0 
          [2] = {Point64} 149500,49500,0 
          [3] = {Point64} 149500,50000,0 
          [4] = {Point64} 50000,50000,0 
          [5] = {Point64} 50000,150000,0 
          [6] = {Point64} 150000,150000,0 
          [7] = {Point64} 150000,50500,0 
          [8] = {Point64} 200000,50500,0 
          [9] = {Point64} 200000,200000,0 
          [10] = {Point64} 0,200000,0 
          [11] = {Point64} 0,0,0 
          [12] = {Point64} 200000,0,0 
           */
        
        Paths unionResNZc = GeoWrangler.close(unionResNZ);
        
        /* Expected output
        unionResNZc = {List<List<Point64>>} Count = 2
         [0] = {List<Point64>} Count = 6
          [0] = {Point64} 0,0,0 
          [1] = {Point64} 0,200000,0 
          [2] = {Point64} 200000,200000,0 
          [3] = {Point64} 200000,0,0 
          [4] = {Point64} 110000,0,0 
          [5] = {Point64} 0,0,0 
         [1] = {List<Point64>} Count = 6
          [0] = {Point64} 90000,50000,0 
          [1] = {Point64} 150000,50000,0 
          [2] = {Point64} 150000,150000,0 
          [3] = {Point64} 50000,150000,0 
          [4] = {Point64} 50000,50000,0 
          [5] = {Point64} 90000,50000,0 
           */
        
        Paths unionResNZc_kH = GeoWrangler.makeKeyHole(unionResNZc, reverseEval:false, biDirectionalEval:true);
        
        /* Expected result
        unionResNZc_kH = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 13
          [0] = {Point64} 200000,0,0 
          [1] = {Point64} 200000,49500,0 
          [2] = {Point64} 149500,49500,0 
          [3] = {Point64} 149500,50000,0 
          [4] = {Point64} 50000,50000,0 
          [5] = {Point64} 50000,150000,0 
          [6] = {Point64} 150000,150000,0 
          [7] = {Point64} 150000,50500,0 
          [8] = {Point64} 200000,50500,0 
          [9] = {Point64} 200000,200000,0 
          [10] = {Point64} 0,200000,0 
          [11] = {Point64} 0,0,0 
          [12] = {Point64} 200000,0,0 
           */

        Paths simplifyRes = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, simplifyRes);
        simplifyRes = GeoWrangler.stripColinear(simplifyRes);

        /* Expected output
        simplifyRes = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 13
          [0] = {Point64} 110000,50000,0 
          [1] = {Point64} 150000,50000,0 
          [2] = {Point64} 150000,150000,0 
          [3] = {Point64} 50000,150000,0 
          [4] = {Point64} 50000,50000,0 
          [5] = {Point64} 90000,50000,0 
          [6] = {Point64} 90000,0,0 
          [7] = {Point64} 0,0,0 
          [8] = {Point64} 0,200000,0 
          [9] = {Point64} 200000,200000,0 
          [10] = {Point64} 200000,0,0 
          [11] = {Point64} 110000,0,0 
          [12] = {Point64} 110000,50000,0 
           */

        Paths simplifyRes_kH = GeoWrangler.makeKeyHole(simplifyRes, reverseEval:false, biDirectionalEval:true);
        
        /* Expected output
        simplifyRes_kH = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 13
          [0] = {Point64} 0,0,0 
          [1] = {Point64} 0,200000,0 
          [2] = {Point64} 200000,200000,0 
          [3] = {Point64} 200000,0,0 
          [4] = {Point64} 110000,0,0 
          [5] = {Point64} 110000,50000,0 
          [6] = {Point64} 150000,50000,0 
          [7] = {Point64} 150000,150000,0 
          [8] = {Point64} 50000,150000,0 
          [9] = {Point64} 50000,50000,0 
          [10] = {Point64} 90000,50000,0 
          [11] = {Point64} 90000,0,0 
          [12] = {Point64} 0,0,0 
           */
        
        Paths simplifyResc = GeoWrangler.close(simplifyRes);
        
        /* Expected output
        simplifyResc = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 13
          [0] = {Point64} 110000,50000,0 
          [1] = {Point64} 150000,50000,0 
          [2] = {Point64} 150000,150000,0 
          [3] = {Point64} 50000,150000,0 
          [4] = {Point64} 50000,50000,0 
          [5] = {Point64} 90000,50000,0 
          [6] = {Point64} 90000,0,0 
          [7] = {Point64} 0,0,0 
          [8] = {Point64} 0,200000,0 
          [9] = {Point64} 200000,200000,0 
          [10] = {Point64} 200000,0,0 
          [11] = {Point64} 110000,0,0 
          [12] = {Point64} 110000,50000,0 
           */
        
        Paths simplifyResc_kH = GeoWrangler.makeKeyHole(simplifyResc, reverseEval:false, biDirectionalEval:true);

        /* Expected output
        simplifyResc_kH = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 13
          [0] = {Point64} 0,0,0 
          [1] = {Point64} 0,200000,0 
          [2] = {Point64} 200000,200000,0 
          [3] = {Point64} 200000,0,0 
          [4] = {Point64} 110000,0,0 
          [5] = {Point64} 110000,50000,0 
          [6] = {Point64} 150000,50000,0 
          [7] = {Point64} 150000,150000,0 
          [8] = {Point64} 50000,150000,0 
          [9] = {Point64} 50000,50000,0 
          [10] = {Point64} 90000,50000,0 
          [11] = {Point64} 90000,0,0 
          [12] = {Point64} 0,0,0 
           */
        
        Paths simplifyRes2 = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, simplifyRes2);
        Paths simplifyRes2_kH = GeoWrangler.makeKeyHole(simplifyRes2, reverseEval:false, biDirectionalEval:true);
        Paths simplifyRes2c = GeoWrangler.close(simplifyRes2);
        Paths simplifyRes2c_kH = GeoWrangler.makeKeyHole(simplifyRes2c, reverseEval:false, biDirectionalEval:true);

        // no good - no result
        Paths intRes = new();
        c.Execute(ClipType.Intersection, FillRule.EvenOdd, intRes);
        Paths intRes_kH = GeoWrangler.makeKeyHole(intRes, reverseEval:false, biDirectionalEval:true);
        Paths intResc = GeoWrangler.close(intRes);
        Paths intResc_kH = GeoWrangler.makeKeyHole(intResc, reverseEval:false, biDirectionalEval:true);

        // no good - no result
        Paths intResP = new();
        c.Execute(ClipType.Intersection, FillRule.Positive, intResP);
        Paths intResP_kH = GeoWrangler.makeKeyHole(intResP, reverseEval:false, biDirectionalEval:true);
        Paths intResPc = GeoWrangler.close(intResP);
        Paths intResPc_kH = GeoWrangler.makeKeyHole(intResPc, reverseEval:false, biDirectionalEval:true);

        // no good - no result
        Paths intResNZ = new();
        c.Execute(ClipType.Intersection, FillRule.NonZero, intResNZ);
        Paths intResNZ_kH = GeoWrangler.makeKeyHole(intResNZ, reverseEval:false, biDirectionalEval:true);
        Paths intResNZc = GeoWrangler.close(intResNZ);
        Paths intResNZc_kH = GeoWrangler.makeKeyHole(intResNZc, reverseEval:false, biDirectionalEval:true);

        Rect64 bounds = Clipper.GetBounds(new Paths { outer });
        Path bb = new()
        {
            new(bounds.left, bounds.bottom),
            new(bounds.left, bounds.top),
            new(bounds.right, bounds.top),
            new(bounds.right, bounds.bottom)
        };

        c.Clear();
        c.AddSubject(bb);
        c.AddClip(outer);

        Paths intRes2 = new();
        c.Execute(ClipType.Intersection, FillRule.EvenOdd, intRes2);
        
        /* Expected output
         intRes2 = {List<List<Point64>>} Count = 2
          [0] = {List<Point64>} Count = 6
           [0] = {Point64} 0,0,0 
           [1] = {Point64} 0,200000,0 
           [2] = {Point64} 200000,200000,0 
           [3] = {Point64} 200000,0,0 
           [4] = {Point64} 110000,0,0 
           [5] = {Point64} 0,0,0 
          [1] = {List<Point64>} Count = 9
           [0] = {Point64} 90000,0,0 
           [1] = {Point64} 110000,0,0 
           [2] = {Point64} 110000,50000,0 
           [3] = {Point64} 150000,50000,0 
           [4] = {Point64} 150000,150000,0 
           [5] = {Point64} 50000,150000,0 
           [6] = {Point64} 50000,50000,0 
           [7] = {Point64} 90000,50000,0 
           [8] = {Point64} 90000,0,0 
           */
        
        Paths intRes2_kH = GeoWrangler.makeKeyHole(intRes2, reverseEval:false, biDirectionalEval:true);
        
        /* Expected output
         intRes2_kH = {List<List<Point64>>} Count = 1
          [0] = {List<Point64>} Count = 13
           [0] = {Point64} 0,0,0 
           [1] = {Point64} 0,200000,0 
           [2] = {Point64} 200000,200000,0 
           [3] = {Point64} 200000,0,0 
           [4] = {Point64} 110000,0,0 
           [5] = {Point64} 110000,50000,0 
           [6] = {Point64} 150000,50000,0 
           [7] = {Point64} 150000,150000,0 
           [8] = {Point64} 50000,150000,0 
           [9] = {Point64} 50000,50000,0 
           [10] = {Point64} 90000,50000,0 
           [11] = {Point64} 90000,0,0 
           [12] = {Point64} 0,0,0 
           */
        
        Paths intRes2c = GeoWrangler.close(intRes2);
        
        /* Expected output
         intRes2c = {List<List<Point64>>} Count = 2
          [0] = {List<Point64>} Count = 6
           [0] = {Point64} 0,0,0 
           [1] = {Point64} 0,200000,0 
           [2] = {Point64} 200000,200000,0 
           [3] = {Point64} 200000,0,0 
           [4] = {Point64} 110000,0,0 
           [5] = {Point64} 0,0,0 
          [1] = {List<Point64>} Count = 9
           [0] = {Point64} 90000,0,0 
           [1] = {Point64} 110000,0,0 
           [2] = {Point64} 110000,50000,0 
           [3] = {Point64} 150000,50000,0 
           [4] = {Point64} 150000,150000,0 
           [5] = {Point64} 50000,150000,0 
           [6] = {Point64} 50000,50000,0 
           [7] = {Point64} 90000,50000,0 
           [8] = {Point64} 90000,0,0 
           */
        
        Paths intRes2c_kH = GeoWrangler.makeKeyHole(intRes2c, reverseEval:false, biDirectionalEval:true);

        /* Expected output
         intRes2c_kH = {List<List<Point64>>} Count = 1
          [0] = {List<Point64>} Count = 13
           [0] = {Point64} 0,0,0 
           [1] = {Point64} 0,200000,0 
           [2] = {Point64} 200000,200000,0 
           [3] = {Point64} 200000,0,0 
           [4] = {Point64} 110000,0,0 
           [5] = {Point64} 110000,50000,0 
           [6] = {Point64} 150000,50000,0 
           [7] = {Point64} 150000,150000,0 
           [8] = {Point64} 50000,150000,0 
           [9] = {Point64} 50000,50000,0 
           [10] = {Point64} 90000,50000,0 
           [11] = {Point64} 90000,0,0 
           [12] = {Point64} 0,0,0 
           */
        
        Paths intRes2P = new();
        c.Execute(ClipType.Intersection, FillRule.Positive, intRes2P);
        
        /* Expected output
         No geometry in intRes2P
         */
        
        // No keyholes as no geometry.
        Paths intRes2P_kH = GeoWrangler.makeKeyHole(intRes2P, reverseEval:false, biDirectionalEval:true);
        Paths intRes2Pc = GeoWrangler.close(intRes2P);
        Paths intRes2Pc_kH = GeoWrangler.makeKeyHole(intRes2Pc, reverseEval:false, biDirectionalEval:true);

        // seems good - get keyhole
        Paths intRes2NZ = new();
        c.Execute(ClipType.Intersection, FillRule.NonZero, intRes2NZ);
        
        /* Expected output
         intRes2NZ = {List<List<Point64>>} Count = 2
          [0] = {List<Point64>} Count = 6
           [0] = {Point64} 0,0,0 
           [1] = {Point64} 0,200000,0 
           [2] = {Point64} 200000,200000,0 
           [3] = {Point64} 200000,0,0 
           [4] = {Point64} 110000,0,0 
           [5] = {Point64} 0,0,0 
          [1] = {List<Point64>} Count = 6
           [0] = {Point64} 90000,50000,0 
           [1] = {Point64} 150000,50000,0 
           [2] = {Point64} 150000,150000,0 
           [3] = {Point64} 50000,150000,0 
           [4] = {Point64} 50000,50000,0 
           [5] = {Point64} 90000,50000,0 
           */
        
        Paths intRes2NZ_kH = GeoWrangler.makeKeyHole(intRes2NZ, reverseEval:false, biDirectionalEval:true);
        
        /* Expected output
         intRes2NZ_kH = {List<List<Point64>>} Count = 1
          [0] = {List<Point64>} Count = 13
           [0] = {Point64} 200000,0,0 
           [1] = {Point64} 200000,49500,0 
           [2] = {Point64} 149500,49500,0 
           [3] = {Point64} 149500,50000,0 
           [4] = {Point64} 50000,50000,0 
           [5] = {Point64} 50000,150000,0 
           [6] = {Point64} 150000,150000,0 
           [7] = {Point64} 150000,50500,0 
           [8] = {Point64} 200000,50500,0 
           [9] = {Point64} 200000,200000,0 
           [10] = {Point64} 0,200000,0 
           [11] = {Point64} 0,0,0 
           [12] = {Point64} 200000,0,0 
           */
        
        Paths intRes2NZc = GeoWrangler.close(intRes2NZ);
        
        /* Expected output
         intRes2NZc = {List<List<Point64>>} Count = 2
          [0] = {List<Point64>} Count = 6
           [0] = {Point64} 0,0,0 
           [1] = {Point64} 0,200000,0 
           [2] = {Point64} 200000,200000,0 
           [3] = {Point64} 200000,0,0 
           [4] = {Point64} 110000,0,0 
           [5] = {Point64} 0,0,0 
          [1] = {List<Point64>} Count = 6
           [0] = {Point64} 90000,50000,0 
           [1] = {Point64} 150000,50000,0 
           [2] = {Point64} 150000,150000,0 
           [3] = {Point64} 50000,150000,0 
           [4] = {Point64} 50000,50000,0 
           [5] = {Point64} 90000,50000,0 
           */
        
        Paths intRes2NZc_kH = GeoWrangler.makeKeyHole(intRes2NZc, reverseEval:false, biDirectionalEval:true);
        
         /* Expected output
          intRes2NZc_kH = {List<List<Point64>>} Count = 1
           [0] = {List<Point64>} Count = 13
            [0] = {Point64} 200000,0,0 
            [1] = {Point64} 200000,49500,0 
            [2] = {Point64} 149500,49500,0 
            [3] = {Point64} 149500,50000,0 
            [4] = {Point64} 50000,50000,0 
            [5] = {Point64} 50000,150000,0 
            [6] = {Point64} 150000,150000,0 
            [7] = {Point64} 150000,50500,0 
            [8] = {Point64} 200000,50500,0 
            [9] = {Point64} 200000,200000,0 
            [10] = {Point64} 0,200000,0 
            [11] = {Point64} 0,0,0 
            [12] = {Point64} 200000,0,0 
          */
        
    }

    private static void selfOverlapTest_reversed()
    {
        Path outer = new()
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
        Paths dSource = new() {outer};

        /* Expected
         dSource = {List<List<Point64>>} Count = 1
          [0] = {List<Point64>} Count = 13
           [0] = {Point64} 0,0,0 
           [1] = {Point64} 0,200000,0 
           [2] = {Point64} 200000,200000,0 
           [3] = {Point64} 200000,0,0 
           [4] = {Point64} 90000,0,0 
           [5] = {Point64} 90000,50000,0 
           [6] = {Point64} 150000,50000,0 
           [7] = {Point64} 150000,150000,0 
           [8] = {Point64} 50000,150000,0 
           [9] = {Point64} 50000,50000,0 
           [10] = {Point64} 110000,50000,0 
           [11] = {Point64} 110000,0,0 
           [12] = {Point64} 0,0,0 
           */
        
        Paths decomp = GeoWrangler.decompose(dSource);
        
        /* Expected output
         decomp = {List<List<Point64>>} Count = 2
          [0] = {List<Point64>} Count = 6
           [0] = {Point64} 0,0,0 
           [1] = {Point64} 0,200000,0 
           [2] = {Point64} 200000,200000,0 
           [3] = {Point64} 200000,0,0 
           [4] = {Point64} 110000,0,0 
           [5] = {Point64} 0,0,0 
          [1] = {List<Point64>} Count = 6
           [0] = {Point64} 50000,50000,0 
           [1] = {Point64} 90000,50000,0 
           [2] = {Point64} 150000,50000,0 
           [3] = {Point64} 150000,150000,0 
           [4] = {Point64} 50000,150000,0 
           [5] = {Point64} 50000,50000,0 
           */
        
        Paths kHD = GeoWrangler.makeKeyHole(dSource, reverseEval:false, biDirectionalEval:true);
        
        /* Expected output
         kHD = {List<List<Point64>>} Count = 1
          [0] = {List<Point64>} Count = 13
           [0] = {Point64} 200000,0,0 
           [1] = {Point64} 200000,49500,0 
           [2] = {Point64} 149500,49500,0 
           [3] = {Point64} 149500,50000,0 
           [4] = {Point64} 50000,50000,0 
           [5] = {Point64} 50000,150000,0 
           [6] = {Point64} 150000,150000,0 
           [7] = {Point64} 150000,50500,0 
           [8] = {Point64} 200000,50500,0 
           [9] = {Point64} 200000,200000,0 
           [10] = {Point64} 0,200000,0 
           [11] = {Point64} 0,0,0 
           [12] = {Point64} 200000,0,0 
           */

        // keyholer test
        Paths kHSource = new() {outer};
        
        /* Expected
         kHSource = {List<List<Point64>>} Count = 1
          [0] = {List<Point64>} Count = 13
           [0] = {Point64} 0,0,0 
           [1] = {Point64} 0,200000,0 
           [2] = {Point64} 200000,200000,0 
           [3] = {Point64} 200000,0,0 
           [4] = {Point64} 90000,0,0 
           [5] = {Point64} 90000,50000,0 
           [6] = {Point64} 150000,50000,0 
           [7] = {Point64} 150000,150000,0 
           [8] = {Point64} 50000,150000,0 
           [9] = {Point64} 50000,50000,0 
           [10] = {Point64} 110000,50000,0 
           [11] = {Point64} 110000,0,0 
           [12] = {Point64} 0,0,0 
           */
        
        Paths kH = GeoWrangler.makeKeyHole(kHSource, reverseEval:false, biDirectionalEval:true);

        /* Expected output
         kH = {List<List<Point64>>} Count = 1
          [0] = {List<Point64>} Count = 13
           [0] = {Point64} 200000,0,0 
           [1] = {Point64} 200000,49500,0 
           [2] = {Point64} 149500,49500,0 
           [3] = {Point64} 149500,50000,0 
           [4] = {Point64} 50000,50000,0 
           [5] = {Point64} 50000,150000,0 
           [6] = {Point64} 150000,150000,0 
           [7] = {Point64} 150000,50500,0 
           [8] = {Point64} 200000,50500,0 
           [9] = {Point64} 200000,200000,0 
           [10] = {Point64} 0,200000,0 
           [11] = {Point64} 0,0,0 
           [12] = {Point64} 200000,0,0 
           */
        
        Clipper64 c = new();
        c.AddSubject(outer);

        Paths unionRes = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, unionRes);

        /* Expected output
         unionRes = {List<List<Point64>>} Count = 1
          [0] = {List<Point64>} Count = 13
           [0] = {Point64} 90000,0,0 
           [1] = {Point64} 0,0,0 
           [2] = {Point64} 0,200000,0 
           [3] = {Point64} 200000,200000,0 
           [4] = {Point64} 200000,0,0 
           [5] = {Point64} 110000,0,0 
           [6] = {Point64} 110000,50000,0 
           [7] = {Point64} 150000,50000,0 
           [8] = {Point64} 150000,150000,0 
           [9] = {Point64} 50000,150000,0 
           [10] = {Point64} 50000,50000,0 
           [11] = {Point64} 90000,50000,0 
           [12] = {Point64} 90000,0,0 
           */
        
        Paths unionRes_kH = GeoWrangler.makeKeyHole(unionRes, reverseEval:false, biDirectionalEval:true);
        
        /* Expected output
         unionRes_kH = {List<List<Point64>>} Count = 1
          [0] = {List<Point64>} Count = 13
           [0] = {Point64} 0,0,0 
           [1] = {Point64} 0,200000,0 
           [2] = {Point64} 200000,200000,0 
           [3] = {Point64} 200000,0,0 
           [4] = {Point64} 110000,0,0 
           [5] = {Point64} 110000,50000,0 
           [6] = {Point64} 150000,50000,0 
           [7] = {Point64} 150000,150000,0 
           [8] = {Point64} 50000,150000,0 
           [9] = {Point64} 50000,50000,0 
           [10] = {Point64} 90000,50000,0 
           [11] = {Point64} 90000,0,0 
           [12] = {Point64} 0,0,0 
           */
        
        Paths unionResc = GeoWrangler.close(unionRes);
        
        /* Expected output
         unionResc = {List<List<Point64>>} Count = 1
          [0] = {List<Point64>} Count = 13
           [0] = {Point64} 90000,0,0 
           [1] = {Point64} 0,0,0 
           [2] = {Point64} 0,200000,0 
           [3] = {Point64} 200000,200000,0 
           [4] = {Point64} 200000,0,0 
           [5] = {Point64} 110000,0,0 
           [6] = {Point64} 110000,50000,0 
           [7] = {Point64} 150000,50000,0 
           [8] = {Point64} 150000,150000,0 
           [9] = {Point64} 50000,150000,0 
           [10] = {Point64} 50000,50000,0 
           [11] = {Point64} 90000,50000,0 
           [12] = {Point64} 90000,0,0 
           */
        
        Paths unionResc_kH = GeoWrangler.makeKeyHole(unionResc, reverseEval:false, biDirectionalEval:true);

        /* Expected output
         unionResc_kH = {List<List<Point64>>} Count = 1
          [0] = {List<Point64>} Count = 13
           [0] = {Point64} 0,0,0 
           [1] = {Point64} 0,200000,0 
           [2] = {Point64} 200000,200000,0 
           [3] = {Point64} 200000,0,0 
           [4] = {Point64} 110000,0,0 
           [5] = {Point64} 110000,50000,0 
           [6] = {Point64} 150000,50000,0 
           [7] = {Point64} 150000,150000,0 
           [8] = {Point64} 50000,150000,0 
           [9] = {Point64} 50000,50000,0 
           [10] = {Point64} 90000,50000,0 
           [11] = {Point64} 90000,0,0 
           [12] = {Point64} 0,0,0 
           */
        
        Paths unionResP = new();
        c.Execute(ClipType.Union, FillRule.Positive, unionResP);

        /* Expected output
           unionResP = {List<List<Point64>>} Count = 2
            [0] = {List<Point64>} Count = 6
             [0] = {Point64} 0,0,0 
             [1] = {Point64} 0,200000,0 
             [2] = {Point64} 200000,200000,0 
             [3] = {Point64} 200000,0,0 
             [4] = {Point64} 110000,0,0 
             [5] = {Point64} 0,0,0 
            [1] = {List<Point64>} Count = 6
             [0] = {Point64} 90000,50000,0 
             [1] = {Point64} 150000,50000,0 
             [2] = {Point64} 150000,150000,0 
             [3] = {Point64} 50000,150000,0 
             [4] = {Point64} 50000,50000,0 
             [5] = {Point64} 90000,50000,0 
           */
        
        Paths unionResP_kH = GeoWrangler.makeKeyHole(unionResP, reverseEval:false, biDirectionalEval:true);
        
        /* Expected output
         unionResP_kH = {List<List<Point64>>} Count = 1
          [0] = {List<Point64>} Count = 13
           [0] = {Point64} 200000,0,0 
           [1] = {Point64} 200000,49500,0 
           [2] = {Point64} 149500,49500,0 
           [3] = {Point64} 149500,50000,0 
           [4] = {Point64} 50000,50000,0 
           [5] = {Point64} 50000,150000,0 
           [6] = {Point64} 150000,150000,0 
           [7] = {Point64} 150000,50500,0 
           [8] = {Point64} 200000,50500,0 
           [9] = {Point64} 200000,200000,0 
           [10] = {Point64} 0,200000,0 
           [11] = {Point64} 0,0,0 
           [12] = {Point64} 200000,0,0 
         */
        
        Paths unionResPc = GeoWrangler.close(unionResP);
        
        /* Expected output
          unionResPc = {List<List<Point64>>} Count = 2
           [0] = {List<Point64>} Count = 6
            [0] = {Point64} 0,0,0 
            [1] = {Point64} 0,200000,0 
            [2] = {Point64} 200000,200000,0 
            [3] = {Point64} 200000,0,0 
            [4] = {Point64} 110000,0,0 
            [5] = {Point64} 0,0,0 
           [1] = {List<Point64>} Count = 6
            [0] = {Point64} 90000,50000,0 
            [1] = {Point64} 150000,50000,0 
            [2] = {Point64} 150000,150000,0 
            [3] = {Point64} 50000,150000,0 
            [4] = {Point64} 50000,50000,0 
            [5] = {Point64} 90000,50000,0 
           */
        
        Paths unionResPc_kH = GeoWrangler.makeKeyHole(unionResPc, reverseEval:false, biDirectionalEval:true);

        /* Expected output
         unionResPc_kH = {List<List<Point64>>} Count = 1
          [0] = {List<Point64>} Count = 13
           [0] = {Point64} 200000,0,0 
           [1] = {Point64} 200000,49500,0 
           [2] = {Point64} 149500,49500,0 
           [3] = {Point64} 149500,50000,0 
           [4] = {Point64} 50000,50000,0 
           [5] = {Point64} 50000,150000,0 
           [6] = {Point64} 150000,150000,0 
           [7] = {Point64} 150000,50500,0 
           [8] = {Point64} 200000,50500,0 
           [9] = {Point64} 200000,200000,0 
           [10] = {Point64} 0,200000,0 
           [11] = {Point64} 0,0,0 
           [12] = {Point64} 200000,0,0 
           */
        
        // seems good - get keyhole
        Paths unionResNZ = new();
        c.Execute(ClipType.Union, FillRule.NonZero, unionResNZ);
        
        /* Expected output
           unionResNZ = {List<List<Point64>>} Count = 2
            [0] = {List<Point64>} Count = 6
             [0] = {Point64} 0,0,0 
             [1] = {Point64} 0,200000,0 
             [2] = {Point64} 200000,200000,0 
             [3] = {Point64} 200000,0,0 
             [4] = {Point64} 110000,0,0 
             [5] = {Point64} 0,0,0 
            [1] = {List<Point64>} Count = 6
             [0] = {Point64} 90000,50000,0 
             [1] = {Point64} 150000,50000,0 
             [2] = {Point64} 150000,150000,0 
             [3] = {Point64} 50000,150000,0 
             [4] = {Point64} 50000,50000,0 
             [5] = {Point64} 90000,50000,0 
           */
        
        Paths unionResNZ_kH = GeoWrangler.makeKeyHole(unionResNZ, reverseEval:false, biDirectionalEval:true);
        
        /* Expected output
           unionResNZ_kH = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 13
             [0] = {Point64} 200000,0,0 
             [1] = {Point64} 200000,49500,0 
             [2] = {Point64} 149500,49500,0 
             [3] = {Point64} 149500,50000,0 
             [4] = {Point64} 50000,50000,0 
             [5] = {Point64} 50000,150000,0 
             [6] = {Point64} 150000,150000,0 
             [7] = {Point64} 150000,50500,0 
             [8] = {Point64} 200000,50500,0 
             [9] = {Point64} 200000,200000,0 
             [10] = {Point64} 0,200000,0 
             [11] = {Point64} 0,0,0 
             [12] = {Point64} 200000,0,0 
           */
        
        Paths unionResNZc = GeoWrangler.close(unionResNZ);
        
        /* Expected output
           unionResNZc = {List<List<Point64>>} Count = 2
            [0] = {List<Point64>} Count = 6
             [0] = {Point64} 0,0,0 
             [1] = {Point64} 0,200000,0 
             [2] = {Point64} 200000,200000,0 
             [3] = {Point64} 200000,0,0 
             [4] = {Point64} 110000,0,0 
             [5] = {Point64} 0,0,0 
            [1] = {List<Point64>} Count = 6
             [0] = {Point64} 90000,50000,0 
             [1] = {Point64} 150000,50000,0 
             [2] = {Point64} 150000,150000,0 
             [3] = {Point64} 50000,150000,0 
             [4] = {Point64} 50000,50000,0 
             [5] = {Point64} 90000,50000,0 
           */
        
        Paths unionResNZc_kH = GeoWrangler.makeKeyHole(unionResNZc, reverseEval:false, biDirectionalEval:true);

        /* Expected output
           unionResNZc_kH = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 13
             [0] = {Point64} 200000,0,0 
             [1] = {Point64} 200000,49500,0 
             [2] = {Point64} 149500,49500,0 
             [3] = {Point64} 149500,50000,0 
             [4] = {Point64} 50000,50000,0 
             [5] = {Point64} 50000,150000,0 
             [6] = {Point64} 150000,150000,0 
             [7] = {Point64} 150000,50500,0 
             [8] = {Point64} 200000,50500,0 
             [9] = {Point64} 200000,200000,0 
             [10] = {Point64} 0,200000,0 
             [11] = {Point64} 0,0,0 
             [12] = {Point64} 200000,0,0 
           */
        
        // no good - overlap region is a gap.
        Paths simplifyRes = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, simplifyRes);
        simplifyRes = GeoWrangler.stripColinear(simplifyRes);
        
        /* Expected output
           simplifyRes = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 13
             [0] = {Point64} 90000,0,0 
             [1] = {Point64} 0,0,0 
             [2] = {Point64} 0,200000,0 
             [3] = {Point64} 200000,200000,0 
             [4] = {Point64} 200000,0,0 
             [5] = {Point64} 110000,0,0 
             [6] = {Point64} 110000,50000,0 
             [7] = {Point64} 150000,50000,0 
             [8] = {Point64} 150000,150000,0 
             [9] = {Point64} 50000,150000,0 
             [10] = {Point64} 50000,50000,0 
             [11] = {Point64} 90000,50000,0 
             [12] = {Point64} 90000,0,0 
           */
        
        Paths simplifyRes_kH = GeoWrangler.makeKeyHole(simplifyRes, reverseEval:false, biDirectionalEval:true);
        
        /* Expected output
           simplifyRes_kH = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 13
             [0] = {Point64} 0,0,0 
             [1] = {Point64} 0,200000,0 
             [2] = {Point64} 200000,200000,0 
             [3] = {Point64} 200000,0,0 
             [4] = {Point64} 110000,0,0 
             [5] = {Point64} 110000,50000,0 
             [6] = {Point64} 150000,50000,0 
             [7] = {Point64} 150000,150000,0 
             [8] = {Point64} 50000,150000,0 
             [9] = {Point64} 50000,50000,0 
             [10] = {Point64} 90000,50000,0 
             [11] = {Point64} 90000,0,0 
             [12] = {Point64} 0,0,0 
           */
        
        Paths simplifyResc = GeoWrangler.close(simplifyRes);
        
        /* Expected output
           simplifyResc = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 13
             [0] = {Point64} 90000,0,0 
             [1] = {Point64} 0,0,0 
             [2] = {Point64} 0,200000,0 
             [3] = {Point64} 200000,200000,0 
             [4] = {Point64} 200000,0,0 
             [5] = {Point64} 110000,0,0 
             [6] = {Point64} 110000,50000,0 
             [7] = {Point64} 150000,50000,0 
             [8] = {Point64} 150000,150000,0 
             [9] = {Point64} 50000,150000,0 
             [10] = {Point64} 50000,50000,0 
             [11] = {Point64} 90000,50000,0 
             [12] = {Point64} 90000,0,0 
           */
        
        Paths simplifyResc_kH = GeoWrangler.makeKeyHole(simplifyResc, reverseEval:false, biDirectionalEval:true);
        
        /* Expected output
           simplifyResc_kH = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 13
             [0] = {Point64} 0,0,0 
             [1] = {Point64} 0,200000,0 
             [2] = {Point64} 200000,200000,0 
             [3] = {Point64} 200000,0,0 
             [4] = {Point64} 110000,0,0 
             [5] = {Point64} 110000,50000,0 
             [6] = {Point64} 150000,50000,0 
             [7] = {Point64} 150000,150000,0 
             [8] = {Point64} 50000,150000,0 
             [9] = {Point64} 50000,50000,0 
             [10] = {Point64} 90000,50000,0 
             [11] = {Point64} 90000,0,0 
             [12] = {Point64} 0,0,0 
           */

        Paths simplifyRes2 = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, simplifyRes2);
        
        /* Expected output
           simplifyRes2 = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 13
             [0] = {Point64} 90000,0,0 
             [1] = {Point64} 0,0,0 
             [2] = {Point64} 0,200000,0 
             [3] = {Point64} 200000,200000,0 
             [4] = {Point64} 200000,0,0 
             [5] = {Point64} 110000,0,0 
             [6] = {Point64} 110000,50000,0 
             [7] = {Point64} 150000,50000,0 
             [8] = {Point64} 150000,150000,0 
             [9] = {Point64} 50000,150000,0 
             [10] = {Point64} 50000,50000,0 
             [11] = {Point64} 90000,50000,0 
             [12] = {Point64} 90000,0,0 
           */
        
        Paths simplifyRes2_kH = GeoWrangler.makeKeyHole(simplifyRes2, reverseEval:false, biDirectionalEval:true);
        
        /* Expected output
           simplifyRes2_kH = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 13
             [0] = {Point64} 0,0,0 
             [1] = {Point64} 0,200000,0 
             [2] = {Point64} 200000,200000,0 
             [3] = {Point64} 200000,0,0 
             [4] = {Point64} 110000,0,0 
             [5] = {Point64} 110000,50000,0 
             [6] = {Point64} 150000,50000,0 
             [7] = {Point64} 150000,150000,0 
             [8] = {Point64} 50000,150000,0 
             [9] = {Point64} 50000,50000,0 
             [10] = {Point64} 90000,50000,0 
             [11] = {Point64} 90000,0,0 
             [12] = {Point64} 0,0,0 
           */
        
        Paths simplifyRes2c = GeoWrangler.close(simplifyRes2);
        
        /* Expected output
           simplifyRes2c = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 13
             [0] = {Point64} 90000,0,0 
             [1] = {Point64} 0,0,0 
             [2] = {Point64} 0,200000,0 
             [3] = {Point64} 200000,200000,0 
             [4] = {Point64} 200000,0,0 
             [5] = {Point64} 110000,0,0 
             [6] = {Point64} 110000,50000,0 
             [7] = {Point64} 150000,50000,0 
             [8] = {Point64} 150000,150000,0 
             [9] = {Point64} 50000,150000,0 
             [10] = {Point64} 50000,50000,0 
             [11] = {Point64} 90000,50000,0 
             [12] = {Point64} 90000,0,0 
           */
        
        Paths simplifyRes2c_kH = GeoWrangler.makeKeyHole(simplifyRes2c, reverseEval:false, biDirectionalEval:true);
        
        /* Expected output
           simplifyRes2c_kH = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 13
             [0] = {Point64} 0,0,0 
             [1] = {Point64} 0,200000,0 
             [2] = {Point64} 200000,200000,0 
             [3] = {Point64} 200000,0,0 
             [4] = {Point64} 110000,0,0 
             [5] = {Point64} 110000,50000,0 
             [6] = {Point64} 150000,50000,0 
             [7] = {Point64} 150000,150000,0 
             [8] = {Point64} 50000,150000,0 
             [9] = {Point64} 50000,50000,0 
             [10] = {Point64} 90000,50000,0 
             [11] = {Point64} 90000,0,0 
             [12] = {Point64} 0,0,0 
           */

        // no good - no result
        Paths intRes = new();
        c.Execute(ClipType.Intersection, FillRule.EvenOdd, intRes);
        Paths intRes_kH = GeoWrangler.makeKeyHole(intRes, reverseEval:false, biDirectionalEval:true);
        Paths intResc = GeoWrangler.close(intRes);
        Paths intResc_kH = GeoWrangler.makeKeyHole(intResc, reverseEval:false, biDirectionalEval:true);

        // no good - no result
        Paths intResP = new();
        c.Execute(ClipType.Intersection, FillRule.Positive, intResP);
        Paths intResP_kH = GeoWrangler.makeKeyHole(intResP, reverseEval:false, biDirectionalEval:true);
        Paths intResPc = GeoWrangler.close(intResP);
        Paths intResPc_kH = GeoWrangler.makeKeyHole(intResPc, reverseEval:false, biDirectionalEval:true);

        // no good - no result
        Paths intResNZ = new();
        c.Execute(ClipType.Intersection, FillRule.NonZero, intResNZ);
        Paths intResNZ_kH = GeoWrangler.makeKeyHole(intResNZ, reverseEval:false, biDirectionalEval:true);
        Paths intResNZc = GeoWrangler.close(intResNZ);
        Paths intResNZc_kH = GeoWrangler.makeKeyHole(intResNZc, reverseEval:false, biDirectionalEval:true);

        Rect64 bounds = Clipper.GetBounds(new Paths { outer });
        Path bb = new()
        {
            new(bounds.left, bounds.bottom),
            new(bounds.left, bounds.top),
            new(bounds.right, bounds.top),
            new(bounds.right, bounds.bottom)
        };

        c.Clear();
        c.AddSubject(bb);
        c.AddClip(outer);

        Paths intRes2 = new();
        c.Execute(ClipType.Intersection, FillRule.EvenOdd, intRes2);
        
        /* Expected output
           intRes2 = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 13
             [0] = {Point64} 90000,0,0 
             [1] = {Point64} 0,0,0 
             [2] = {Point64} 0,200000,0 
             [3] = {Point64} 200000,200000,0 
             [4] = {Point64} 200000,0,0 
             [5] = {Point64} 110000,0,0 
             [6] = {Point64} 110000,50000,0 
             [7] = {Point64} 150000,50000,0 
             [8] = {Point64} 150000,150000,0 
             [9] = {Point64} 50000,150000,0 
             [10] = {Point64} 50000,50000,0 
             [11] = {Point64} 90000,50000,0 
             [12] = {Point64} 90000,0,0 
           */
        
        Paths intRes2_kH = GeoWrangler.makeKeyHole(intRes2, reverseEval:false, biDirectionalEval:true);
        
        /* Expected output
           intRes2_kH = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 13
             [0] = {Point64} 0,0,0 
             [1] = {Point64} 0,200000,0 
             [2] = {Point64} 200000,200000,0 
             [3] = {Point64} 200000,0,0 
             [4] = {Point64} 110000,0,0 
             [5] = {Point64} 110000,50000,0 
             [6] = {Point64} 150000,50000,0 
             [7] = {Point64} 150000,150000,0 
             [8] = {Point64} 50000,150000,0 
             [9] = {Point64} 50000,50000,0 
             [10] = {Point64} 90000,50000,0 
             [11] = {Point64} 90000,0,0 
             [12] = {Point64} 0,0,0 
           */
        
        Paths intRes2c = GeoWrangler.close(intRes2);
        
        /* Expected output
           intRes2c = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 13
             [0] = {Point64} 90000,0,0 
             [1] = {Point64} 0,0,0 
             [2] = {Point64} 0,200000,0 
             [3] = {Point64} 200000,200000,0 
             [4] = {Point64} 200000,0,0 
             [5] = {Point64} 110000,0,0 
             [6] = {Point64} 110000,50000,0 
             [7] = {Point64} 150000,50000,0 
             [8] = {Point64} 150000,150000,0 
             [9] = {Point64} 50000,150000,0 
             [10] = {Point64} 50000,50000,0 
             [11] = {Point64} 90000,50000,0 
             [12] = {Point64} 90000,0,0 
           */
        
        Paths intRes2c_kH = GeoWrangler.makeKeyHole(intRes2c, reverseEval:false, biDirectionalEval:true);
        
        /* Expected output
           intRes2c_kH = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 13
             [0] = {Point64} 0,0,0 
             [1] = {Point64} 0,200000,0 
             [2] = {Point64} 200000,200000,0 
             [3] = {Point64} 200000,0,0 
             [4] = {Point64} 110000,0,0 
             [5] = {Point64} 110000,50000,0 
             [6] = {Point64} 150000,50000,0 
             [7] = {Point64} 150000,150000,0 
             [8] = {Point64} 50000,150000,0 
             [9] = {Point64} 50000,50000,0 
             [10] = {Point64} 90000,50000,0 
             [11] = {Point64} 90000,0,0 
             [12] = {Point64} 0,0,0 
           */

        // no good - no result
        Paths intRes2P = new();
        c.Execute(ClipType.Intersection, FillRule.Positive, intRes2P);
        
        /* Expected output
         No geometry in intRes2P
         */
        
        // No results - no geometry
        Paths intRes2P_kH = GeoWrangler.makeKeyHole(intRes2P, reverseEval:false, biDirectionalEval:true);
        Paths intRes2Pc = GeoWrangler.close(intRes2P);
        Paths intRes2Pc_kH = GeoWrangler.makeKeyHole(intRes2Pc, reverseEval:false, biDirectionalEval:true);

        // seems good - get keyhole
        Paths intRes2NZ = new();
        c.Execute(ClipType.Intersection, FillRule.NonZero, intRes2NZ);
        
        /* Expected output
           intRes2NZ = {List<List<Point64>>} Count = 2
            [0] = {List<Point64>} Count = 6
             [0] = {Point64} 0,0,0 
             [1] = {Point64} 0,200000,0 
             [2] = {Point64} 200000,200000,0 
             [3] = {Point64} 200000,0,0 
             [4] = {Point64} 110000,0,0 
             [5] = {Point64} 0,0,0 
            [1] = {List<Point64>} Count = 6
             [0] = {Point64} 90000,50000,0 
             [1] = {Point64} 150000,50000,0 
             [2] = {Point64} 150000,150000,0 
             [3] = {Point64} 50000,150000,0 
             [4] = {Point64} 50000,50000,0 
             [5] = {Point64} 90000,50000,0 
           */
        
        Paths intRes2NZ_kH = GeoWrangler.makeKeyHole(intRes2NZ, reverseEval:false, biDirectionalEval:true);
        
        /* Expected output
           intRes2NZ_kH = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 13
             [0] = {Point64} 200000,0,0 
             [1] = {Point64} 200000,49500,0 
             [2] = {Point64} 149500,49500,0 
             [3] = {Point64} 149500,50000,0 
             [4] = {Point64} 50000,50000,0 
             [5] = {Point64} 50000,150000,0 
             [6] = {Point64} 150000,150000,0 
             [7] = {Point64} 150000,50500,0 
             [8] = {Point64} 200000,50500,0 
             [9] = {Point64} 200000,200000,0 
             [10] = {Point64} 0,200000,0 
             [11] = {Point64} 0,0,0 
             [12] = {Point64} 200000,0,0 
           */
        
        Paths intRes2NZc = GeoWrangler.close(intRes2NZ);
        
        /* Expected output
           intRes2NZc = {List<List<Point64>>} Count = 2
            [0] = {List<Point64>} Count = 6
             [0] = {Point64} 0,0,0 
             [1] = {Point64} 0,200000,0 
             [2] = {Point64} 200000,200000,0 
             [3] = {Point64} 200000,0,0 
             [4] = {Point64} 110000,0,0 
             [5] = {Point64} 0,0,0 
            [1] = {List<Point64>} Count = 6
             [0] = {Point64} 90000,50000,0 
             [1] = {Point64} 150000,50000,0 
             [2] = {Point64} 150000,150000,0 
             [3] = {Point64} 50000,150000,0 
             [4] = {Point64} 50000,50000,0 
             [5] = {Point64} 90000,50000,0 
           */
        
        Paths intRes2NZc_kH = GeoWrangler.makeKeyHole(intRes2NZc, reverseEval:false, biDirectionalEval:true);
        
        /* Expected output
           intRes2NZc_kH = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 13
             [0] = {Point64} 200000,0,0 
             [1] = {Point64} 200000,49500,0 
             [2] = {Point64} 149500,49500,0 
             [3] = {Point64} 149500,50000,0 
             [4] = {Point64} 50000,50000,0 
             [5] = {Point64} 50000,150000,0 
             [6] = {Point64} 150000,150000,0 
             [7] = {Point64} 150000,50500,0 
             [8] = {Point64} 200000,50500,0 
             [9] = {Point64} 200000,200000,0 
             [10] = {Point64} 0,200000,0 
             [11] = {Point64} 0,0,0 
             [12] = {Point64} 200000,0,0 
           */
        
    }

    private static void comboTest()
    {
        // Manually create the sliver.
        Path outer = new()
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

        Path inner1 = new()
        {
            new(100000, -100000),
            new(100000, 100000),
            new(200000, 100000),
            new(200000, -100000),
            new(100000, -100000)
        };

        // Segment the paths to match real-world case.
        /*
        Fragmenter f = new(10000);
        Path outer_f = f.fragmentPath(outer);

        Path inner1_f = f.fragmentPath(inner1);
        */
        Paths kHSource = new()
        {
            outer,
            inner1
        };

        /* Expected
           kHSource = {List<List<Point64>>} Count = 2
            [0] = {List<Point64>} Count = 9
             [0] = {Point64} -200000,-200000,0 
             [1] = {Point64} 300000,-200000,0 
             [2] = {Point64} 300000,200000,0 
             [3] = {Point64} -200000,200000,0 
             [4] = {Point64} -200000,-99900,0 
             [5] = {Point64} -300000,-99900,0 
             [6] = {Point64} -300000,-100100,0 
             [7] = {Point64} -200000,-100100,0 
             [8] = {Point64} -200000,-200000,0 
            [1] = {List<Point64>} Count = 5
             [0] = {Point64} 100000,-100000,0 
             [1] = {Point64} 100000,100000,0 
             [2] = {Point64} 200000,100000,0 
             [3] = {Point64} 200000,-100000,0 
             [4] = {Point64} 100000,-100000,0 
           */
        
        Paths kH = GeoWrangler.makeKeyHole(kHSource, reverseEval:false, biDirectionalEval:true);

        /* Expected output
           kH = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 17
             [0] = {Point64} -300000,-100100,0 
             [1] = {Point64} -200000,-100100,0 
             [2] = {Point64} -200000,-200000,0 
             [3] = {Point64} 300000,-200000,0 
             [4] = {Point64} 300000,-100500,0 
             [5] = {Point64} 199500,-100500,0 
             [6] = {Point64} 199500,-100000,0 
             [7] = {Point64} 100000,-100000,0 
             [8] = {Point64} 100000,100000,0 
             [9] = {Point64} 200000,100000,0 
             [10] = {Point64} 200000,-99500,0 
             [11] = {Point64} 300000,-99500,0 
             [12] = {Point64} 300000,200000,0 
             [13] = {Point64} -200000,200000,0 
             [14] = {Point64} -200000,-99900,0 
             [15] = {Point64} -300000,-99900,0 
             [16] = {Point64} -300000,-100100,0 
           */
        
        // Gap removal test
        Paths gR = GeoWrangler.gapRemoval(kH, 100);

        /* Expected output
           gR = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 17
             [0] = {Point64} -300000,-100100,0 
             [1] = {Point64} -200000,-100100,0 
             [2] = {Point64} -200000,-200000,0 
             [3] = {Point64} 300000,-200000,0 
             [4] = {Point64} 300000,-100500,0 
             [5] = {Point64} 199500,-100500,0 
             [6] = {Point64} 199500,-100000,0 
             [7] = {Point64} 100000,-100000,0 
             [8] = {Point64} 100000,100000,0 
             [9] = {Point64} 200000,100000,0 
             [10] = {Point64} 200000,-99500,0 
             [11] = {Point64} 300000,-99500,0 
             [12] = {Point64} 300000,200000,0 
             [13] = {Point64} -200000,200000,0 
             [14] = {Point64} -200000,-99900,0 
             [15] = {Point64} -300000,-99900,0 
             [16] = {Point64} -300000,-100100,0 
           */
        
        // Sliver removal test
        Paths sR = GeoWrangler.gapRemoval(kH, -100);
        
        /* Expected output
           sR = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 13
             [0] = {Point64} 300000,-200000,0 
             [1] = {Point64} 300000,-100500,0 
             [2] = {Point64} 199500,-100500,0 
             [3] = {Point64} 199500,-100000,0 
             [4] = {Point64} 100000,-100000,0 
             [5] = {Point64} 100000,100000,0 
             [6] = {Point64} 200000,100000,0 
             [7] = {Point64} 200000,-99500,0 
             [8] = {Point64} 300000,-99500,0 
             [9] = {Point64} 300000,200000,0 
             [10] = {Point64} -200000,200000,0 
             [11] = {Point64} -200000,-200000,0 
             [12] = {Point64} 300000,-200000,0 
           */
        
    }

    private static void simple_islandTest()
    {
        Path outer1 = new()
        {
            new(-200000, -200000),
            new(200000, -200000),
            new(200000, 200000),
            new(-200000, 200000),
            new(-200000, -200000)
        };

        Path inner1 = new()
        {
            new(-100000, -100000),
            new(-100000, 100000),
            new(100000, 100000),
            new(100000, -100000),
            new(-100000, -100000)
        };

        Path outer2 = new()
        {
            new(-200000, 400000),
            new(200000, 400000),
            new(200000, 800000),
            new(-200000, 800000),
            new(-200000, 400000)
        };

        Path inner2 = new()
        {
            new(-100000, 500000),
            new(-100000, 700000),
            new(100000, 700000),
            new(100000, 500000),
            new(-100000, 500000)
        };

        Paths kHSource = new()
        {
            outer1,
            inner1,
            outer2,
            inner2
        };
        
        /* Expected
           kHSource = {List<List<Point64>>} Count = 4
            [0] = {List<Point64>} Count = 5
             [0] = {Point64} -200000,-200000,0 
             [1] = {Point64} 200000,-200000,0 
             [2] = {Point64} 200000,200000,0 
             [3] = {Point64} -200000,200000,0 
             [4] = {Point64} -200000,-200000,0 
            [1] = {List<Point64>} Count = 5
             [0] = {Point64} -100000,-100000,0 
             [1] = {Point64} -100000,100000,0 
             [2] = {Point64} 100000,100000,0 
             [3] = {Point64} 100000,-100000,0 
             [4] = {Point64} -100000,-100000,0 
            [2] = {List<Point64>} Count = 5
             [0] = {Point64} -200000,400000,0 
             [1] = {Point64} 200000,400000,0 
             [2] = {Point64} 200000,800000,0 
             [3] = {Point64} -200000,800000,0 
             [4] = {Point64} -200000,400000,0 
            [3] = {List<Point64>} Count = 5
             [0] = {Point64} -100000,500000,0 
             [1] = {Point64} -100000,700000,0 
             [2] = {Point64} 100000,700000,0 
             [3] = {Point64} 100000,500000,0 
             [4] = {Point64} -100000,500000,0 
           */

        // Generate keyholed geometry
        Paths kH = GeoWrangler.makeKeyHole(kHSource, reverseEval:false, biDirectionalEval:true);

        /* Expected output
           kH = {List<List<Point64>>} Count = 2
            [0] = {List<Point64>} Count = 13
             [0] = {Point64} 200000,400000,0 
             [1] = {Point64} 200000,499500,0 
             [2] = {Point64} 99500,499500,0 
             [3] = {Point64} 99500,500000,0 
             [4] = {Point64} -100000,500000,0 
             [5] = {Point64} -100000,700000,0 
             [6] = {Point64} 100000,700000,0 
             [7] = {Point64} 100000,500500,0 
             [8] = {Point64} 200000,500500,0 
             [9] = {Point64} 200000,800000,0 
             [10] = {Point64} -200000,800000,0 
             [11] = {Point64} -200000,400000,0 
             [12] = {Point64} 200000,400000,0 
            [1] = {List<Point64>} Count = 13
             [0] = {Point64} 200000,-200000,0 
             [1] = {Point64} 200000,-100500,0 
             [2] = {Point64} 99500,-100500,0 
             [3] = {Point64} 99500,-100000,0 
             [4] = {Point64} -100000,-100000,0 
             [5] = {Point64} -100000,100000,0 
             [6] = {Point64} 100000,100000,0 
             [7] = {Point64} 100000,-99500,0 
             [8] = {Point64} 200000,-99500,0 
             [9] = {Point64} 200000,200000,0 
             [10] = {Point64} -200000,200000,0 
             [11] = {Point64} -200000,-200000,0 
             [12] = {Point64} 200000,-200000,0 
           */

        // Gap removal test
        Paths gR = GeoWrangler.gapRemoval(kH, 100);

        /* Expected output
           gR = {List<List<Point64>>} Count = 2
            [0] = {List<Point64>} Count = 13
             [0] = {Point64} 200000,400000,0 
             [1] = {Point64} 200000,499500,0 
             [2] = {Point64} 99500,499500,0 
             [3] = {Point64} 99500,500000,0 
             [4] = {Point64} -100000,500000,0 
             [5] = {Point64} -100000,700000,0 
             [6] = {Point64} 100000,700000,0 
             [7] = {Point64} 100000,500500,0 
             [8] = {Point64} 200000,500500,0 
             [9] = {Point64} 200000,800000,0 
             [10] = {Point64} -200000,800000,0 
             [11] = {Point64} -200000,400000,0 
             [12] = {Point64} 200000,400000,0 
            [1] = {List<Point64>} Count = 13
             [0] = {Point64} 200000,-200000,0 
             [1] = {Point64} 200000,-100500,0 
             [2] = {Point64} 99500,-100500,0 
             [3] = {Point64} 99500,-100000,0 
             [4] = {Point64} -100000,-100000,0 
             [5] = {Point64} -100000,100000,0 
             [6] = {Point64} 100000,100000,0 
             [7] = {Point64} 100000,-99500,0 
             [8] = {Point64} 200000,-99500,0 
             [9] = {Point64} 200000,200000,0 
             [10] = {Point64} -200000,200000,0 
             [11] = {Point64} -200000,-200000,0 
             [12] = {Point64} 200000,-200000,0 
           */
        
        // Generate sliver geometry.
        Paths sL = new();
        Clipper64 c = new();
        c.AddSubject(outer1);
        c.AddSubject(outer2);
        c.AddClip(kH);
        c.Execute(ClipType.Difference, FillRule.EvenOdd, sL);
        
        /* Expected output
           sL = {List<List<Point64>>} Count = 2
            [0] = {List<Point64>} Count = 8
             [0] = {Point64} 99500,499500,0 
             [1] = {Point64} 99500,500000,0 
             [2] = {Point64} -100000,500000,0 
             [3] = {Point64} -100000,700000,0 
             [4] = {Point64} 100000,700000,0 
             [5] = {Point64} 100000,500500,0 
             [6] = {Point64} 200000,500500,0 
             [7] = {Point64} 200000,499500,0 
            [1] = {List<Point64>} Count = 8
             [0] = {Point64} 99500,-100500,0 
             [1] = {Point64} 99500,-100000,0 
             [2] = {Point64} -100000,-100000,0 
             [3] = {Point64} -100000,100000,0 
             [4] = {Point64} 100000,100000,0 
             [5] = {Point64} 100000,-99500,0 
             [6] = {Point64} 200000,-99500,0 
             [7] = {Point64} 200000,-100500,0 
           */
        
        // Sliver removal test
        Paths sR = GeoWrangler.gapRemoval(sL, -100);
        
        /* Expected output
           sR = {List<List<Point64>>} Count = 2
            [0] = {List<Point64>} Count = 9
             [0] = {Point64} -100000,500000,0 
             [1] = {Point64} -100000,700000,0 
             [2] = {Point64} 100000,700000,0 
             [3] = {Point64} 100000,500500,0 
             [4] = {Point64} 200000,500500,0 
             [5] = {Point64} 200000,499500,0 
             [6] = {Point64} 99500,499500,0 
             [7] = {Point64} 99500,500000,0 
             [8] = {Point64} -100000,500000,0 
            [1] = {List<Point64>} Count = 9
             [0] = {Point64} -100000,-100000,0 
             [1] = {Point64} -100000,100000,0 
             [2] = {Point64} 100000,100000,0 
             [3] = {Point64} 100000,-99500,0 
             [4] = {Point64} 200000,-99500,0 
             [5] = {Point64} 200000,-100500,0 
             [6] = {Point64} 99500,-100500,0 
             [7] = {Point64} 99500,-100000,0 
             [8] = {Point64} -100000,-100000,0 
           */
    }

    private static void complex_islandTest()
    {
        // Island 1 - mix of sliver and gap at the end.
        // Manually create the sliver.
        Path outer = new()
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

        Path inner1 = new()
        {
            new(100000, -100000),
            new(100000, 100000),
            new(200000, 100000),
            new(200000, -100000),
            new(100000, -100000)
        };

        Paths kHSource = new()
        {
            outer,
            inner1
        };

        // Island 2 - simple single hole
        Path outer2 = new()
        {
            new(-200000, 400000),
            new(200000, 400000),
            new(200000, 800000),
            new(-200000, 800000),
            new(-200000, 400000)
        };

        Path inner2 = new()
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
        Path outer3 = new()
        {
            new(-300000, 1200000),
            new(300000, 1200000),
            new(300000, 1600000),
            new(-300000, 1600000),
            new(-300000, 1200000)
        };

        Path inner3_1 = new()
        {
            new(-200000, 1300000),
            new(-200000, 1500000),
            new(-100000, 1500000),
            new(-100000, 1300000),
            new(-200000, 1300000)
        };

        Path inner3_2 = new()
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

        /* Expected
           kHSource = {List<List<Point64>>} Count = 7
            [0] = {List<Point64>} Count = 9
             [0] = {Point64} -200000,-200000,0 
             [1] = {Point64} 300000,-200000,0 
             [2] = {Point64} 300000,200000,0 
             [3] = {Point64} -200000,200000,0 
             [4] = {Point64} -200000,-99900,0 
             [5] = {Point64} -300000,-99900,0 
             [6] = {Point64} -300000,-100100,0 
             [7] = {Point64} -200000,-100100,0 
             [8] = {Point64} -200000,-200000,0 
            [1] = {List<Point64>} Count = 5
             [0] = {Point64} 100000,-100000,0 
             [1] = {Point64} 100000,100000,0 
             [2] = {Point64} 200000,100000,0 
             [3] = {Point64} 200000,-100000,0 
             [4] = {Point64} 100000,-100000,0 
            [2] = {List<Point64>} Count = 5
             [0] = {Point64} -200000,400000,0 
             [1] = {Point64} 200000,400000,0 
             [2] = {Point64} 200000,800000,0 
             [3] = {Point64} -200000,800000,0 
             [4] = {Point64} -200000,400000,0 
            [3] = {List<Point64>} Count = 5
             [0] = {Point64} -100000,500000,0 
             [1] = {Point64} -100000,700000,0 
             [2] = {Point64} 100000,700000,0 
             [3] = {Point64} 100000,500000,0 
             [4] = {Point64} -100000,500000,0 
            [4] = {List<Point64>} Count = 5
             [0] = {Point64} -300000,1200000,0 
             [1] = {Point64} 300000,1200000,0 
             [2] = {Point64} 300000,1600000,0 
             [3] = {Point64} -300000,1600000,0 
             [4] = {Point64} -300000,1200000,0 
            [5] = {List<Point64>} Count = 5
             [0] = {Point64} -200000,1300000,0 
             [1] = {Point64} -200000,1500000,0 
             [2] = {Point64} -100000,1500000,0 
             [3] = {Point64} -100000,1300000,0 
             [4] = {Point64} -200000,1300000,0 
            [6] = {List<Point64>} Count = 5
             [0] = {Point64} 100000,1300000,0 
             [1] = {Point64} 100000,1500000,0 
             [2] = {Point64} 200000,1500000,0 
             [3] = {Point64} 200000,1300000,0 
             [4] = {Point64} 100000,1300000,0 
           */
        
        Paths kH = GeoWrangler.makeKeyHole(kHSource, reverseEval:false, biDirectionalEval:true);

        /* Expected output
           kH = {List<List<Point64>>} Count = 3
            [0] = {List<Point64>} Count = 21
             [0] = {Point64} 300000,1200000,0 
             [1] = {Point64} 300000,1299500,0 
             [2] = {Point64} 199500,1299500,0 
             [3] = {Point64} 199500,1300000,0 
             [4] = {Point64} 100000,1300000,0 
             [5] = {Point64} 100000,1500000,0 
             [6] = {Point64} 200000,1500000,0 
             [7] = {Point64} 200000,1300500,0 
             [8] = {Point64} 300000,1300500,0 
             [9] = {Point64} 300000,1600000,0 
             [10] = {Point64} -300000,1600000,0 
             [11] = {Point64} -300000,1500500,0 
             [12] = {Point64} -199500,1500500,0 
             [13] = {Point64} -199500,1500000,0 
             [14] = {Point64} -100000,1500000,0 
             [15] = {Point64} -100000,1300000,0 
             [16] = {Point64} -200000,1300000,0 
             [17] = {Point64} -200000,1499500,0 
             [18] = {Point64} -300000,1499500,0 
             [19] = {Point64} -300000,1200000,0 
             [20] = {Point64} 300000,1200000,0 
            [1] = {List<Point64>} Count = 13
             [0] = {Point64} 200000,400000,0 
             [1] = {Point64} 200000,499500,0 
             [2] = {Point64} 99500,499500,0 
             [3] = {Point64} 99500,500000,0 
             [4] = {Point64} -100000,500000,0 
             [5] = {Point64} -100000,700000,0 
             [6] = {Point64} 100000,700000,0 
             [7] = {Point64} 100000,500500,0 
             [8] = {Point64} 200000,500500,0 
             [9] = {Point64} 200000,800000,0 
             [10] = {Point64} -200000,800000,0 
             [11] = {Point64} -200000,400000,0 
             [12] = {Point64} 200000,400000,0 
            [2] = {List<Point64>} Count = 17
             [0] = {Point64} -300000,-100100,0 
             [1] = {Point64} -200000,-100100,0 
             [2] = {Point64} -200000,-200000,0 
             [3] = {Point64} 300000,-200000,0 
             [4] = {Point64} 300000,-100500,0 
             [5] = {Point64} 199500,-100500,0 
             [6] = {Point64} 199500,-100000,0 
             [7] = {Point64} 100000,-100000,0 
             [8] = {Point64} 100000,100000,0 
             [9] = {Point64} 200000,100000,0 
             [10] = {Point64} 200000,-99500,0 
             [11] = {Point64} 300000,-99500,0 
             [12] = {Point64} 300000,200000,0 
             [13] = {Point64} -200000,200000,0 
             [14] = {Point64} -200000,-99900,0 
             [15] = {Point64} -300000,-99900,0 
             [16] = {Point64} -300000,-100100,0 
           */
        
        // Gap removal test
        Paths gR = GeoWrangler.gapRemoval(kH, 100);

        /* Expected output
           gR = {List<List<Point64>>} Count = 3
            [0] = {List<Point64>} Count = 21
             [0] = {Point64} 300000,1200000,0 
             [1] = {Point64} 300000,1299500,0 
             [2] = {Point64} 199500,1299500,0 
             [3] = {Point64} 199500,1300000,0 
             [4] = {Point64} 100000,1300000,0 
             [5] = {Point64} 100000,1500000,0 
             [6] = {Point64} 200000,1500000,0 
             [7] = {Point64} 200000,1300500,0 
             [8] = {Point64} 300000,1300500,0 
             [9] = {Point64} 300000,1600000,0 
             [10] = {Point64} -300000,1600000,0 
             [11] = {Point64} -300000,1500500,0 
             [12] = {Point64} -199500,1500500,0 
             [13] = {Point64} -199500,1500000,0 
             [14] = {Point64} -100000,1500000,0 
             [15] = {Point64} -100000,1300000,0 
             [16] = {Point64} -200000,1300000,0 
             [17] = {Point64} -200000,1499500,0 
             [18] = {Point64} -300000,1499500,0 
             [19] = {Point64} -300000,1200000,0 
             [20] = {Point64} 300000,1200000,0 
            [1] = {List<Point64>} Count = 13
             [0] = {Point64} 200000,400000,0 
             [1] = {Point64} 200000,499500,0 
             [2] = {Point64} 99500,499500,0 
             [3] = {Point64} 99500,500000,0 
             [4] = {Point64} -100000,500000,0 
             [5] = {Point64} -100000,700000,0 
             [6] = {Point64} 100000,700000,0 
             [7] = {Point64} 100000,500500,0 
             [8] = {Point64} 200000,500500,0 
             [9] = {Point64} 200000,800000,0 
             [10] = {Point64} -200000,800000,0 
             [11] = {Point64} -200000,400000,0 
             [12] = {Point64} 200000,400000,0 
            [2] = {List<Point64>} Count = 17
             [0] = {Point64} -300000,-100100,0 
             [1] = {Point64} -200000,-100100,0 
             [2] = {Point64} -200000,-200000,0 
             [3] = {Point64} 300000,-200000,0 
             [4] = {Point64} 300000,-100500,0 
             [5] = {Point64} 199500,-100500,0 
             [6] = {Point64} 199500,-100000,0 
             [7] = {Point64} 100000,-100000,0 
             [8] = {Point64} 100000,100000,0 
             [9] = {Point64} 200000,100000,0 
             [10] = {Point64} 200000,-99500,0 
             [11] = {Point64} 300000,-99500,0 
             [12] = {Point64} 300000,200000,0 
             [13] = {Point64} -200000,200000,0 
             [14] = {Point64} -200000,-99900,0 
             [15] = {Point64} -300000,-99900,0 
             [16] = {Point64} -300000,-100100,0 
           */
        
        // Sliver removal test
        Paths sR = GeoWrangler.gapRemoval(kH, -100);
        
        /* Expected output
           sR = {List<List<Point64>>} Count = 3
            [0] = {List<Point64>} Count = 21
             [0] = {Point64} 300000,1200000,0 
             [1] = {Point64} 300000,1299500,0 
             [2] = {Point64} 199500,1299500,0 
             [3] = {Point64} 199500,1300000,0 
             [4] = {Point64} 100000,1300000,0 
             [5] = {Point64} 100000,1500000,0 
             [6] = {Point64} 200000,1500000,0 
             [7] = {Point64} 200000,1300500,0 
             [8] = {Point64} 300000,1300500,0 
             [9] = {Point64} 300000,1600000,0 
             [10] = {Point64} -300000,1600000,0 
             [11] = {Point64} -300000,1500500,0 
             [12] = {Point64} -199500,1500500,0 
             [13] = {Point64} -199500,1500000,0 
             [14] = {Point64} -100000,1500000,0 
             [15] = {Point64} -100000,1300000,0 
             [16] = {Point64} -200000,1300000,0 
             [17] = {Point64} -200000,1499500,0 
             [18] = {Point64} -300000,1499500,0 
             [19] = {Point64} -300000,1200000,0 
             [20] = {Point64} 300000,1200000,0 
            [1] = {List<Point64>} Count = 13
             [0] = {Point64} 200000,400000,0 
             [1] = {Point64} 200000,499500,0 
             [2] = {Point64} 99500,499500,0 
             [3] = {Point64} 99500,500000,0 
             [4] = {Point64} -100000,500000,0 
             [5] = {Point64} -100000,700000,0 
             [6] = {Point64} 100000,700000,0 
             [7] = {Point64} 100000,500500,0 
             [8] = {Point64} 200000,500500,0 
             [9] = {Point64} 200000,800000,0 
             [10] = {Point64} -200000,800000,0 
             [11] = {Point64} -200000,400000,0 
             [12] = {Point64} 200000,400000,0 
            [2] = {List<Point64>} Count = 13
             [0] = {Point64} 300000,-200000,0 
             [1] = {Point64} 300000,-100500,0 
             [2] = {Point64} 199500,-100500,0 
             [3] = {Point64} 199500,-100000,0 
             [4] = {Point64} 100000,-100000,0 
             [5] = {Point64} 100000,100000,0 
             [6] = {Point64} 200000,100000,0 
             [7] = {Point64} 200000,-99500,0 
             [8] = {Point64} 300000,-99500,0 
             [9] = {Point64} 300000,200000,0 
             [10] = {Point64} -200000,200000,0 
             [11] = {Point64} -200000,-200000,0 
             [12] = {Point64} 300000,-200000,0 
           */
        
    }

    private static void multiHoleTest()
    {
        Paths paths = new();
        Path path_1 = new()
        {
            new(50000, 0),
            new(0, 0),
            new(0, 155000),
            new(50000, 155000),
            new(50000, 0)
        };
        paths.Add(path_1);

        Path path_2 = new()
        {
            new(35000, 135000),
            new(35000, 150000),
            new(5000, 150000),
            new(5000, 135000),
            new(35000, 135000)
        };
        paths.Add(path_2);

        Path path_3 = new()
        {
            new(22000, 95000),
            new(22000, 125000),
            new(5000, 125000),
            new(5000, 95000),
            new(22000, 95000)
        };
        paths.Add(path_3);

        Path path_4 = new()
        {
            new(35000, 45000),
            new(35000, 75000),
            new(5000, 75000),
            new(5000, 45000),
            new(35000, 45000)
        };
        paths.Add(path_4);

        Path path_5 = new()
        {
            new(35000, 5000),
            new(35000, 35000),
            new(5000, 35000),
            new(5000, 5000),
            new(35000, 5000)
        };
        paths.Add(path_5);

        /* Expected
           paths = {List<List<Point64>>} Count = 5
            [0] = {List<Point64>} Count = 5
             [0] = {Point64} 50000,0,0 
             [1] = {Point64} 0,0,0 
             [2] = {Point64} 0,155000,0 
             [3] = {Point64} 50000,155000,0 
             [4] = {Point64} 50000,0,0 
            [1] = {List<Point64>} Count = 5
             [0] = {Point64} 35000,135000,0 
             [1] = {Point64} 35000,150000,0 
             [2] = {Point64} 5000,150000,0 
             [3] = {Point64} 5000,135000,0 
             [4] = {Point64} 35000,135000,0 
            [2] = {List<Point64>} Count = 5
             [0] = {Point64} 22000,95000,0 
             [1] = {Point64} 22000,125000,0 
             [2] = {Point64} 5000,125000,0 
             [3] = {Point64} 5000,95000,0 
             [4] = {Point64} 22000,95000,0 
            [3] = {List<Point64>} Count = 5
             [0] = {Point64} 35000,45000,0 
             [1] = {Point64} 35000,75000,0 
             [2] = {Point64} 5000,75000,0 
             [3] = {Point64} 5000,45000,0 
             [4] = {Point64} 35000,45000,0 
            [4] = {List<Point64>} Count = 5
             [0] = {Point64} 35000,5000,0 
             [1] = {Point64} 35000,35000,0 
             [2] = {Point64} 5000,35000,0 
             [3] = {Point64} 5000,5000,0 
             [4] = {Point64} 35000,5000,0 
           */
        
        // Generate keyholed geometry
        Paths kH = GeoWrangler.makeKeyHole(paths, reverseEval:false, biDirectionalEval:true);
        
        /* Expected output
           kH = {List<List<Point64>>} Count = 1
            [0] = {List<Point64>} Count = 37
             [0] = {Point64} 50000,0,0 
             [1] = {Point64} 50000,155000,0 
             [2] = {Point64} 0,155000,0 
             [3] = {Point64} 0,150500,0 
             [4] = {Point64} 5500,150500,0 
             [5] = {Point64} 5500,150000,0 
             [6] = {Point64} 35000,150000,0 
             [7] = {Point64} 35000,135000,0 
             [8] = {Point64} 5000,135000,0 
             [9] = {Point64} 5000,149500,0 
             [10] = {Point64} 0,149500,0 
             [11] = {Point64} 0,125500,0 
             [12] = {Point64} 5500,125500,0 
             [13] = {Point64} 5500,125000,0 
             [14] = {Point64} 22000,125000,0 
             [15] = {Point64} 22000,95000,0 
             [16] = {Point64} 5000,95000,0 
             [17] = {Point64} 5000,124500,0 
             [18] = {Point64} 0,124500,0 
             [19] = {Point64} 0,75500,0 
             [20] = {Point64} 5500,75500,0 
             [21] = {Point64} 5500,75000,0 
             [22] = {Point64} 35000,75000,0 
             [23] = {Point64} 35000,45000,0 
             [24] = {Point64} 5000,45000,0 
             [25] = {Point64} 5000,74500,0 
             [26] = {Point64} 0,74500,0 
             [27] = {Point64} 0,35500,0 
             [28] = {Point64} 5500,35500,0 
             [29] = {Point64} 5500,35000,0 
             [30] = {Point64} 35000,35000,0 
             [31] = {Point64} 35000,5000,0 
             [32] = {Point64} 5000,5000,0 
             [33] = {Point64} 5000,34500,0 
             [34] = {Point64} 0,34500,0 
             [35] = {Point64} 0,0,0 
             [36] = {Point64} 50000,0,0 
           */
    }
}