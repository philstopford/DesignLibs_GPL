using Clipper2Lib;
using geoWrangler;
using geoCoreLib;
using PartitionTestGeometrySource;

namespace partitionTest;

internal class Program
{
    private static void Main(string[] args)
    {
        Console.WriteLine("Part One");
        partOne();

        Console.WriteLine("Part Two");
        partTwo();

        Console.WriteLine("Part Three (takes a while)");
        partThree();

        Console.WriteLine("Part Four (takes even longer)");
        partFour();
        
        Console.WriteLine("Part Five");
        partFive();
    }
    
    private static void partOne()
    {
        // L
        PathD L = TestGeometry.getL();
        PathD rL = TestGeometry.getRL();

        // Reversed orientation.
        PathD L_ccw = new (L);
        L_ccw.Reverse();

        // U
        PathD U = TestGeometry.getU();

        // T
        PathD T = TestGeometry.getT();

        // X
        PathD X = TestGeometry.getX();

        // S
        PathD S = TestGeometry.getS();

        // Negative S
        PathD nS = TestGeometry.getnegS();

        // Complex 1
        PathD C1 = TestGeometry.getComplex1();


        // Complex 2
        PathD C2 = TestGeometry.getComplex2();

        // Complex 3
        PathD C3 = TestGeometry.getComplex3();

        // C3 = GeoWrangler.clockwiseAndReorder(C3);

        bool orth = GeoWrangler.orthogonal(C2, angularTolerance: 0);
        bool orth2 = GeoWrangler.orthogonal(C3, angularTolerance: 0);

        // Complex 10, rot 15
        PathD C10R15 = TestGeometry.getComplex10rot15();


        // Staircase
        PathD S1 = TestGeometry.getStaircase();

        // Rectangular decomposition will return non-orthogonal polygons in the output when encountered.

        int rayLength = 1000;

        bool abort = false;

        PathsD l = GeoWrangler.rectangular_decomposition(ref abort, L, maxRayLength: rayLength);

        /* Expected output
            l = {List<Path>} Count = 2
             [0] = {Path} GeoLibPoint[4]
              [0] = GeoLibPoint
               X = {int} 0
               Y = {int} 0
               tag = {int} 0
              [1] = GeoLibPoint
               X = {int} 0
               Y = {int} 50
               tag = {int} 0
              [2] = GeoLibPoint
               X = {int} 10
               Y = {int} 50
               tag = {int} 0
              [3] = GeoLibPoint
               X = {int} 10
               Y = {int} 0
               tag = {int} 0
             [1] = {Path} GeoLibPoint[4]
              [0] = GeoLibPoint
               X = {int} 10
               Y = {int} 0
               tag = {int} 0
              [1] = GeoLibPoint
               X = {int} 10
               Y = {int} 20
               tag = {int} 0
              [2] = GeoLibPoint
               X = {int} 60
               Y = {int} 20
               tag = {int} 0
              [3] = GeoLibPoint
               X = {int} 60
               Y = {int} 0
               tag = {int} 0         
         */

        writeToLayout("l", L, l);

        PathD L_small = Clipper.ScalePath(L, 0.01);
        PathsD lsmall = GeoWrangler.rectangular_decomposition(ref abort, L_small, maxRayLength: rayLength);
        writeToLayout("l_small", L_small, lsmall);

        PathsD lccw = GeoWrangler.rectangular_decomposition(ref abort, L_ccw, maxRayLength: rayLength);

        /* Expected output
           lccw = {List<Path>} Count = 2
            [0] = {Path} GeoLibPoint[4]
             [0] = GeoLibPoint
              X = {int} 0
              Y = {int} 0
              tag = {int} 0
             [1] = GeoLibPoint
              X = {int} 0
              Y = {int} 50
              tag = {int} 0
             [2] = GeoLibPoint
              X = {int} 10
              Y = {int} 50
              tag = {int} 0
             [3] = GeoLibPoint
              X = {int} 10
              Y = {int} 0
              tag = {int} 0
            [1] = {Path} GeoLibPoint[4]
             [0] = GeoLibPoint
              X = {int} 10
              Y = {int} 0
              tag = {int} 0
             [1] = GeoLibPoint
              X = {int} 10
              Y = {int} 20
              tag = {int} 0
             [2] = GeoLibPoint
              X = {int} 60
              Y = {int} 20
              tag = {int} 0
             [3] = GeoLibPoint
              X = {int} 60
              Y = {int} 0
              tag = {int} 0
            */

        writeToLayout("lccw", L_ccw, lccw);

        PathsD rl = GeoWrangler.rectangular_decomposition(ref abort, rL, maxRayLength: rayLength);

        /* Expected output
           rl = {List<Path>} Count = 2
            [0] = {Path} GeoLibPoint[4]
             [0] = GeoLibPoint
              X = {int} 10
              Y = {int} 0
              tag = {int} 0
             [1] = GeoLibPoint
              X = {int} 10
              Y = {int} 50
              tag = {int} 0
             [2] = GeoLibPoint
              X = {int} 60
              Y = {int} 50
              tag = {int} 0
             [3] = GeoLibPoint
              X = {int} 60
              Y = {int} 0
              tag = {int} 0
            [1] = {Path} GeoLibPoint[4]
             [0] = GeoLibPoint
              X = {int} 0
              Y = {int} 0
              tag = {int} 0
             [1] = GeoLibPoint
              X = {int} 0
              Y = {int} 20
              tag = {int} 0
             [2] = GeoLibPoint
              X = {int} 10
              Y = {int} 20
              tag = {int} 0
             [3] = GeoLibPoint
              X = {int} 10
              Y = {int} 0
              tag = {int} 0
            */

        writeToLayout("rl", rL, rl);

        PathsD u = GeoWrangler.rectangular_decomposition(ref abort, U, maxRayLength: rayLength);

        /* Expected output
           u = {List<Path>} Count = 3
            [0] = {Path} GeoLibPoint[4]
             [0] = GeoLibPoint
              X = {int} 0
              Y = {int} 0
              tag = {int} 0
             [1] = GeoLibPoint
              X = {int} 0
              Y = {int} 50
              tag = {int} 0
             [2] = GeoLibPoint
              X = {int} 10
              Y = {int} 50
              tag = {int} 0
             [3] = GeoLibPoint
              X = {int} 10
              Y = {int} 0
              tag = {int} 0
            [1] = {Path} GeoLibPoint[4]
             [0] = GeoLibPoint
              X = {int} 60
              Y = {int} 0
              tag = {int} 0
             [1] = GeoLibPoint
              X = {int} 60
              Y = {int} 80
              tag = {int} 0
             [2] = GeoLibPoint
              X = {int} 120
              Y = {int} 80
              tag = {int} 0
             [3] = GeoLibPoint
              X = {int} 120
              Y = {int} 0
              tag = {int} 0
            [2] = {Path} GeoLibPoint[4]
             [0] = GeoLibPoint
              X = {int} 10
              Y = {int} 0
              tag = {int} 0
             [1] = GeoLibPoint
              X = {int} 10
              Y = {int} 20
              tag = {int} 0
             [2] = GeoLibPoint
              X = {int} 60
              Y = {int} 20
              tag = {int} 0
             [3] = GeoLibPoint
              X = {int} 60
              Y = {int} 0
              tag = {int} 0
          */

        writeToLayout("u", U, u);

        PathsD t = GeoWrangler.rectangular_decomposition(ref abort, T, maxRayLength: rayLength);

        /* Expected output
           t = {List<Path>} Count = 3
            [0] = {Path} GeoLibPoint[4]
             [0] = GeoLibPoint
              X = {int} 60
              Y = {int} 50
              tag = {int} 0
             [1] = GeoLibPoint
              X = {int} 60
              Y = {int} 80
              tag = {int} 0
             [2] = GeoLibPoint
              X = {int} 80
              Y = {int} 80
              tag = {int} 0
             [3] = GeoLibPoint
              X = {int} 80
              Y = {int} 50
              tag = {int} 0
            [1] = {Path} GeoLibPoint[4]
             [0] = GeoLibPoint
              X = {int} 40
              Y = {int} 0
              tag = {int} 0
             [1] = GeoLibPoint
              X = {int} 40
              Y = {int} 80
              tag = {int} 0
             [2] = GeoLibPoint
              X = {int} 60
              Y = {int} 80
              tag = {int} 0
             [3] = GeoLibPoint
              X = {int} 60
              Y = {int} 0
              tag = {int} 0
            [2] = {Path} GeoLibPoint[4]
             [0] = GeoLibPoint
              X = {int} 0
              Y = {int} 50
              tag = {int} 0
             [1] = GeoLibPoint
              X = {int} 0
              Y = {int} 80
              tag = {int} 0
             [2] = GeoLibPoint
              X = {int} 40
              Y = {int} 80
              tag = {int} 0
             [3] = GeoLibPoint
              X = {int} 40
              Y = {int} 50
              tag = {int} 0
         */

        writeToLayout("t", T, t);

        PathsD x = GeoWrangler.rectangular_decomposition(ref abort, X, maxRayLength: rayLength);

        /* Expected output
           x = {List<Path>} Count = 3
            [0] = {Path} GeoLibPoint[4]
             [0] = GeoLibPoint
              X = {int} 0
              Y = {int} 50
              tag = {int} 0
             [1] = GeoLibPoint
              X = {int} 0
              Y = {int} 80
              tag = {int} 0
             [2] = GeoLibPoint
              X = {int} 60
              Y = {int} 80
              tag = {int} 0
             [3] = GeoLibPoint
              X = {int} 60
              Y = {int} 50
              tag = {int} 0
            [1] = {Path} GeoLibPoint[4]
             [0] = GeoLibPoint
              X = {int} 60
              Y = {int} 20
              tag = {int} 0
             [1] = GeoLibPoint
              X = {int} 60
              Y = {int} 100
              tag = {int} 0
             [2] = GeoLibPoint
              X = {int} 80
              Y = {int} 100
              tag = {int} 0
             [3] = GeoLibPoint
              X = {int} 80
              Y = {int} 20
              tag = {int} 0
            [2] = {Path} GeoLibPoint[4]
             [0] = GeoLibPoint
              X = {int} 80
              Y = {int} 50
              tag = {int} 0
             [1] = GeoLibPoint
              X = {int} 80
              Y = {int} 80
              tag = {int} 0
             [2] = GeoLibPoint
              X = {int} 100
              Y = {int} 80
              tag = {int} 0
             [3] = GeoLibPoint
              X = {int} 100
              Y = {int} 50
              tag = {int} 0
         */

        writeToLayout("x", X, x);

        PathsD s = GeoWrangler.rectangular_decomposition(ref abort, S, maxRayLength: rayLength);

        /* Expected output
           s = {List<Path>} Count = 5
            [0] = {Path} GeoLibPoint[4]
             [0] = GeoLibPoint
              X = {int} 0
              Y = {int} 50
              tag = {int} 0
             [1] = GeoLibPoint
              X = {int} 0
              Y = {int} 110
              tag = {int} 0
             [2] = GeoLibPoint
              X = {int} 20
              Y = {int} 110
              tag = {int} 0
             [3] = GeoLibPoint
              X = {int} 20
              Y = {int} 50
              tag = {int} 0
            [1] = {Path} GeoLibPoint[4]
             [0] = GeoLibPoint
              X = {int} 0
              Y = {int} 0
              tag = {int} 0
             [1] = GeoLibPoint
              X = {int} 0
              Y = {int} 20
              tag = {int} 0
             [2] = GeoLibPoint
              X = {int} 20
              Y = {int} 20
              tag = {int} 0
             [3] = GeoLibPoint
              X = {int} 20
              Y = {int} 0
              tag = {int} 0
            [2] = {Path} GeoLibPoint[4]
             [0] = GeoLibPoint
              X = {int} 80
              Y = {int} 0
              tag = {int} 0
             [1] = GeoLibPoint
              X = {int} 80
              Y = {int} 60
              tag = {int} 0
             [2] = GeoLibPoint
              X = {int} 100
              Y = {int} 60
              tag = {int} 0
             [3] = GeoLibPoint
              X = {int} 100
              Y = {int} 0
              tag = {int} 0
            [3] = {Path} GeoLibPoint[4]
             [0] = GeoLibPoint
              X = {int} 80
              Y = {int} 80
              tag = {int} 0
             [1] = GeoLibPoint
              X = {int} 80
              Y = {int} 110
              tag = {int} 0
             [2] = GeoLibPoint
              X = {int} 100
              Y = {int} 110
              tag = {int} 0
             [3] = GeoLibPoint
              X = {int} 100
              Y = {int} 80
              tag = {int} 0
            [4] = {Path} GeoLibPoint[4]
             [0] = GeoLibPoint
              X = {int} 20
              Y = {int} 0
              tag = {int} 0
             [1] = GeoLibPoint
              X = {int} 20
              Y = {int} 110
              tag = {int} 0
             [2] = GeoLibPoint
              X = {int} 80
              Y = {int} 110
              tag = {int} 0
             [3] = GeoLibPoint
              X = {int} 80
              Y = {int} 0
              tag = {int} 0
         */

        writeToLayout("s", S, s);

        PathsD ns = GeoWrangler.rectangular_decomposition(ref abort, nS, maxRayLength: rayLength);

        /* Expected output
           ns = {List<Path>} Count = 5
            [0] = {Path} GeoLibPoint[4]
             [0] = GeoLibPoint
              X = {int} 0
              Y = {int} -150
              tag = {int} 0
             [1] = GeoLibPoint
              X = {int} 0
              Y = {int} -90
              tag = {int} 0
             [2] = GeoLibPoint
              X = {int} 20
              Y = {int} -90
              tag = {int} 0
             [3] = GeoLibPoint
              X = {int} 20
              Y = {int} -150
              tag = {int} 0
            [1] = {Path} GeoLibPoint[4]
             [0] = GeoLibPoint
              X = {int} 0
              Y = {int} -200
              tag = {int} 0
             [1] = GeoLibPoint
              X = {int} 0
              Y = {int} -180
              tag = {int} 0
             [2] = GeoLibPoint
              X = {int} 20
              Y = {int} -180
              tag = {int} 0
             [3] = GeoLibPoint
              X = {int} 20
              Y = {int} -200
              tag = {int} 0
            [2] = {Path} GeoLibPoint[4]
             [0] = GeoLibPoint
              X = {int} 80
              Y = {int} -200
              tag = {int} 0
             [1] = GeoLibPoint
              X = {int} 80
              Y = {int} -140
              tag = {int} 0
             [2] = GeoLibPoint
              X = {int} 100
              Y = {int} -140
              tag = {int} 0
             [3] = GeoLibPoint
              X = {int} 100
              Y = {int} -200
              tag = {int} 0
            [3] = {Path} GeoLibPoint[4]
             [0] = GeoLibPoint
              X = {int} 80
              Y = {int} -120
              tag = {int} 0
             [1] = GeoLibPoint
              X = {int} 80
              Y = {int} -90
              tag = {int} 0
             [2] = GeoLibPoint
              X = {int} 100
              Y = {int} -90
              tag = {int} 0
             [3] = GeoLibPoint
              X = {int} 100
              Y = {int} -120
              tag = {int} 0
            [4] = {Path} GeoLibPoint[4]
             [0] = GeoLibPoint
              X = {int} 20
              Y = {int} -200
              tag = {int} 0
             [1] = GeoLibPoint
              X = {int} 20
              Y = {int} -90
              tag = {int} 0
             [2] = GeoLibPoint
              X = {int} 80
              Y = {int} -90
              tag = {int} 0
             [3] = GeoLibPoint
              X = {int} 80
              Y = {int} -200
              tag = {int} 0
         */

        writeToLayout("ns", nS, ns);

        PathsD c1 = GeoWrangler.rectangular_decomposition(ref abort, C1, maxRayLength: rayLength);

        // Expect 17 quads

        writeToLayout("c1", C1, c1);

        PathsD c2 = GeoWrangler.rectangular_decomposition(ref abort, C2, maxRayLength: rayLength);

        // Expect 81 quads

        writeToLayout("c2", C2, c2);

        PathsD c3 = GeoWrangler.rectangular_decomposition(ref abort, C3, maxRayLength: rayLength);

        // Expect 13 quads

        writeToLayout("c3", C3, c3);

        PathsD c10r15 = GeoWrangler.rectangular_decomposition(ref abort, C10R15, maxRayLength: rayLength);

        // Expect 1 irregular polygon

        writeToLayout("c10r15", C10R15, c10r15);

        PathsD s1 = GeoWrangler.rectangular_decomposition(ref abort, S1, maxRayLength: rayLength);

        /* Expected output
           s1 = {List<Path>} Count = 4
            [0] = {Path} GeoLibPoint[4]
             [0] = GeoLibPoint
              X = {int} -50
              Y = {int} -50
              tag = {int} 0
             [1] = GeoLibPoint
              X = {int} -50
              Y = {int} 0
              tag = {int} 0
             [2] = GeoLibPoint
              X = {int} 0
              Y = {int} 0
              tag = {int} 0
             [3] = GeoLibPoint
              X = {int} 0
              Y = {int} -50
              tag = {int} 0
            [1] = {Path} GeoLibPoint[4]
             [0] = GeoLibPoint
              X = {int} 0
              Y = {int} -50
              tag = {int} 0
             [1] = GeoLibPoint
              X = {int} 0
              Y = {int} 120
              tag = {int} 0
             [2] = GeoLibPoint
              X = {int} 100
              Y = {int} 120
              tag = {int} 0
             [3] = GeoLibPoint
              X = {int} 100
              Y = {int} -50
              tag = {int} 0
            [2] = {Path} GeoLibPoint[4]
             [0] = GeoLibPoint
              X = {int} 150
              Y = {int} -50
              tag = {int} 0
             [1] = GeoLibPoint
              X = {int} 150
              Y = {int} 300
              tag = {int} 0
             [2] = GeoLibPoint
              X = {int} 200
              Y = {int} 300
              tag = {int} 0
             [3] = GeoLibPoint
              X = {int} 200
              Y = {int} -50
              tag = {int} 0
            [3] = {Path} GeoLibPoint[4]
             [0] = GeoLibPoint
              X = {int} 100
              Y = {int} -50
              tag = {int} 0
             [1] = GeoLibPoint
              X = {int} 100
              Y = {int} 200
              tag = {int} 0
             [2] = GeoLibPoint
              X = {int} 150
              Y = {int} 200
              tag = {int} 0
             [3] = GeoLibPoint
              X = {int} 150
              Y = {int} -50
              tag = {int} 0
         */

        writeToLayout("s1", S1, s1);
    }

    private static void partTwo()
    {
        PathsD incoming = new();
        PathD lPieces = new ()
        {
            new(0.00000, 0.00000),
            new(0.00000, 0.05000),
            new(0.01000, 0.05000),
            new(0.01000, 0.00000),
            new(0.00000, 0.00000)
        };


        PathD lPiece2 = new ()
        {
            new(0.01000, 0.00000),
            new(0.01000, 0.02000),
            new(0.06000, 0.02000),
            new(0.06000, 0.00000),
            new(0.01000, 0.00000)
        };

        lPieces.Reverse();
        lPiece2.Reverse();
        incoming.Add(lPieces);
        incoming.Add(lPiece2);
        
        ClipperD c = new();
        c.AddSubject(incoming);
        PathsD ret = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, ret);
    }

    private static void partThree()
    {
        System.Diagnostics.Stopwatch sw = new();
        int rayLength = 1000;

        Console.WriteLine("  Preparing....");
        sw.Start();

        PathD poly = TestGeometry.verycomplex();

        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Geo clean-up....");
        sw.Restart();
        poly = GeoWrangler.clockwiseAndReorderXY(poly);
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Conversion....");
        sw.Restart();
        PathsD done = new() { new(poly) }; // this was originally scaled up by 1000 in the integer pipeline.
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        bool abort = false;

        Console.WriteLine("  Decomposition (vertical)....");
        sw.Restart();
        PathsD ns = GeoWrangler.rectangular_decomposition(ref abort, done, maxRayLength: rayLength);

        // Expect 721 quads.

        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Writing....");
        writeToLayout("complex", done[0], ns);

        Console.WriteLine("  Decomposition (horizontal)....");
        sw.Restart();
        ns = GeoWrangler.rectangular_decomposition(ref abort, done, maxRayLength: rayLength, vertical: false);

        // Expect 721 quads.

        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Writing....");
        writeToLayout("complex_horizontal", done[0], ns);

        Console.WriteLine("  Done.");
    }

    private static void partFour()
    {
        System.Diagnostics.Stopwatch sw = new();
        Console.WriteLine(" Part 1....");
        Console.WriteLine("  Preparing....");
        sw.Start();
        PathD points_1 = TestGeometry.ortho_fractal_1();
        points_1 = GeoWrangler.close(points_1);
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        partFour_do(points_1, "complex_loop");

        Console.WriteLine(" Part 2....");
        Console.WriteLine("  Preparing....");
        sw.Start();
        PathD points_2 = TestGeometry.ortho_fractal_2();
        points_2 = GeoWrangler.close(points_2);
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");
        partFour_do(points_2, "complex_loop_rot");
    }

    private static void partFour_do(PathD points, string baseString)
    {
        System.Diagnostics.Stopwatch sw = new();

        bool vertical = true;
        bool abort = false;

        Console.WriteLine("  Keyhole....");
        // Give the keyholder a whirl:
        sw.Restart();
        PathD toDecomp = GeoWrangler.makeKeyHole(points, reverseEval:false, biDirectionalEval:true)[0];
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Query....");
        sw.Restart();
        PathD bounds = GeoWrangler.getBounds(toDecomp);
        PointD dist = GeoWrangler.distanceBetweenPoints_point(bounds[0], bounds[1]);
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Decomposition (vertical)....");
        sw.Restart();
        PathsD decompOut = GeoWrangler.rectangular_decomposition(ref abort, toDecomp,
            maxRayLength: (long)Math.Max(Math.Abs(dist.x), Math.Abs(dist.y)) * 1, vertical: vertical);
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Writing....");
        writeToLayout(baseString, points, decompOut);

        Console.WriteLine("  Decomposition (horizontal)....");
        sw.Restart();
        decompOut = GeoWrangler.rectangular_decomposition(ref abort, toDecomp,
            maxRayLength: (long)Math.Max(Math.Abs(dist.x), Math.Abs(dist.y)) * 1, vertical: !vertical);
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Writing....");
        writeToLayout(baseString + "_horizontal", points, decompOut);
        Console.WriteLine("  Done.");

        sw.Stop();
    }

    private static void partFive()
    {
        bool vertical = true;
        bool abort = false;

        PathsD polydata = new()
        {
            new()
            {
                new (-40, -30),
                new (-40, 70),
                new (0, 70),
                new (0, 50),
                new (-20, 50),
                new (-20, -10),
                new (0, -10),
                new (0, -30),
                new (-40, -30),
            },
            new()
            {
                new (0, -30),
                new (0, -10),
                new (20, -10),
                new (20, 50),
                new (0, 50),
                new (0, 70),
                new (40, 70),
                new (40, -30),
                new (0, -30),
            },
            new()
            {
                new (-80, -60),
                new (-80, 100),
                new (0, 100),
                new (0, 80),
                new (-60, 80),
                new (-60, -50),
                new (0, -50),
                new (0, -60),
                new (-80, -60),
            },
            new()
            {
                new (0, -60),
                new (0, -50),
                new (60, -50),
                new (60, 80),
                new (0, 80),
                new (0, 100),
                new (80, 100),
                new (80, -60),
                new (0, -60),
            },
            new()
            {
                new (-8, -27),
                new (-8, 40),
                new (9, 40),
                new (9, -27),
                new (-8, -27)
            },
            new()
            {
                new (-14, -4),
                new (-14, 15),
                new (-7, 15),
                new (-7, -4),
                new (-14, -4),
            },
            new()
            {
                new (10, 9),
                new (10, 20),
                new (13, 20),
                new (13, 9),
                new (10, 9),
            },
            new()
            {
                new (48, -1),
                new (48, 31),
                new (55, 31),
                new (55, -1),
                new (48, -1),
            },
            new()
            {
                new (-11, -44),
                new (-11, -39),
                new (16, -39),
                new (16, -44),
                new (-11, -44),
            },
            new()
            {
                new (-51, 3),
                new (-51, 23),
                new (-47, 23),
                new (-47, 3),
                new (-51, 3),
            },
            new()
            {
                new (-16, 76),
                new (-16, 77),
                new (-3, 77),
                new (-3, 76),
                new (-16, 76),
            },
        };
        
        int scaleFactorForOperation = 1000;

        PathsD out_decomp = new();
        for (int i = 0; i < polydata.Count; i++)
        {
         PathD points = new (polydata[i]);
         points = GeoWrangler.removeDuplicates(points);
         points = GeoWrangler.stripColinear(points);
         points = GeoWrangler.clockwiseAndReorderXY(points);
         
         PathD toDecomp = GeoWrangler.makeKeyHole(GeoWrangler.sliverGapRemoval(points), reverseEval:false, biDirectionalEval:false)[0];
         PathD  bounds = GeoWrangler.getBounds(toDecomp);
         PointD dist = GeoWrangler.distanceBetweenPoints_point(bounds[0], bounds[1]);

         PathsD decompOut = GeoWrangler.rectangular_decomposition(ref abort, toDecomp,
          maxRayLength: (long) Math.Max(Math.Abs(dist.x), Math.Abs(dist.y)), vertical: vertical);
         
         out_decomp.AddRange(decompOut);
        }
    }


    private static void writeToLayout(string filename, PathD orig, PathsD decomped)
    {
        // Can the system define geometry and write it correctly to Oasis and GDS files.
        GeoCore g = new();
        g.reset();
        GCDrawingfield drawing_ = new("")
        {
            accyear = 2018,
            accmonth = 12,
            accday = 5,
            acchour = 2,
            accmin = 10,
            accsec = 10,
            modyear = 2018,
            modmonth = 12,
            modday = 5,
            modhour = 2,
            modmin = 10,
            modsec = 10,
            databaseunits = 1000,
            userunits = 1E-3, // 0.001 / 1E-6;
            libname = "noname"
        };

        GCCell gcell = drawing_.addCell();
        gcell.accyear = 2018;
        gcell.accmonth = 12;
        gcell.accday = 5;
        gcell.acchour = 2;
        gcell.accmin = 10;
        gcell.accsec = 10;
        gcell.modyear = 2018;
        gcell.modmonth = 12;
        gcell.modday = 5;
        gcell.modhour = 2;
        gcell.modmin = 10;
        gcell.modsec = 10;

        gcell.cellName = "test";

        gcell.addPolygon(GeoWrangler.path64FromPathD(orig), 1, 0);

        for (int i = 0; i < decomped.Count; i++)
        {
            gcell.addBox(GeoWrangler.path64FromPathD(decomped[i]), 1, 1);
        }

        g.setDrawing(drawing_);
        g.setValid(true);

        gds.gdsWriter gw = new(g, "../../../../../decomp_out/" + filename + "_partitiontest.gds");
        gw.save();

        oasis.oasWriter ow = new(g, "../../../../../decomp_out/" + filename + "_partitiontest.oas");
        ow.save();
    }
}