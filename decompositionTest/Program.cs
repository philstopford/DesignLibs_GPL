using Clipper2Lib;
using geoLib;
using geoWrangler;
using geoCoreLib;
using PartitionTestGeometrySource;

namespace partitionTest;

internal class Program
{
    private static void Main(string[] args)
    {
     debug();
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

    private static void debug()
    {
        Path64 lPoly = Clipper.MakePath(new int[] {
         0,0,
         0,500000,
         100000,500000,
         100000,200000, 
         600000,200000, 
         600000,0,
         0,0
           }
        );
        Path64 newEdge = Clipper.MakePath(new int[]
        {
         100000,200000,
         100000,0
        });

        Paths64 newEdges = new() { newEdge };
        
        
        ClipperOffset co = new() {PreserveCollinear = false};
        co.AddPaths(newEdges, JoinType.Miter, EndType.Square);

        Paths64 cutters = co.Execute(2.0);

        Clipper64 c = new();

        c.AddSubject(lPoly);

        // Take first cutter only - we only cut once, no matter how many potential cutters we have.
        c.AddClip(cutters[0]);
        Paths64 f = new();
        c.Execute(ClipType.Difference, FillRule.EvenOdd, f);

    }
    
    private static void partOne()
    {
        // L
        Path64 L = TestGeometry.getL();
        Path64 rL = TestGeometry.getRL();

        // Reversed orientation.
        Path64 L_ccw = new (L);
        L_ccw.Reverse();

        // U
        Path64 U = TestGeometry.getU();

        // T
        Path64 T = TestGeometry.getT();

        // X
        Path64 X = TestGeometry.getX();

        // S
        Path64 S = TestGeometry.getS();

        // Negative S
        Path64 nS = TestGeometry.getnegS();

        // Complex 1
        Path64 C1 = TestGeometry.getComplex1();


        // Complex 2
        Path64 C2 = TestGeometry.getComplex2();

        // Complex 3
        Path64 C3 = TestGeometry.getComplex3();

        // C3 = GeoWrangler.clockwiseAndReorder(C3);

        bool orth = GeoWrangler.orthogonal(C2, angularTolerance: 0);
        bool orth2 = GeoWrangler.orthogonal(C3, angularTolerance: 0);

        // Complex 10, rot 15
        Path64 C10R15 = TestGeometry.getComplex10rot15();


        // Staircase
        Path64 S1 = TestGeometry.getStaircase();

        // Rectangular decomposition will return non-orthogonal polygons in the output when encountered.

        int rayLength = 1000;

        bool abort = false;

        Paths64 l = GeoWrangler.rectangular_decomposition(ref abort, L, maxRayLength: rayLength);

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

        Paths64 lccw = GeoWrangler.rectangular_decomposition(ref abort, L_ccw, maxRayLength: rayLength);

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

        Paths64 rl = GeoWrangler.rectangular_decomposition(ref abort, rL, maxRayLength: rayLength);

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

        Paths64 u = GeoWrangler.rectangular_decomposition(ref abort, U, maxRayLength: rayLength);

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

        Paths64 t = GeoWrangler.rectangular_decomposition(ref abort, T, maxRayLength: rayLength);

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

        Paths64 x = GeoWrangler.rectangular_decomposition(ref abort, X, maxRayLength: rayLength);

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

        Paths64 s = GeoWrangler.rectangular_decomposition(ref abort, S, maxRayLength: rayLength);

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

        Paths64 ns = GeoWrangler.rectangular_decomposition(ref abort, nS, maxRayLength: rayLength);

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

        Paths64 c1 = GeoWrangler.rectangular_decomposition(ref abort, C1, maxRayLength: rayLength);

        // Expect 17 quads

        writeToLayout("c1", C1, c1);

        Paths64 c2 = GeoWrangler.rectangular_decomposition(ref abort, C2, maxRayLength: rayLength);

        // Expect 81 quads

        writeToLayout("c2", C2, c2);

        Paths64 c3 = GeoWrangler.rectangular_decomposition(ref abort, C3, maxRayLength: rayLength);

        // Expect 13 quads

        writeToLayout("c3", C3, c3);

        Paths64 c10r15 = GeoWrangler.rectangular_decomposition(ref abort, C10R15, maxRayLength: rayLength);

        // Expect 1 irregular polygon

        writeToLayout("c10r15", C10R15, c10r15);

        Paths64 s1 = GeoWrangler.rectangular_decomposition(ref abort, S1, maxRayLength: rayLength);

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

        Paths64 paths = GeoWrangler.paths64FromPathsD(incoming);

        Clipper64 c = new();
        c.AddSubject(paths);
        Paths64 ret = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, ret);

        PathsD done = GeoWrangler.pathsDFromPaths64(ret);
    }

    private static void partThree()
    {
        System.Diagnostics.Stopwatch sw = new();
        int rayLength = 10000;

        Console.WriteLine("  Preparing....");
        sw.Start();

        PathD poly = new ()
        {
            new(0.01000, -0.21300),
            new(0.01000, -0.20200),
            new(0.00000, -0.20200),
            new(0.00000, -0.21200),
            new(-0.01000, -0.21200),
            new(-0.01000, -0.18200),
            new(0.00000, -0.18200),
            new(0.00000, -0.19200),
            new(0.01000, -0.19200),
            new(0.01000, -0.15200),
            new(0.00000, -0.15200),
            new(0.00000, -0.14200),
            new(-0.06000, -0.14200),
            new(-0.06000, -0.13200),
            new(-0.05000, -0.13200),
            new(-0.05000, -0.12200),
            new(-0.06000, -0.12200),
            new(-0.06000, -0.11200),
            new(-0.03000, -0.11200),
            new(-0.03000, -0.12200),
            new(-0.04000, -0.12200),
            new(-0.04000, -0.13200),
            new(0.00000, -0.13200),
            new(0.00000, -0.12200),
            new(0.01000, -0.12200),
            new(0.01000, -0.05000),
            new(0.00000, -0.05000),
            new(0.00000, -0.06000),
            new(-0.01000, -0.06000),
            new(-0.01000, -0.03000),
            new(0.00000, -0.03000),
            new(0.00000, -0.04000),
            new(0.01000, -0.04000),
            new(0.01000, 0.00000),
            new(0.00000, 0.00000),
            new(0.00000, 0.01000),
            new(-0.07000, 0.01000),
            new(-0.07000, 0.00000),
            new(-0.06000, 0.00000),
            new(-0.06000, -0.01000),
            new(-0.09000, -0.01000),
            new(-0.09000, 0.00000),
            new(-0.08000, 0.00000),
            new(-0.08000, 0.01000),
            new(-0.12000, 0.01000),
            new(-0.12000, 0.00000),
            new(-0.13000, 0.00000),
            new(-0.13000, -0.06100),
            new(-0.14000, -0.06100),
            new(-0.14000, -0.05000),
            new(-0.15000, -0.05000),
            new(-0.15000, -0.06000),
            new(-0.16000, -0.06000),
            new(-0.16000, -0.03000),
            new(-0.15000, -0.03000),
            new(-0.15000, -0.04000),
            new(-0.14000, -0.04000),
            new(-0.14000, 0.00000),
            new(-0.15000, 0.00000),
            new(-0.15000, 0.01000),
            new(-0.21000, 0.01000),
            new(-0.21000, 0.02000),
            new(-0.20000, 0.02000),
            new(-0.20000, 0.03000),
            new(-0.21000, 0.03000),
            new(-0.21000, 0.04000),
            new(-0.18000, 0.04000),
            new(-0.18000, 0.03000),
            new(-0.19000, 0.03000),
            new(-0.19000, 0.02000),
            new(-0.15000, 0.02000),
            new(-0.15000, 0.03000),
            new(-0.14000, 0.03000),
            new(-0.14000, 0.09100),
            new(-0.13000, 0.09100),
            new(-0.13000, 0.08000),
            new(-0.12000, 0.08000),
            new(-0.12000, 0.09000),
            new(-0.11000, 0.09000),
            new(-0.11000, 0.06000),
            new(-0.12000, 0.06000),
            new(-0.12000, 0.07000),
            new(-0.13000, 0.07000),
            new(-0.13000, 0.03000),
            new(-0.12000, 0.03000),
            new(-0.12000, 0.02000),
            new(-0.05000, 0.02000),
            new(-0.05000, 0.03000),
            new(-0.06000, 0.03000),
            new(-0.06000, 0.04000),
            new(-0.03000, 0.04000),
            new(-0.03000, 0.03000),
            new(-0.04000, 0.03000),
            new(-0.04000, 0.02000),
            new(0.00000, 0.02000),
            new(0.00000, 0.03000),
            new(0.01000, 0.03000),
            new(0.01000, 0.10200),
            new(0.00000, 0.10200),
            new(0.00000, 0.09200),
            new(-0.01000, 0.09200),
            new(-0.01000, 0.12200),
            new(0.00000, 0.12200),
            new(0.00000, 0.11200),
            new(0.01000, 0.11200),
            new(0.01000, 0.15200),
            new(0.00000, 0.15200),
            new(0.00000, 0.16200),
            new(-0.06000, 0.16200),
            new(-0.06000, 0.17200),
            new(-0.05000, 0.17200),
            new(-0.05000, 0.18200),
            new(-0.06000, 0.18200),
            new(-0.06000, 0.19200),
            new(-0.03000, 0.19200),
            new(-0.03000, 0.18200),
            new(-0.04000, 0.18200),
            new(-0.04000, 0.17200),
            new(0.00000, 0.17200),
            new(0.00000, 0.18200),
            new(0.01000, 0.18200),
            new(0.01000, 0.25200),
            new(0.00000, 0.25200),
            new(0.00000, 0.24200),
            new(-0.01000, 0.24200),
            new(-0.01000, 0.27200),
            new(0.00000, 0.27200),
            new(0.00000, 0.26200),
            new(0.01000, 0.26200),
            new(0.01000, 0.30200),
            new(0.00000, 0.30200),
            new(0.00000, 0.31200),
            new(-0.06000, 0.31200),
            new(-0.06000, 0.32200),
            new(-0.05000, 0.32200),
            new(-0.05000, 0.33200),
            new(-0.06000, 0.33200),
            new(-0.06000, 0.34200),
            new(-0.03000, 0.34200),
            new(-0.03000, 0.33200),
            new(-0.04000, 0.33200),
            new(-0.04000, 0.32200),
            new(0.00000, 0.32200),
            new(0.00000, 0.33200),
            new(0.01000, 0.33200),
            new(0.01000, 0.40400),
            new(0.00000, 0.40400),
            new(0.00000, 0.39400),
            new(-0.01000, 0.39400),
            new(-0.01000, 0.42400),
            new(0.00000, 0.42400),
            new(0.00000, 0.41400),
            new(0.01000, 0.41400),
            new(0.01000, 0.45400),
            new(0.00000, 0.45400),
            new(0.00000, 0.46400),
            new(-0.07000, 0.46400),
            new(-0.07000, 0.45400),
            new(-0.06000, 0.45400),
            new(-0.06000, 0.44400),
            new(-0.09000, 0.44400),
            new(-0.09000, 0.45400),
            new(-0.08000, 0.45400),
            new(-0.08000, 0.46400),
            new(-0.12000, 0.46400),
            new(-0.12000, 0.45400),
            new(-0.13000, 0.45400),
            new(-0.13000, 0.39300),
            new(-0.14000, 0.39300),
            new(-0.14000, 0.40400),
            new(-0.15000, 0.40400),
            new(-0.15000, 0.39400),
            new(-0.16000, 0.39400),
            new(-0.16000, 0.42400),
            new(-0.15000, 0.42400),
            new(-0.15000, 0.41400),
            new(-0.14000, 0.41400),
            new(-0.14000, 0.45400),
            new(-0.15000, 0.45400),
            new(-0.15000, 0.46400),
            new(-0.21000, 0.46400),
            new(-0.21000, 0.47400),
            new(-0.20000, 0.47400),
            new(-0.20000, 0.48400),
            new(-0.21000, 0.48400),
            new(-0.21000, 0.49400),
            new(-0.18000, 0.49400),
            new(-0.18000, 0.48400),
            new(-0.19000, 0.48400),
            new(-0.19000, 0.47400),
            new(-0.15000, 0.47400),
            new(-0.15000, 0.48400),
            new(-0.14000, 0.48400),
            new(-0.14000, 0.54500),
            new(-0.13000, 0.54500),
            new(-0.13000, 0.53400),
            new(-0.12000, 0.53400),
            new(-0.12000, 0.54400),
            new(-0.11000, 0.54400),
            new(-0.11000, 0.51400),
            new(-0.12000, 0.51400),
            new(-0.12000, 0.52400),
            new(-0.13000, 0.52400),
            new(-0.13000, 0.48400),
            new(-0.12000, 0.48400),
            new(-0.12000, 0.47400),
            new(-0.05000, 0.47400),
            new(-0.05000, 0.48400),
            new(-0.06000, 0.48400),
            new(-0.06000, 0.49400),
            new(-0.03000, 0.49400),
            new(-0.03000, 0.48400),
            new(-0.04000, 0.48400),
            new(-0.04000, 0.47400),
            new(0.00000, 0.47400),
            new(0.00000, 0.48400),
            new(0.01000, 0.48400),
            new(0.01000, 0.55600),
            new(0.00000, 0.55600),
            new(0.00000, 0.54600),
            new(-0.01000, 0.54600),
            new(-0.01000, 0.57600),
            new(0.00000, 0.57600),
            new(0.00000, 0.56600),
            new(0.01000, 0.56600),
            new(0.01000, 0.60600),
            new(0.00000, 0.60600),
            new(0.00000, 0.61600),
            new(-0.06000, 0.61600),
            new(-0.06000, 0.62600),
            new(-0.05000, 0.62600),
            new(-0.05000, 0.63600),
            new(-0.06000, 0.63600),
            new(-0.06000, 0.64600),
            new(-0.03000, 0.64600),
            new(-0.03000, 0.63600),
            new(-0.04000, 0.63600),
            new(-0.04000, 0.62600),
            new(0.00000, 0.62600),
            new(0.00000, 0.63600),
            new(0.01000, 0.63600),
            new(0.01000, 0.70600),
            new(0.00000, 0.70600),
            new(0.00000, 0.69600),
            new(-0.01000, 0.69600),
            new(-0.01000, 0.72600),
            new(0.00000, 0.72600),
            new(0.00000, 0.71600),
            new(0.01000, 0.71600),
            new(0.01000, 0.75600),
            new(0.00000, 0.75600),
            new(0.00000, 0.76600),
            new(-0.06000, 0.76600),
            new(-0.06000, 0.77600),
            new(-0.05000, 0.77600),
            new(-0.05000, 0.78600),
            new(-0.06000, 0.78600),
            new(-0.06000, 0.79600),
            new(-0.03000, 0.79600),
            new(-0.03000, 0.78600),
            new(-0.04000, 0.78600),
            new(-0.04000, 0.77600),
            new(0.00000, 0.77600),
            new(0.00000, 0.78600),
            new(0.01000, 0.78600),
            new(0.01000, 0.85800),
            new(0.00000, 0.85800),
            new(0.00000, 0.84800),
            new(-0.01000, 0.84800),
            new(-0.01000, 0.87800),
            new(0.00000, 0.87800),
            new(0.00000, 0.86800),
            new(0.01000, 0.86800),
            new(0.01000, 0.90800),
            new(0.00000, 0.90800),
            new(0.00000, 0.91800),
            new(-0.07000, 0.91800),
            new(-0.07000, 0.90800),
            new(-0.06000, 0.90800),
            new(-0.06000, 0.89800),
            new(-0.09000, 0.89800),
            new(-0.09000, 0.90800),
            new(-0.08000, 0.90800),
            new(-0.08000, 0.91800),
            new(-0.12000, 0.91800),
            new(-0.12000, 0.90800),
            new(-0.13000, 0.90800),
            new(-0.13000, 0.84700),
            new(-0.14000, 0.84700),
            new(-0.14000, 0.85800),
            new(-0.15000, 0.85800),
            new(-0.15000, 0.84800),
            new(-0.16000, 0.84800),
            new(-0.16000, 0.87800),
            new(-0.15000, 0.87800),
            new(-0.15000, 0.86800),
            new(-0.14000, 0.86800),
            new(-0.14000, 0.90800),
            new(-0.15000, 0.90800),
            new(-0.15000, 0.91800),
            new(-0.21000, 0.91800),
            new(-0.21000, 0.92800),
            new(-0.20000, 0.92800),
            new(-0.20000, 0.93800),
            new(-0.21000, 0.93800),
            new(-0.21000, 0.94800),
            new(-0.18000, 0.94800),
            new(-0.18000, 0.93800),
            new(-0.19000, 0.93800),
            new(-0.19000, 0.92800),
            new(-0.15000, 0.92800),
            new(-0.15000, 0.93800),
            new(-0.14000, 0.93800),
            new(-0.14000, 0.99900),
            new(-0.13000, 0.99900),
            new(-0.13000, 0.98800),
            new(-0.12000, 0.98800),
            new(-0.12000, 0.99800),
            new(-0.11000, 0.99800),
            new(-0.11000, 0.96800),
            new(-0.12000, 0.96800),
            new(-0.12000, 0.97800),
            new(-0.13000, 0.97800),
            new(-0.13000, 0.93800),
            new(-0.12000, 0.93800),
            new(-0.12000, 0.92800),
            new(-0.05000, 0.92800),
            new(-0.05000, 0.93800),
            new(-0.06000, 0.93800),
            new(-0.06000, 0.94800),
            new(-0.03000, 0.94800),
            new(-0.03000, 0.93800),
            new(-0.04000, 0.93800),
            new(-0.04000, 0.92800),
            new(0.00000, 0.92800),
            new(0.00000, 0.93800),
            new(0.01000, 0.93800),
            new(0.01000, 1.01000),
            new(0.00000, 1.01000),
            new(0.00000, 1.00000),
            new(-0.01000, 1.00000),
            new(-0.01000, 1.03000),
            new(0.00000, 1.03000),
            new(0.00000, 1.02000),
            new(0.01000, 1.02000),
            new(0.01000, 1.06000),
            new(0.00000, 1.06000),
            new(0.00000, 1.07000),
            new(-0.06000, 1.07000),
            new(-0.06000, 1.08000),
            new(-0.05000, 1.08000),
            new(-0.05000, 1.09000),
            new(-0.06000, 1.09000),
            new(-0.06000, 1.10000),
            new(-0.03000, 1.10000),
            new(-0.03000, 1.09000),
            new(-0.04000, 1.09000),
            new(-0.04000, 1.08000),
            new(0.00000, 1.08000),
            new(0.00000, 1.09000),
            new(0.01000, 1.09000),
            new(0.01000, 1.16000),
            new(0.00000, 1.16000),
            new(0.00000, 1.15000),
            new(-0.01000, 1.15000),
            new(-0.01000, 1.18000),
            new(0.00000, 1.18000),
            new(0.00000, 1.17000),
            new(0.01000, 1.17000),
            new(0.01000, 1.21000),
            new(0.00000, 1.21000),
            new(0.00000, 1.22000),
            new(-0.06000, 1.22000),
            new(-0.06000, 1.23000),
            new(-0.05000, 1.23000),
            new(-0.05000, 1.24000),
            new(-0.06000, 1.24000),
            new(-0.06000, 1.25000),
            new(-0.03000, 1.25000),
            new(-0.03000, 1.24000),
            new(-0.04000, 1.24000),
            new(-0.04000, 1.23000),
            new(0.00000, 1.23000),
            new(0.00000, 1.24000),
            new(0.01000, 1.24000),
            new(0.01000, 1.31200),
            new(0.00000, 1.31200),
            new(0.00000, 1.30200),
            new(-0.01000, 1.30200),
            new(-0.01000, 1.33200),
            new(0.00000, 1.33200),
            new(0.00000, 1.32200),
            new(0.01000, 1.32200),
            new(0.01000, 1.36200),
            new(0.00000, 1.36200),
            new(0.00000, 1.37200),
            new(-0.07000, 1.37200),
            new(-0.07000, 1.36200),
            new(-0.06000, 1.36200),
            new(-0.06000, 1.35200),
            new(-0.09000, 1.35200),
            new(-0.09000, 1.36200),
            new(-0.08000, 1.36200),
            new(-0.08000, 1.37200),
            new(-0.12000, 1.37200),
            new(-0.12000, 1.36200),
            new(-0.13000, 1.36200),
            new(-0.13000, 1.30100),
            new(-0.14000, 1.30100),
            new(-0.14000, 1.31200),
            new(-0.15000, 1.31200),
            new(-0.15000, 1.30200),
            new(-0.16000, 1.30200),
            new(-0.16000, 1.33200),
            new(-0.15000, 1.33200),
            new(-0.15000, 1.32200),
            new(-0.14000, 1.32200),
            new(-0.14000, 1.36200),
            new(-0.15000, 1.36200),
            new(-0.15000, 1.37200),
            new(-0.21000, 1.37200),
            new(-0.21000, 1.38200),
            new(-0.20000, 1.38200),
            new(-0.20000, 1.39200),
            new(-0.21000, 1.39200),
            new(-0.21000, 1.40200),
            new(-0.18000, 1.40200),
            new(-0.18000, 1.39200),
            new(-0.19000, 1.39200),
            new(-0.19000, 1.38200),
            new(-0.15000, 1.38200),
            new(-0.15000, 1.39200),
            new(-0.14000, 1.39200),
            new(-0.14000, 1.45300),
            new(-0.13000, 1.45300),
            new(-0.13000, 1.44200),
            new(-0.12000, 1.44200),
            new(-0.12000, 1.45200),
            new(-0.11000, 1.45200),
            new(-0.11000, 1.42200),
            new(-0.12000, 1.42200),
            new(-0.12000, 1.43200),
            new(-0.13000, 1.43200),
            new(-0.13000, 1.39200),
            new(-0.12000, 1.39200),
            new(-0.12000, 1.38200),
            new(-0.05000, 1.38200),
            new(-0.05000, 1.39200),
            new(-0.06000, 1.39200),
            new(-0.06000, 1.40200),
            new(-0.03000, 1.40200),
            new(-0.03000, 1.39200),
            new(-0.04000, 1.39200),
            new(-0.04000, 1.38200),
            new(0.00000, 1.38200),
            new(0.00000, 1.39200),
            new(0.01000, 1.39200),
            new(0.01000, 1.46400),
            new(0.00000, 1.46400),
            new(0.00000, 1.45400),
            new(-0.01000, 1.45400),
            new(-0.01000, 1.48400),
            new(0.00000, 1.48400),
            new(0.00000, 1.47400),
            new(0.01000, 1.47400),
            new(0.01000, 1.51400),
            new(0.00000, 1.51400),
            new(0.00000, 1.52400),
            new(-0.06000, 1.52400),
            new(-0.06000, 1.53400),
            new(-0.05000, 1.53400),
            new(-0.05000, 1.54400),
            new(-0.06000, 1.54400),
            new(-0.06000, 1.55400),
            new(-0.03000, 1.55400),
            new(-0.03000, 1.54400),
            new(-0.04000, 1.54400),
            new(-0.04000, 1.53400),
            new(0.00000, 1.53400),
            new(0.00000, 1.54400),
            new(0.01000, 1.54400),
            new(0.01000, 1.61400),
            new(0.00000, 1.61400),
            new(0.00000, 1.60400),
            new(-0.01000, 1.60400),
            new(-0.01000, 1.63400),
            new(0.00000, 1.63400),
            new(0.00000, 1.62400),
            new(0.01000, 1.62400),
            new(0.01000, 1.66400),
            new(0.00000, 1.66400),
            new(0.00000, 1.67400),
            new(-0.06000, 1.67400),
            new(-0.06000, 1.68400),
            new(-0.05000, 1.68400),
            new(-0.05000, 1.69400),
            new(-0.06000, 1.69400),
            new(-0.06000, 1.70400),
            new(-0.03000, 1.70400),
            new(-0.03000, 1.69400),
            new(-0.04000, 1.69400),
            new(-0.04000, 1.68400),
            new(0.00000, 1.68400),
            new(0.00000, 1.69400),
            new(0.01000, 1.69400),
            new(0.01000, 1.76600),
            new(0.00000, 1.76600),
            new(0.00000, 1.75600),
            new(-0.01000, 1.75600),
            new(-0.01000, 1.78600),
            new(0.00000, 1.78600),
            new(0.00000, 1.77600),
            new(0.01000, 1.77600),
            new(0.01000, 1.81600),
            new(0.00000, 1.81600),
            new(0.00000, 1.82600),
            new(-0.07000, 1.82600),
            new(-0.07000, 1.81600),
            new(-0.06000, 1.81600),
            new(-0.06000, 1.80600),
            new(-0.09000, 1.80600),
            new(-0.09000, 1.81600),
            new(-0.08000, 1.81600),
            new(-0.08000, 1.82600),
            new(-0.12000, 1.82600),
            new(-0.12000, 1.81600),
            new(-0.13000, 1.81600),
            new(-0.13000, 1.75500),
            new(-0.14000, 1.75500),
            new(-0.14000, 1.76600),
            new(-0.15000, 1.76600),
            new(-0.15000, 1.75600),
            new(-0.16000, 1.75600),
            new(-0.16000, 1.78600),
            new(-0.15000, 1.78600),
            new(-0.15000, 1.77600),
            new(-0.14000, 1.77600),
            new(-0.14000, 1.81600),
            new(-0.15000, 1.81600),
            new(-0.15000, 1.82600),
            new(-0.21000, 1.82600),
            new(-0.21000, 1.83600),
            new(-0.20000, 1.83600),
            new(-0.20000, 1.84600),
            new(-0.21000, 1.84600),
            new(-0.21000, 1.85600),
            new(-0.18000, 1.85600),
            new(-0.18000, 1.84600),
            new(-0.19000, 1.84600),
            new(-0.19000, 1.83600),
            new(-0.15000, 1.83600),
            new(-0.15000, 1.84600),
            new(-0.14000, 1.84600),
            new(-0.14000, 1.90700),
            new(-0.13000, 1.90700),
            new(-0.13000, 1.89600),
            new(-0.12000, 1.89600),
            new(-0.12000, 1.90600),
            new(-0.11000, 1.90600),
            new(-0.11000, 1.87600),
            new(-0.12000, 1.87600),
            new(-0.12000, 1.88600),
            new(-0.13000, 1.88600),
            new(-0.13000, 1.84600),
            new(-0.12000, 1.84600),
            new(-0.12000, 1.83600),
            new(-0.05000, 1.83600),
            new(-0.05000, 1.84600),
            new(-0.06000, 1.84600),
            new(-0.06000, 1.85600),
            new(-0.03000, 1.85600),
            new(-0.03000, 1.84600),
            new(-0.04000, 1.84600),
            new(-0.04000, 1.83600),
            new(0.00000, 1.83600),
            new(0.00000, 1.84600),
            new(0.01000, 1.84600),
            new(0.01000, 1.91800),
            new(0.00000, 1.91800),
            new(0.00000, 1.90800),
            new(-0.01000, 1.90800),
            new(-0.01000, 1.93800),
            new(0.00000, 1.93800),
            new(0.00000, 1.92800),
            new(0.01000, 1.92800),
            new(0.01000, 1.96800),
            new(0.00000, 1.96800),
            new(0.00000, 1.97800),
            new(-0.06000, 1.97800),
            new(-0.06000, 1.98800),
            new(-0.05000, 1.98800),
            new(-0.05000, 1.99800),
            new(-0.06000, 1.99800),
            new(-0.06000, 2.00800),
            new(-0.03000, 2.00800),
            new(-0.03000, 1.99800),
            new(-0.04000, 1.99800),
            new(-0.04000, 1.98800),
            new(0.00000, 1.98800),
            new(0.00000, 1.99800),
            new(0.01000, 1.99800),
            new(0.01000, 2.05900),
            new(0.02000, 2.05900),
            new(0.02000, 2.04800),
            new(0.03000, 2.04800),
            new(0.03000, 2.05800),
            new(0.04000, 2.05800),
            new(0.04000, 2.02800),
            new(0.03000, 2.02800),
            new(0.03000, 2.03800),
            new(0.02000, 2.03800),
            new(0.02000, 1.99800),
            new(0.03000, 1.99800),
            new(0.03000, 1.98800),
            new(0.09000, 1.98800),
            new(0.09000, 1.97800),
            new(0.08000, 1.97800),
            new(0.08000, 1.96800),
            new(0.09000, 1.96800),
            new(0.09000, 1.95800),
            new(0.06000, 1.95800),
            new(0.06000, 1.96800),
            new(0.07000, 1.96800),
            new(0.07000, 1.97800),
            new(0.03000, 1.97800),
            new(0.03000, 1.96800),
            new(0.02000, 1.96800),
            new(0.02000, 1.89600),
            new(0.03000, 1.89600),
            new(0.03000, 1.90600),
            new(0.04000, 1.90600),
            new(0.04000, 1.87600),
            new(0.03000, 1.87600),
            new(0.03000, 1.88600),
            new(0.02000, 1.88600),
            new(0.02000, 1.84600),
            new(0.03000, 1.84600),
            new(0.03000, 1.83600),
            new(0.10000, 1.83600),
            new(0.10000, 1.84600),
            new(0.09000, 1.84600),
            new(0.09000, 1.85600),
            new(0.12000, 1.85600),
            new(0.12000, 1.84600),
            new(0.11000, 1.84600),
            new(0.11000, 1.83600),
            new(0.15000, 1.83600),
            new(0.15000, 1.84600),
            new(0.16000, 1.84600),
            new(0.16000, 1.90700),
            new(0.17000, 1.90700),
            new(0.17000, 1.89600),
            new(0.18000, 1.89600),
            new(0.18000, 1.90600),
            new(0.19000, 1.90600),
            new(0.19000, 1.87600),
            new(0.18000, 1.87600),
            new(0.18000, 1.88600),
            new(0.17000, 1.88600),
            new(0.17000, 1.84600),
            new(0.18000, 1.84600),
            new(0.18000, 1.83600),
            new(0.24000, 1.83600),
            new(0.24000, 1.82600),
            new(0.23000, 1.82600),
            new(0.23000, 1.81600),
            new(0.24000, 1.81600),
            new(0.24000, 1.80600),
            new(0.21000, 1.80600),
            new(0.21000, 1.81600),
            new(0.22000, 1.81600),
            new(0.22000, 1.82600),
            new(0.18000, 1.82600),
            new(0.18000, 1.81600),
            new(0.17000, 1.81600),
            new(0.17000, 1.75500),
            new(0.16000, 1.75500),
            new(0.16000, 1.76600),
            new(0.15000, 1.76600),
            new(0.15000, 1.75600),
            new(0.14000, 1.75600),
            new(0.14000, 1.78600),
            new(0.15000, 1.78600),
            new(0.15000, 1.77600),
            new(0.16000, 1.77600),
            new(0.16000, 1.81600),
            new(0.15000, 1.81600),
            new(0.15000, 1.82600),
            new(0.08000, 1.82600),
            new(0.08000, 1.81600),
            new(0.09000, 1.81600),
            new(0.09000, 1.80600),
            new(0.06000, 1.80600),
            new(0.06000, 1.81600),
            new(0.07000, 1.81600),
            new(0.07000, 1.82600),
            new(0.03000, 1.82600),
            new(0.03000, 1.81600),
            new(0.02000, 1.81600),
            new(0.02000, 1.74400),
            new(0.03000, 1.74400),
            new(0.03000, 1.75400),
            new(0.04000, 1.75400),
            new(0.04000, 1.72400),
            new(0.03000, 1.72400),
            new(0.03000, 1.73400),
            new(0.02000, 1.73400),
            new(0.02000, 1.69400),
            new(0.03000, 1.69400),
            new(0.03000, 1.68400),
            new(0.09000, 1.68400),
            new(0.09000, 1.67400),
            new(0.08000, 1.67400),
            new(0.08000, 1.66400),
            new(0.09000, 1.66400),
            new(0.09000, 1.65400),
            new(0.06000, 1.65400),
            new(0.06000, 1.66400),
            new(0.07000, 1.66400),
            new(0.07000, 1.67400),
            new(0.03000, 1.67400),
            new(0.03000, 1.66400),
            new(0.02000, 1.66400),
            new(0.02000, 1.59400),
            new(0.03000, 1.59400),
            new(0.03000, 1.60400),
            new(0.04000, 1.60400),
            new(0.04000, 1.57400),
            new(0.03000, 1.57400),
            new(0.03000, 1.58400),
            new(0.02000, 1.58400),
            new(0.02000, 1.54400),
            new(0.03000, 1.54400),
            new(0.03000, 1.53400),
            new(0.09000, 1.53400),
            new(0.09000, 1.52400),
            new(0.08000, 1.52400),
            new(0.08000, 1.51400),
            new(0.09000, 1.51400),
            new(0.09000, 1.50400),
            new(0.06000, 1.50400),
            new(0.06000, 1.51400),
            new(0.07000, 1.51400),
            new(0.07000, 1.52400),
            new(0.03000, 1.52400),
            new(0.03000, 1.51400),
            new(0.02000, 1.51400),
            new(0.02000, 1.44200),
            new(0.03000, 1.44200),
            new(0.03000, 1.45200),
            new(0.04000, 1.45200),
            new(0.04000, 1.42200),
            new(0.03000, 1.42200),
            new(0.03000, 1.43200),
            new(0.02000, 1.43200),
            new(0.02000, 1.39200),
            new(0.03000, 1.39200),
            new(0.03000, 1.38200),
            new(0.10000, 1.38200),
            new(0.10000, 1.39200),
            new(0.09000, 1.39200),
            new(0.09000, 1.40200),
            new(0.12000, 1.40200),
            new(0.12000, 1.39200),
            new(0.11000, 1.39200),
            new(0.11000, 1.38200),
            new(0.15000, 1.38200),
            new(0.15000, 1.39200),
            new(0.16000, 1.39200),
            new(0.16000, 1.45300),
            new(0.17000, 1.45300),
            new(0.17000, 1.44200),
            new(0.18000, 1.44200),
            new(0.18000, 1.45200),
            new(0.19000, 1.45200),
            new(0.19000, 1.42200),
            new(0.18000, 1.42200),
            new(0.18000, 1.43200),
            new(0.17000, 1.43200),
            new(0.17000, 1.39200),
            new(0.18000, 1.39200),
            new(0.18000, 1.38200),
            new(0.24000, 1.38200),
            new(0.24000, 1.37200),
            new(0.23000, 1.37200),
            new(0.23000, 1.36200),
            new(0.24000, 1.36200),
            new(0.24000, 1.35200),
            new(0.21000, 1.35200),
            new(0.21000, 1.36200),
            new(0.22000, 1.36200),
            new(0.22000, 1.37200),
            new(0.18000, 1.37200),
            new(0.18000, 1.36200),
            new(0.17000, 1.36200),
            new(0.17000, 1.30100),
            new(0.16000, 1.30100),
            new(0.16000, 1.31200),
            new(0.15000, 1.31200),
            new(0.15000, 1.30200),
            new(0.14000, 1.30200),
            new(0.14000, 1.33200),
            new(0.15000, 1.33200),
            new(0.15000, 1.32200),
            new(0.16000, 1.32200),
            new(0.16000, 1.36200),
            new(0.15000, 1.36200),
            new(0.15000, 1.37200),
            new(0.08000, 1.37200),
            new(0.08000, 1.36200),
            new(0.09000, 1.36200),
            new(0.09000, 1.35200),
            new(0.06000, 1.35200),
            new(0.06000, 1.36200),
            new(0.07000, 1.36200),
            new(0.07000, 1.37200),
            new(0.03000, 1.37200),
            new(0.03000, 1.36200),
            new(0.02000, 1.36200),
            new(0.02000, 1.29000),
            new(0.03000, 1.29000),
            new(0.03000, 1.30000),
            new(0.04000, 1.30000),
            new(0.04000, 1.27000),
            new(0.03000, 1.27000),
            new(0.03000, 1.28000),
            new(0.02000, 1.28000),
            new(0.02000, 1.24000),
            new(0.03000, 1.24000),
            new(0.03000, 1.23000),
            new(0.09000, 1.23000),
            new(0.09000, 1.22000),
            new(0.08000, 1.22000),
            new(0.08000, 1.21000),
            new(0.09000, 1.21000),
            new(0.09000, 1.20000),
            new(0.06000, 1.20000),
            new(0.06000, 1.21000),
            new(0.07000, 1.21000),
            new(0.07000, 1.22000),
            new(0.03000, 1.22000),
            new(0.03000, 1.21000),
            new(0.02000, 1.21000),
            new(0.02000, 1.14000),
            new(0.03000, 1.14000),
            new(0.03000, 1.15000),
            new(0.04000, 1.15000),
            new(0.04000, 1.12000),
            new(0.03000, 1.12000),
            new(0.03000, 1.13000),
            new(0.02000, 1.13000),
            new(0.02000, 1.09000),
            new(0.03000, 1.09000),
            new(0.03000, 1.08000),
            new(0.09000, 1.08000),
            new(0.09000, 1.07000),
            new(0.08000, 1.07000),
            new(0.08000, 1.06000),
            new(0.09000, 1.06000),
            new(0.09000, 1.05000),
            new(0.06000, 1.05000),
            new(0.06000, 1.06000),
            new(0.07000, 1.06000),
            new(0.07000, 1.07000),
            new(0.03000, 1.07000),
            new(0.03000, 1.06000),
            new(0.02000, 1.06000),
            new(0.02000, 0.98800),
            new(0.03000, 0.98800),
            new(0.03000, 0.99800),
            new(0.04000, 0.99800),
            new(0.04000, 0.96800),
            new(0.03000, 0.96800),
            new(0.03000, 0.97800),
            new(0.02000, 0.97800),
            new(0.02000, 0.93800),
            new(0.03000, 0.93800),
            new(0.03000, 0.92800),
            new(0.10000, 0.92800),
            new(0.10000, 0.93800),
            new(0.09000, 0.93800),
            new(0.09000, 0.94800),
            new(0.12000, 0.94800),
            new(0.12000, 0.93800),
            new(0.11000, 0.93800),
            new(0.11000, 0.92800),
            new(0.15000, 0.92800),
            new(0.15000, 0.93800),
            new(0.16000, 0.93800),
            new(0.16000, 0.99900),
            new(0.17000, 0.99900),
            new(0.17000, 0.98800),
            new(0.18000, 0.98800),
            new(0.18000, 0.99800),
            new(0.19000, 0.99800),
            new(0.19000, 0.96800),
            new(0.18000, 0.96800),
            new(0.18000, 0.97800),
            new(0.17000, 0.97800),
            new(0.17000, 0.93800),
            new(0.18000, 0.93800),
            new(0.18000, 0.92800),
            new(0.24000, 0.92800),
            new(0.24000, 0.91800),
            new(0.23000, 0.91800),
            new(0.23000, 0.90800),
            new(0.24000, 0.90800),
            new(0.24000, 0.89800),
            new(0.21000, 0.89800),
            new(0.21000, 0.90800),
            new(0.22000, 0.90800),
            new(0.22000, 0.91800),
            new(0.18000, 0.91800),
            new(0.18000, 0.90800),
            new(0.17000, 0.90800),
            new(0.17000, 0.84700),
            new(0.16000, 0.84700),
            new(0.16000, 0.85800),
            new(0.15000, 0.85800),
            new(0.15000, 0.84800),
            new(0.14000, 0.84800),
            new(0.14000, 0.87800),
            new(0.15000, 0.87800),
            new(0.15000, 0.86800),
            new(0.16000, 0.86800),
            new(0.16000, 0.90800),
            new(0.15000, 0.90800),
            new(0.15000, 0.91800),
            new(0.08000, 0.91800),
            new(0.08000, 0.90800),
            new(0.09000, 0.90800),
            new(0.09000, 0.89800),
            new(0.06000, 0.89800),
            new(0.06000, 0.90800),
            new(0.07000, 0.90800),
            new(0.07000, 0.91800),
            new(0.03000, 0.91800),
            new(0.03000, 0.90800),
            new(0.02000, 0.90800),
            new(0.02000, 0.83600),
            new(0.03000, 0.83600),
            new(0.03000, 0.84600),
            new(0.04000, 0.84600),
            new(0.04000, 0.81600),
            new(0.03000, 0.81600),
            new(0.03000, 0.82600),
            new(0.02000, 0.82600),
            new(0.02000, 0.78600),
            new(0.03000, 0.78600),
            new(0.03000, 0.77600),
            new(0.09000, 0.77600),
            new(0.09000, 0.76600),
            new(0.08000, 0.76600),
            new(0.08000, 0.75600),
            new(0.09000, 0.75600),
            new(0.09000, 0.74600),
            new(0.06000, 0.74600),
            new(0.06000, 0.75600),
            new(0.07000, 0.75600),
            new(0.07000, 0.76600),
            new(0.03000, 0.76600),
            new(0.03000, 0.75600),
            new(0.02000, 0.75600),
            new(0.02000, 0.68600),
            new(0.03000, 0.68600),
            new(0.03000, 0.69600),
            new(0.04000, 0.69600),
            new(0.04000, 0.66600),
            new(0.03000, 0.66600),
            new(0.03000, 0.67600),
            new(0.02000, 0.67600),
            new(0.02000, 0.63600),
            new(0.03000, 0.63600),
            new(0.03000, 0.62600),
            new(0.09000, 0.62600),
            new(0.09000, 0.61600),
            new(0.08000, 0.61600),
            new(0.08000, 0.60600),
            new(0.09000, 0.60600),
            new(0.09000, 0.59600),
            new(0.06000, 0.59600),
            new(0.06000, 0.60600),
            new(0.07000, 0.60600),
            new(0.07000, 0.61600),
            new(0.03000, 0.61600),
            new(0.03000, 0.60600),
            new(0.02000, 0.60600),
            new(0.02000, 0.53400),
            new(0.03000, 0.53400),
            new(0.03000, 0.54400),
            new(0.04000, 0.54400),
            new(0.04000, 0.51400),
            new(0.03000, 0.51400),
            new(0.03000, 0.52400),
            new(0.02000, 0.52400),
            new(0.02000, 0.48400),
            new(0.03000, 0.48400),
            new(0.03000, 0.47400),
            new(0.10000, 0.47400),
            new(0.10000, 0.48400),
            new(0.09000, 0.48400),
            new(0.09000, 0.49400),
            new(0.12000, 0.49400),
            new(0.12000, 0.48400),
            new(0.11000, 0.48400),
            new(0.11000, 0.47400),
            new(0.15000, 0.47400),
            new(0.15000, 0.48400),
            new(0.16000, 0.48400),
            new(0.16000, 0.54500),
            new(0.17000, 0.54500),
            new(0.17000, 0.53400),
            new(0.18000, 0.53400),
            new(0.18000, 0.54400),
            new(0.19000, 0.54400),
            new(0.19000, 0.51400),
            new(0.18000, 0.51400),
            new(0.18000, 0.52400),
            new(0.17000, 0.52400),
            new(0.17000, 0.48400),
            new(0.18000, 0.48400),
            new(0.18000, 0.47400),
            new(0.24000, 0.47400),
            new(0.24000, 0.46400),
            new(0.23000, 0.46400),
            new(0.23000, 0.45400),
            new(0.24000, 0.45400),
            new(0.24000, 0.44400),
            new(0.21000, 0.44400),
            new(0.21000, 0.45400),
            new(0.22000, 0.45400),
            new(0.22000, 0.46400),
            new(0.18000, 0.46400),
            new(0.18000, 0.45400),
            new(0.17000, 0.45400),
            new(0.17000, 0.39300),
            new(0.16000, 0.39300),
            new(0.16000, 0.40400),
            new(0.15000, 0.40400),
            new(0.15000, 0.39400),
            new(0.14000, 0.39400),
            new(0.14000, 0.42400),
            new(0.15000, 0.42400),
            new(0.15000, 0.41400),
            new(0.16000, 0.41400),
            new(0.16000, 0.45400),
            new(0.15000, 0.45400),
            new(0.15000, 0.46400),
            new(0.08000, 0.46400),
            new(0.08000, 0.45400),
            new(0.09000, 0.45400),
            new(0.09000, 0.44400),
            new(0.06000, 0.44400),
            new(0.06000, 0.45400),
            new(0.07000, 0.45400),
            new(0.07000, 0.46400),
            new(0.03000, 0.46400),
            new(0.03000, 0.45400),
            new(0.02000, 0.45400),
            new(0.02000, 0.38200),
            new(0.03000, 0.38200),
            new(0.03000, 0.39200),
            new(0.04000, 0.39200),
            new(0.04000, 0.36200),
            new(0.03000, 0.36200),
            new(0.03000, 0.37200),
            new(0.02000, 0.37200),
            new(0.02000, 0.33200),
            new(0.03000, 0.33200),
            new(0.03000, 0.32200),
            new(0.09000, 0.32200),
            new(0.09000, 0.31200),
            new(0.08000, 0.31200),
            new(0.08000, 0.30200),
            new(0.09000, 0.30200),
            new(0.09000, 0.29200),
            new(0.06000, 0.29200),
            new(0.06000, 0.30200),
            new(0.07000, 0.30200),
            new(0.07000, 0.31200),
            new(0.03000, 0.31200),
            new(0.03000, 0.30200),
            new(0.02000, 0.30200),
            new(0.02000, 0.23200),
            new(0.03000, 0.23200),
            new(0.03000, 0.24200),
            new(0.04000, 0.24200),
            new(0.04000, 0.21200),
            new(0.03000, 0.21200),
            new(0.03000, 0.22200),
            new(0.02000, 0.22200),
            new(0.02000, 0.18200),
            new(0.03000, 0.18200),
            new(0.03000, 0.17200),
            new(0.09000, 0.17200),
            new(0.09000, 0.16200),
            new(0.08000, 0.16200),
            new(0.08000, 0.15200),
            new(0.09000, 0.15200),
            new(0.09000, 0.14200),
            new(0.06000, 0.14200),
            new(0.06000, 0.15200),
            new(0.07000, 0.15200),
            new(0.07000, 0.16200),
            new(0.03000, 0.16200),
            new(0.03000, 0.15200),
            new(0.02000, 0.15200),
            new(0.02000, 0.08000),
            new(0.03000, 0.08000),
            new(0.03000, 0.09000),
            new(0.04000, 0.09000),
            new(0.04000, 0.06000),
            new(0.03000, 0.06000),
            new(0.03000, 0.07000),
            new(0.02000, 0.07000),
            new(0.02000, 0.03000),
            new(0.03000, 0.03000),
            new(0.03000, 0.02000),
            new(0.10000, 0.02000),
            new(0.10000, 0.03000),
            new(0.09000, 0.03000),
            new(0.09000, 0.04000),
            new(0.12000, 0.04000),
            new(0.12000, 0.03000),
            new(0.11000, 0.03000),
            new(0.11000, 0.02000),
            new(0.15000, 0.02000),
            new(0.15000, 0.03000),
            new(0.16000, 0.03000),
            new(0.16000, 0.09100),
            new(0.17000, 0.09100),
            new(0.17000, 0.08000),
            new(0.18000, 0.08000),
            new(0.18000, 0.09000),
            new(0.19000, 0.09000),
            new(0.19000, 0.06000),
            new(0.18000, 0.06000),
            new(0.18000, 0.07000),
            new(0.17000, 0.07000),
            new(0.17000, 0.03000),
            new(0.18000, 0.03000),
            new(0.18000, 0.02000),
            new(0.25000, 0.02000),
            new(0.25000, 0.03000),
            new(0.24000, 0.03000),
            new(0.24000, 0.04000),
            new(0.27000, 0.04000),
            new(0.27000, 0.03000),
            new(0.26000, 0.03000),
            new(0.26000, 0.02000),
            new(0.30000, 0.02000),
            new(0.30000, 0.03000),
            new(0.31000, 0.03000),
            new(0.31000, 0.09100),
            new(0.32000, 0.09100),
            new(0.32000, 0.08000),
            new(0.33000, 0.08000),
            new(0.33000, 0.09000),
            new(0.34000, 0.09000),
            new(0.34000, 0.06000),
            new(0.33000, 0.06000),
            new(0.33000, 0.07000),
            new(0.32000, 0.07000),
            new(0.32000, 0.03000),
            new(0.33000, 0.03000),
            new(0.33000, 0.02000),
            new(0.40000, 0.02000),
            new(0.40000, 0.03000),
            new(0.39000, 0.03000),
            new(0.39000, 0.04000),
            new(0.42000, 0.04000),
            new(0.42000, 0.03000),
            new(0.41000, 0.03000),
            new(0.41000, 0.02000),
            new(0.45000, 0.02000),
            new(0.45000, 0.03000),
            new(0.46000, 0.03000),
            new(0.46000, 0.10200),
            new(0.45000, 0.10200),
            new(0.45000, 0.09200),
            new(0.44000, 0.09200),
            new(0.44000, 0.12200),
            new(0.45000, 0.12200),
            new(0.45000, 0.11200),
            new(0.46000, 0.11200),
            new(0.46000, 0.15200),
            new(0.45000, 0.15200),
            new(0.45000, 0.16200),
            new(0.39000, 0.16200),
            new(0.39000, 0.17200),
            new(0.40000, 0.17200),
            new(0.40000, 0.18200),
            new(0.39000, 0.18200),
            new(0.39000, 0.19200),
            new(0.42000, 0.19200),
            new(0.42000, 0.18200),
            new(0.41000, 0.18200),
            new(0.41000, 0.17200),
            new(0.45000, 0.17200),
            new(0.45000, 0.18200),
            new(0.46000, 0.18200),
            new(0.46000, 0.24300),
            new(0.47000, 0.24300),
            new(0.47000, 0.23200),
            new(0.48000, 0.23200),
            new(0.48000, 0.24200),
            new(0.49000, 0.24200),
            new(0.49000, 0.21200),
            new(0.48000, 0.21200),
            new(0.48000, 0.22200),
            new(0.47000, 0.22200),
            new(0.47000, 0.18200),
            new(0.48000, 0.18200),
            new(0.48000, 0.17200),
            new(0.54000, 0.17200),
            new(0.54000, 0.16200),
            new(0.53000, 0.16200),
            new(0.53000, 0.15200),
            new(0.54000, 0.15200),
            new(0.54000, 0.14200),
            new(0.51000, 0.14200),
            new(0.51000, 0.15200),
            new(0.52000, 0.15200),
            new(0.52000, 0.16200),
            new(0.48000, 0.16200),
            new(0.48000, 0.15200),
            new(0.47000, 0.15200),
            new(0.47000, 0.08000),
            new(0.48000, 0.08000),
            new(0.48000, 0.09000),
            new(0.49000, 0.09000),
            new(0.49000, 0.06000),
            new(0.48000, 0.06000),
            new(0.48000, 0.07000),
            new(0.47000, 0.07000),
            new(0.47000, 0.03000),
            new(0.48000, 0.03000),
            new(0.48000, 0.02000),
            new(0.55000, 0.02000),
            new(0.55000, 0.03000),
            new(0.54000, 0.03000),
            new(0.54000, 0.04000),
            new(0.57000, 0.04000),
            new(0.57000, 0.03000),
            new(0.56000, 0.03000),
            new(0.56000, 0.02000),
            new(0.60000, 0.02000),
            new(0.60000, 0.03000),
            new(0.61000, 0.03000),
            new(0.61000, 0.09100),
            new(0.62000, 0.09100),
            new(0.62000, 0.08000),
            new(0.63000, 0.08000),
            new(0.63000, 0.09000),
            new(0.64000, 0.09000),
            new(0.64000, 0.06000),
            new(0.63000, 0.06000),
            new(0.63000, 0.07000),
            new(0.62000, 0.07000),
            new(0.62000, 0.03000),
            new(0.63000, 0.03000),
            new(0.63000, 0.02000),
            new(0.70000, 0.02000),
            new(0.70000, 0.03000),
            new(0.69000, 0.03000),
            new(0.69000, 0.04000),
            new(0.72000, 0.04000),
            new(0.72000, 0.03000),
            new(0.71000, 0.03000),
            new(0.71000, 0.02000),
            new(0.75000, 0.02000),
            new(0.75000, 0.03000),
            new(0.76000, 0.03000),
            new(0.76000, 0.09100),
            new(0.77000, 0.09100),
            new(0.77000, 0.08000),
            new(0.78000, 0.08000),
            new(0.78000, 0.09000),
            new(0.79000, 0.09000),
            new(0.79000, 0.06000),
            new(0.78000, 0.06000),
            new(0.78000, 0.07000),
            new(0.77000, 0.07000),
            new(0.77000, 0.03000),
            new(0.78000, 0.03000),
            new(0.78000, 0.02000),
            new(0.85000, 0.02000),
            new(0.85000, 0.03000),
            new(0.84000, 0.03000),
            new(0.84000, 0.04000),
            new(0.87000, 0.04000),
            new(0.87000, 0.03000),
            new(0.86000, 0.03000),
            new(0.86000, 0.02000),
            new(0.90000, 0.02000),
            new(0.90000, 0.03000),
            new(0.91000, 0.03000),
            new(0.91000, 0.10200),
            new(0.90000, 0.10200),
            new(0.90000, 0.09200),
            new(0.89000, 0.09200),
            new(0.89000, 0.12200),
            new(0.90000, 0.12200),
            new(0.90000, 0.11200),
            new(0.91000, 0.11200),
            new(0.91000, 0.15200),
            new(0.90000, 0.15200),
            new(0.90000, 0.16200),
            new(0.84000, 0.16200),
            new(0.84000, 0.17200),
            new(0.85000, 0.17200),
            new(0.85000, 0.18200),
            new(0.84000, 0.18200),
            new(0.84000, 0.19200),
            new(0.87000, 0.19200),
            new(0.87000, 0.18200),
            new(0.86000, 0.18200),
            new(0.86000, 0.17200),
            new(0.90000, 0.17200),
            new(0.90000, 0.18200),
            new(0.91000, 0.18200),
            new(0.91000, 0.24300),
            new(0.92000, 0.24300),
            new(0.92000, 0.23200),
            new(0.93000, 0.23200),
            new(0.93000, 0.24200),
            new(0.94000, 0.24200),
            new(0.94000, 0.21200),
            new(0.93000, 0.21200),
            new(0.93000, 0.22200),
            new(0.92000, 0.22200),
            new(0.92000, 0.18200),
            new(0.93000, 0.18200),
            new(0.93000, 0.17200),
            new(0.99000, 0.17200),
            new(0.99000, 0.16200),
            new(0.98000, 0.16200),
            new(0.98000, 0.15200),
            new(0.99000, 0.15200),
            new(0.99000, 0.14200),
            new(0.96000, 0.14200),
            new(0.96000, 0.15200),
            new(0.97000, 0.15200),
            new(0.97000, 0.16200),
            new(0.93000, 0.16200),
            new(0.93000, 0.15200),
            new(0.92000, 0.15200),
            new(0.92000, 0.08000),
            new(0.93000, 0.08000),
            new(0.93000, 0.09000),
            new(0.94000, 0.09000),
            new(0.94000, 0.06000),
            new(0.93000, 0.06000),
            new(0.93000, 0.07000),
            new(0.92000, 0.07000),
            new(0.92000, 0.03000),
            new(0.93000, 0.03000),
            new(0.93000, 0.02000),
            new(1.00000, 0.02000),
            new(1.00000, 0.03000),
            new(0.99000, 0.03000),
            new(0.99000, 0.04000),
            new(1.02000, 0.04000),
            new(1.02000, 0.03000),
            new(1.01000, 0.03000),
            new(1.01000, 0.02000),
            new(1.05000, 0.02000),
            new(1.05000, 0.03000),
            new(1.06000, 0.03000),
            new(1.06000, 0.09100),
            new(1.07000, 0.09100),
            new(1.07000, 0.08000),
            new(1.08000, 0.08000),
            new(1.08000, 0.09000),
            new(1.09000, 0.09000),
            new(1.09000, 0.06000),
            new(1.08000, 0.06000),
            new(1.08000, 0.07000),
            new(1.07000, 0.07000),
            new(1.07000, 0.03000),
            new(1.08000, 0.03000),
            new(1.08000, 0.02000),
            new(1.15000, 0.02000),
            new(1.15000, 0.03000),
            new(1.14000, 0.03000),
            new(1.14000, 0.04000),
            new(1.17000, 0.04000),
            new(1.17000, 0.03000),
            new(1.16000, 0.03000),
            new(1.16000, 0.02000),
            new(1.20000, 0.02000),
            new(1.20000, 0.03000),
            new(1.21000, 0.03000),
            new(1.21000, 0.09100),
            new(1.22000, 0.09100),
            new(1.22000, 0.08000),
            new(1.23000, 0.08000),
            new(1.23000, 0.09000),
            new(1.24000, 0.09000),
            new(1.24000, 0.06000),
            new(1.23000, 0.06000),
            new(1.23000, 0.07000),
            new(1.22000, 0.07000),
            new(1.22000, 0.03000),
            new(1.23000, 0.03000),
            new(1.23000, 0.02000),
            new(1.30000, 0.02000),
            new(1.30000, 0.03000),
            new(1.29000, 0.03000),
            new(1.29000, 0.04000),
            new(1.32000, 0.04000),
            new(1.32000, 0.03000),
            new(1.31000, 0.03000),
            new(1.31000, 0.02000),
            new(1.35000, 0.02000),
            new(1.35000, 0.03000),
            new(1.36000, 0.03000),
            new(1.36000, 0.10200),
            new(1.35000, 0.10200),
            new(1.35000, 0.09200),
            new(1.34000, 0.09200),
            new(1.34000, 0.12200),
            new(1.35000, 0.12200),
            new(1.35000, 0.11200),
            new(1.36000, 0.11200),
            new(1.36000, 0.15200),
            new(1.35000, 0.15200),
            new(1.35000, 0.16200),
            new(1.29000, 0.16200),
            new(1.29000, 0.17200),
            new(1.30000, 0.17200),
            new(1.30000, 0.18200),
            new(1.29000, 0.18200),
            new(1.29000, 0.19200),
            new(1.32000, 0.19200),
            new(1.32000, 0.18200),
            new(1.31000, 0.18200),
            new(1.31000, 0.17200),
            new(1.35000, 0.17200),
            new(1.35000, 0.18200),
            new(1.36000, 0.18200),
            new(1.36000, 0.24300),
            new(1.37000, 0.24300),
            new(1.37000, 0.23200),
            new(1.38000, 0.23200),
            new(1.38000, 0.24200),
            new(1.39000, 0.24200),
            new(1.39000, 0.21200),
            new(1.38000, 0.21200),
            new(1.38000, 0.22200),
            new(1.37000, 0.22200),
            new(1.37000, 0.18200),
            new(1.38000, 0.18200),
            new(1.38000, 0.17200),
            new(1.44000, 0.17200),
            new(1.44000, 0.16200),
            new(1.43000, 0.16200),
            new(1.43000, 0.15200),
            new(1.44000, 0.15200),
            new(1.44000, 0.14200),
            new(1.41000, 0.14200),
            new(1.41000, 0.15200),
            new(1.42000, 0.15200),
            new(1.42000, 0.16200),
            new(1.38000, 0.16200),
            new(1.38000, 0.15200),
            new(1.37000, 0.15200),
            new(1.37000, 0.08000),
            new(1.38000, 0.08000),
            new(1.38000, 0.09000),
            new(1.39000, 0.09000),
            new(1.39000, 0.06000),
            new(1.38000, 0.06000),
            new(1.38000, 0.07000),
            new(1.37000, 0.07000),
            new(1.37000, 0.03000),
            new(1.38000, 0.03000),
            new(1.38000, 0.02000),
            new(1.45000, 0.02000),
            new(1.45000, 0.03000),
            new(1.44000, 0.03000),
            new(1.44000, 0.04000),
            new(1.47000, 0.04000),
            new(1.47000, 0.03000),
            new(1.46000, 0.03000),
            new(1.46000, 0.02000),
            new(1.50000, 0.02000),
            new(1.50000, 0.03000),
            new(1.51000, 0.03000),
            new(1.51000, 0.09100),
            new(1.52000, 0.09100),
            new(1.52000, 0.08000),
            new(1.53000, 0.08000),
            new(1.53000, 0.09000),
            new(1.54000, 0.09000),
            new(1.54000, 0.06000),
            new(1.53000, 0.06000),
            new(1.53000, 0.07000),
            new(1.52000, 0.07000),
            new(1.52000, 0.03000),
            new(1.53000, 0.03000),
            new(1.53000, 0.02000),
            new(1.60000, 0.02000),
            new(1.60000, 0.03000),
            new(1.59000, 0.03000),
            new(1.59000, 0.04000),
            new(1.62000, 0.04000),
            new(1.62000, 0.03000),
            new(1.61000, 0.03000),
            new(1.61000, 0.02000),
            new(1.65000, 0.02000),
            new(1.65000, 0.03000),
            new(1.66000, 0.03000),
            new(1.66000, 0.09100),
            new(1.67000, 0.09100),
            new(1.67000, 0.08000),
            new(1.68000, 0.08000),
            new(1.68000, 0.09000),
            new(1.69000, 0.09000),
            new(1.69000, 0.06000),
            new(1.68000, 0.06000),
            new(1.68000, 0.07000),
            new(1.67000, 0.07000),
            new(1.67000, 0.03000),
            new(1.68000, 0.03000),
            new(1.68000, 0.02000),
            new(1.75000, 0.02000),
            new(1.75000, 0.03000),
            new(1.74000, 0.03000),
            new(1.74000, 0.04000),
            new(1.77000, 0.04000),
            new(1.77000, 0.03000),
            new(1.76000, 0.03000),
            new(1.76000, 0.02000),
            new(1.80000, 0.02000),
            new(1.80000, 0.03000),
            new(1.81000, 0.03000),
            new(1.81000, 0.10200),
            new(1.80000, 0.10200),
            new(1.80000, 0.09200),
            new(1.79000, 0.09200),
            new(1.79000, 0.12200),
            new(1.80000, 0.12200),
            new(1.80000, 0.11200),
            new(1.81000, 0.11200),
            new(1.81000, 0.15200),
            new(1.80000, 0.15200),
            new(1.80000, 0.16200),
            new(1.74000, 0.16200),
            new(1.74000, 0.17200),
            new(1.75000, 0.17200),
            new(1.75000, 0.18200),
            new(1.74000, 0.18200),
            new(1.74000, 0.19200),
            new(1.77000, 0.19200),
            new(1.77000, 0.18200),
            new(1.76000, 0.18200),
            new(1.76000, 0.17200),
            new(1.80000, 0.17200),
            new(1.80000, 0.18200),
            new(1.81000, 0.18200),
            new(1.81000, 0.24300),
            new(1.82000, 0.24300),
            new(1.82000, 0.23200),
            new(1.83000, 0.23200),
            new(1.83000, 0.24200),
            new(1.84000, 0.24200),
            new(1.84000, 0.21200),
            new(1.83000, 0.21200),
            new(1.83000, 0.22200),
            new(1.82000, 0.22200),
            new(1.82000, 0.18200),
            new(1.83000, 0.18200),
            new(1.83000, 0.17200),
            new(1.89000, 0.17200),
            new(1.89000, 0.16200),
            new(1.88000, 0.16200),
            new(1.88000, 0.15200),
            new(1.89000, 0.15200),
            new(1.89000, 0.14200),
            new(1.86000, 0.14200),
            new(1.86000, 0.15200),
            new(1.87000, 0.15200),
            new(1.87000, 0.16200),
            new(1.83000, 0.16200),
            new(1.83000, 0.15200),
            new(1.82000, 0.15200),
            new(1.82000, 0.08000),
            new(1.83000, 0.08000),
            new(1.83000, 0.09000),
            new(1.84000, 0.09000),
            new(1.84000, 0.06000),
            new(1.83000, 0.06000),
            new(1.83000, 0.07000),
            new(1.82000, 0.07000),
            new(1.82000, 0.03000),
            new(1.83000, 0.03000),
            new(1.83000, 0.02000),
            new(1.90000, 0.02000),
            new(1.90000, 0.03000),
            new(1.89000, 0.03000),
            new(1.89000, 0.04000),
            new(1.92000, 0.04000),
            new(1.92000, 0.03000),
            new(1.91000, 0.03000),
            new(1.91000, 0.02000),
            new(1.95000, 0.02000),
            new(1.95000, 0.03000),
            new(1.96000, 0.03000),
            new(1.96000, 0.09100),
            new(1.97000, 0.09100),
            new(1.97000, 0.08000),
            new(1.98000, 0.08000),
            new(1.98000, 0.09000),
            new(1.99000, 0.09000),
            new(1.99000, 0.06000),
            new(1.98000, 0.06000),
            new(1.98000, 0.07000),
            new(1.97000, 0.07000),
            new(1.97000, 0.03000),
            new(1.98000, 0.03000),
            new(1.98000, 0.02000),
            new(2.04000, 0.02000),
            new(2.04000, 0.01000),
            new(2.03000, 0.01000),
            new(2.03000, 0.00000),
            new(2.04000, 0.00000),
            new(2.04000, -0.01000),
            new(2.01000, -0.01000),
            new(2.01000, 0.00000),
            new(2.02000, 0.00000),
            new(2.02000, 0.01000),
            new(1.98000, 0.01000),
            new(1.98000, 0.00000),
            new(1.97000, 0.00000),
            new(1.97000, -0.06100),
            new(1.96000, -0.06100),
            new(1.96000, -0.05000),
            new(1.95000, -0.05000),
            new(1.95000, -0.06000),
            new(1.94000, -0.06000),
            new(1.94000, -0.03000),
            new(1.95000, -0.03000),
            new(1.95000, -0.04000),
            new(1.96000, -0.04000),
            new(1.96000, 0.00000),
            new(1.95000, 0.00000),
            new(1.95000, 0.01000),
            new(1.88000, 0.01000),
            new(1.88000, 0.00000),
            new(1.89000, 0.00000),
            new(1.89000, -0.01000),
            new(1.86000, -0.01000),
            new(1.86000, 0.00000),
            new(1.87000, 0.00000),
            new(1.87000, 0.01000),
            new(1.83000, 0.01000),
            new(1.83000, 0.00000),
            new(1.82000, 0.00000),
            new(1.82000, -0.07200),
            new(1.83000, -0.07200),
            new(1.83000, -0.06200),
            new(1.84000, -0.06200),
            new(1.84000, -0.09200),
            new(1.83000, -0.09200),
            new(1.83000, -0.08200),
            new(1.82000, -0.08200),
            new(1.82000, -0.12200),
            new(1.83000, -0.12200),
            new(1.83000, -0.13200),
            new(1.89000, -0.13200),
            new(1.89000, -0.14200),
            new(1.88000, -0.14200),
            new(1.88000, -0.15200),
            new(1.89000, -0.15200),
            new(1.89000, -0.16200),
            new(1.86000, -0.16200),
            new(1.86000, -0.15200),
            new(1.87000, -0.15200),
            new(1.87000, -0.14200),
            new(1.83000, -0.14200),
            new(1.83000, -0.15200),
            new(1.82000, -0.15200),
            new(1.82000, -0.21300),
            new(1.81000, -0.21300),
            new(1.81000, -0.20200),
            new(1.80000, -0.20200),
            new(1.80000, -0.21200),
            new(1.79000, -0.21200),
            new(1.79000, -0.18200),
            new(1.80000, -0.18200),
            new(1.80000, -0.19200),
            new(1.81000, -0.19200),
            new(1.81000, -0.15200),
            new(1.80000, -0.15200),
            new(1.80000, -0.14200),
            new(1.74000, -0.14200),
            new(1.74000, -0.13200),
            new(1.75000, -0.13200),
            new(1.75000, -0.12200),
            new(1.74000, -0.12200),
            new(1.74000, -0.11200),
            new(1.77000, -0.11200),
            new(1.77000, -0.12200),
            new(1.76000, -0.12200),
            new(1.76000, -0.13200),
            new(1.80000, -0.13200),
            new(1.80000, -0.12200),
            new(1.81000, -0.12200),
            new(1.81000, -0.05000),
            new(1.80000, -0.05000),
            new(1.80000, -0.06000),
            new(1.79000, -0.06000),
            new(1.79000, -0.03000),
            new(1.80000, -0.03000),
            new(1.80000, -0.04000),
            new(1.81000, -0.04000),
            new(1.81000, 0.00000),
            new(1.80000, 0.00000),
            new(1.80000, 0.01000),
            new(1.73000, 0.01000),
            new(1.73000, 0.00000),
            new(1.74000, 0.00000),
            new(1.74000, -0.01000),
            new(1.71000, -0.01000),
            new(1.71000, 0.00000),
            new(1.72000, 0.00000),
            new(1.72000, 0.01000),
            new(1.68000, 0.01000),
            new(1.68000, 0.00000),
            new(1.67000, 0.00000),
            new(1.67000, -0.06100),
            new(1.66000, -0.06100),
            new(1.66000, -0.05000),
            new(1.65000, -0.05000),
            new(1.65000, -0.06000),
            new(1.64000, -0.06000),
            new(1.64000, -0.03000),
            new(1.65000, -0.03000),
            new(1.65000, -0.04000),
            new(1.66000, -0.04000),
            new(1.66000, 0.00000),
            new(1.65000, 0.00000),
            new(1.65000, 0.01000),
            new(1.58000, 0.01000),
            new(1.58000, 0.00000),
            new(1.59000, 0.00000),
            new(1.59000, -0.01000),
            new(1.56000, -0.01000),
            new(1.56000, 0.00000),
            new(1.57000, 0.00000),
            new(1.57000, 0.01000),
            new(1.53000, 0.01000),
            new(1.53000, 0.00000),
            new(1.52000, 0.00000),
            new(1.52000, -0.06100),
            new(1.51000, -0.06100),
            new(1.51000, -0.05000),
            new(1.50000, -0.05000),
            new(1.50000, -0.06000),
            new(1.49000, -0.06000),
            new(1.49000, -0.03000),
            new(1.50000, -0.03000),
            new(1.50000, -0.04000),
            new(1.51000, -0.04000),
            new(1.51000, 0.00000),
            new(1.50000, 0.00000),
            new(1.50000, 0.01000),
            new(1.43000, 0.01000),
            new(1.43000, 0.00000),
            new(1.44000, 0.00000),
            new(1.44000, -0.01000),
            new(1.41000, -0.01000),
            new(1.41000, 0.00000),
            new(1.42000, 0.00000),
            new(1.42000, 0.01000),
            new(1.38000, 0.01000),
            new(1.38000, 0.00000),
            new(1.37000, 0.00000),
            new(1.37000, -0.07200),
            new(1.38000, -0.07200),
            new(1.38000, -0.06200),
            new(1.39000, -0.06200),
            new(1.39000, -0.09200),
            new(1.38000, -0.09200),
            new(1.38000, -0.08200),
            new(1.37000, -0.08200),
            new(1.37000, -0.12200),
            new(1.38000, -0.12200),
            new(1.38000, -0.13200),
            new(1.44000, -0.13200),
            new(1.44000, -0.14200),
            new(1.43000, -0.14200),
            new(1.43000, -0.15200),
            new(1.44000, -0.15200),
            new(1.44000, -0.16200),
            new(1.41000, -0.16200),
            new(1.41000, -0.15200),
            new(1.42000, -0.15200),
            new(1.42000, -0.14200),
            new(1.38000, -0.14200),
            new(1.38000, -0.15200),
            new(1.37000, -0.15200),
            new(1.37000, -0.21300),
            new(1.36000, -0.21300),
            new(1.36000, -0.20200),
            new(1.35000, -0.20200),
            new(1.35000, -0.21200),
            new(1.34000, -0.21200),
            new(1.34000, -0.18200),
            new(1.35000, -0.18200),
            new(1.35000, -0.19200),
            new(1.36000, -0.19200),
            new(1.36000, -0.15200),
            new(1.35000, -0.15200),
            new(1.35000, -0.14200),
            new(1.29000, -0.14200),
            new(1.29000, -0.13200),
            new(1.30000, -0.13200),
            new(1.30000, -0.12200),
            new(1.29000, -0.12200),
            new(1.29000, -0.11200),
            new(1.32000, -0.11200),
            new(1.32000, -0.12200),
            new(1.31000, -0.12200),
            new(1.31000, -0.13200),
            new(1.35000, -0.13200),
            new(1.35000, -0.12200),
            new(1.36000, -0.12200),
            new(1.36000, -0.05000),
            new(1.35000, -0.05000),
            new(1.35000, -0.06000),
            new(1.34000, -0.06000),
            new(1.34000, -0.03000),
            new(1.35000, -0.03000),
            new(1.35000, -0.04000),
            new(1.36000, -0.04000),
            new(1.36000, 0.00000),
            new(1.35000, 0.00000),
            new(1.35000, 0.01000),
            new(1.28000, 0.01000),
            new(1.28000, 0.00000),
            new(1.29000, 0.00000),
            new(1.29000, -0.01000),
            new(1.26000, -0.01000),
            new(1.26000, 0.00000),
            new(1.27000, 0.00000),
            new(1.27000, 0.01000),
            new(1.23000, 0.01000),
            new(1.23000, 0.00000),
            new(1.22000, 0.00000),
            new(1.22000, -0.06100),
            new(1.21000, -0.06100),
            new(1.21000, -0.05000),
            new(1.20000, -0.05000),
            new(1.20000, -0.06000),
            new(1.19000, -0.06000),
            new(1.19000, -0.03000),
            new(1.20000, -0.03000),
            new(1.20000, -0.04000),
            new(1.21000, -0.04000),
            new(1.21000, 0.00000),
            new(1.20000, 0.00000),
            new(1.20000, 0.01000),
            new(1.13000, 0.01000),
            new(1.13000, 0.00000),
            new(1.14000, 0.00000),
            new(1.14000, -0.01000),
            new(1.11000, -0.01000),
            new(1.11000, 0.00000),
            new(1.12000, 0.00000),
            new(1.12000, 0.01000),
            new(1.08000, 0.01000),
            new(1.08000, 0.00000),
            new(1.07000, 0.00000),
            new(1.07000, -0.06100),
            new(1.06000, -0.06100),
            new(1.06000, -0.05000),
            new(1.05000, -0.05000),
            new(1.05000, -0.06000),
            new(1.04000, -0.06000),
            new(1.04000, -0.03000),
            new(1.05000, -0.03000),
            new(1.05000, -0.04000),
            new(1.06000, -0.04000),
            new(1.06000, 0.00000),
            new(1.05000, 0.00000),
            new(1.05000, 0.01000),
            new(0.98000, 0.01000),
            new(0.98000, 0.00000),
            new(0.99000, 0.00000),
            new(0.99000, -0.01000),
            new(0.96000, -0.01000),
            new(0.96000, 0.00000),
            new(0.97000, 0.00000),
            new(0.97000, 0.01000),
            new(0.93000, 0.01000),
            new(0.93000, 0.00000),
            new(0.92000, 0.00000),
            new(0.92000, -0.07200),
            new(0.93000, -0.07200),
            new(0.93000, -0.06200),
            new(0.94000, -0.06200),
            new(0.94000, -0.09200),
            new(0.93000, -0.09200),
            new(0.93000, -0.08200),
            new(0.92000, -0.08200),
            new(0.92000, -0.12200),
            new(0.93000, -0.12200),
            new(0.93000, -0.13200),
            new(0.99000, -0.13200),
            new(0.99000, -0.14200),
            new(0.98000, -0.14200),
            new(0.98000, -0.15200),
            new(0.99000, -0.15200),
            new(0.99000, -0.16200),
            new(0.96000, -0.16200),
            new(0.96000, -0.15200),
            new(0.97000, -0.15200),
            new(0.97000, -0.14200),
            new(0.93000, -0.14200),
            new(0.93000, -0.15200),
            new(0.92000, -0.15200),
            new(0.92000, -0.21300),
            new(0.91000, -0.21300),
            new(0.91000, -0.20200),
            new(0.90000, -0.20200),
            new(0.90000, -0.21200),
            new(0.89000, -0.21200),
            new(0.89000, -0.18200),
            new(0.90000, -0.18200),
            new(0.90000, -0.19200),
            new(0.91000, -0.19200),
            new(0.91000, -0.15200),
            new(0.90000, -0.15200),
            new(0.90000, -0.14200),
            new(0.84000, -0.14200),
            new(0.84000, -0.13200),
            new(0.85000, -0.13200),
            new(0.85000, -0.12200),
            new(0.84000, -0.12200),
            new(0.84000, -0.11200),
            new(0.87000, -0.11200),
            new(0.87000, -0.12200),
            new(0.86000, -0.12200),
            new(0.86000, -0.13200),
            new(0.90000, -0.13200),
            new(0.90000, -0.12200),
            new(0.91000, -0.12200),
            new(0.91000, -0.05000),
            new(0.90000, -0.05000),
            new(0.90000, -0.06000),
            new(0.89000, -0.06000),
            new(0.89000, -0.03000),
            new(0.90000, -0.03000),
            new(0.90000, -0.04000),
            new(0.91000, -0.04000),
            new(0.91000, 0.00000),
            new(0.90000, 0.00000),
            new(0.90000, 0.01000),
            new(0.83000, 0.01000),
            new(0.83000, 0.00000),
            new(0.84000, 0.00000),
            new(0.84000, -0.01000),
            new(0.81000, -0.01000),
            new(0.81000, 0.00000),
            new(0.82000, 0.00000),
            new(0.82000, 0.01000),
            new(0.78000, 0.01000),
            new(0.78000, 0.00000),
            new(0.77000, 0.00000),
            new(0.77000, -0.06100),
            new(0.76000, -0.06100),
            new(0.76000, -0.05000),
            new(0.75000, -0.05000),
            new(0.75000, -0.06000),
            new(0.74000, -0.06000),
            new(0.74000, -0.03000),
            new(0.75000, -0.03000),
            new(0.75000, -0.04000),
            new(0.76000, -0.04000),
            new(0.76000, 0.00000),
            new(0.75000, 0.00000),
            new(0.75000, 0.01000),
            new(0.68000, 0.01000),
            new(0.68000, 0.00000),
            new(0.69000, 0.00000),
            new(0.69000, -0.01000),
            new(0.66000, -0.01000),
            new(0.66000, 0.00000),
            new(0.67000, 0.00000),
            new(0.67000, 0.01000),
            new(0.63000, 0.01000),
            new(0.63000, 0.00000),
            new(0.62000, 0.00000),
            new(0.62000, -0.06100),
            new(0.61000, -0.06100),
            new(0.61000, -0.05000),
            new(0.60000, -0.05000),
            new(0.60000, -0.06000),
            new(0.59000, -0.06000),
            new(0.59000, -0.03000),
            new(0.60000, -0.03000),
            new(0.60000, -0.04000),
            new(0.61000, -0.04000),
            new(0.61000, 0.00000),
            new(0.60000, 0.00000),
            new(0.60000, 0.01000),
            new(0.53000, 0.01000),
            new(0.53000, 0.00000),
            new(0.54000, 0.00000),
            new(0.54000, -0.01000),
            new(0.51000, -0.01000),
            new(0.51000, 0.00000),
            new(0.52000, 0.00000),
            new(0.52000, 0.01000),
            new(0.48000, 0.01000),
            new(0.48000, 0.00000),
            new(0.47000, 0.00000),
            new(0.47000, -0.07200),
            new(0.48000, -0.07200),
            new(0.48000, -0.06200),
            new(0.49000, -0.06200),
            new(0.49000, -0.09200),
            new(0.48000, -0.09200),
            new(0.48000, -0.08200),
            new(0.47000, -0.08200),
            new(0.47000, -0.12200),
            new(0.48000, -0.12200),
            new(0.48000, -0.13200),
            new(0.54000, -0.13200),
            new(0.54000, -0.14200),
            new(0.53000, -0.14200),
            new(0.53000, -0.15200),
            new(0.54000, -0.15200),
            new(0.54000, -0.16200),
            new(0.51000, -0.16200),
            new(0.51000, -0.15200),
            new(0.52000, -0.15200),
            new(0.52000, -0.14200),
            new(0.48000, -0.14200),
            new(0.48000, -0.15200),
            new(0.47000, -0.15200),
            new(0.47000, -0.21300),
            new(0.46000, -0.21300),
            new(0.46000, -0.20200),
            new(0.45000, -0.20200),
            new(0.45000, -0.21200),
            new(0.44000, -0.21200),
            new(0.44000, -0.18200),
            new(0.45000, -0.18200),
            new(0.45000, -0.19200),
            new(0.46000, -0.19200),
            new(0.46000, -0.15200),
            new(0.45000, -0.15200),
            new(0.45000, -0.14200),
            new(0.39000, -0.14200),
            new(0.39000, -0.13200),
            new(0.40000, -0.13200),
            new(0.40000, -0.12200),
            new(0.39000, -0.12200),
            new(0.39000, -0.11200),
            new(0.42000, -0.11200),
            new(0.42000, -0.12200),
            new(0.41000, -0.12200),
            new(0.41000, -0.13200),
            new(0.45000, -0.13200),
            new(0.45000, -0.12200),
            new(0.46000, -0.12200),
            new(0.46000, -0.05000),
            new(0.45000, -0.05000),
            new(0.45000, -0.06000),
            new(0.44000, -0.06000),
            new(0.44000, -0.03000),
            new(0.45000, -0.03000),
            new(0.45000, -0.04000),
            new(0.46000, -0.04000),
            new(0.46000, 0.00000),
            new(0.45000, 0.00000),
            new(0.45000, 0.01000),
            new(0.38000, 0.01000),
            new(0.38000, 0.00000),
            new(0.39000, 0.00000),
            new(0.39000, -0.01000),
            new(0.36000, -0.01000),
            new(0.36000, 0.00000),
            new(0.37000, 0.00000),
            new(0.37000, 0.01000),
            new(0.33000, 0.01000),
            new(0.33000, 0.00000),
            new(0.32000, 0.00000),
            new(0.32000, -0.06100),
            new(0.31000, -0.06100),
            new(0.31000, -0.05000),
            new(0.30000, -0.05000),
            new(0.30000, -0.06000),
            new(0.29000, -0.06000),
            new(0.29000, -0.03000),
            new(0.30000, -0.03000),
            new(0.30000, -0.04000),
            new(0.31000, -0.04000),
            new(0.31000, 0.00000),
            new(0.30000, 0.00000),
            new(0.30000, 0.01000),
            new(0.23000, 0.01000),
            new(0.23000, 0.00000),
            new(0.24000, 0.00000),
            new(0.24000, -0.01000),
            new(0.21000, -0.01000),
            new(0.21000, 0.00000),
            new(0.22000, 0.00000),
            new(0.22000, 0.01000),
            new(0.18000, 0.01000),
            new(0.18000, 0.00000),
            new(0.17000, 0.00000),
            new(0.17000, -0.06100),
            new(0.16000, -0.06100),
            new(0.16000, -0.05000),
            new(0.15000, -0.05000),
            new(0.15000, -0.06000),
            new(0.14000, -0.06000),
            new(0.14000, -0.03000),
            new(0.15000, -0.03000),
            new(0.15000, -0.04000),
            new(0.16000, -0.04000),
            new(0.16000, 0.00000),
            new(0.15000, 0.00000),
            new(0.15000, 0.01000),
            new(0.08000, 0.01000),
            new(0.08000, 0.00000),
            new(0.09000, 0.00000),
            new(0.09000, -0.01000),
            new(0.06000, -0.01000),
            new(0.06000, 0.00000),
            new(0.07000, 0.00000),
            new(0.07000, 0.01000),
            new(0.03000, 0.01000),
            new(0.03000, 0.00000),
            new(0.02000, 0.00000),
            new(0.02000, -0.07200),
            new(0.03000, -0.07200),
            new(0.03000, -0.06200),
            new(0.04000, -0.06200),
            new(0.04000, -0.09200),
            new(0.03000, -0.09200),
            new(0.03000, -0.08200),
            new(0.02000, -0.08200),
            new(0.02000, -0.12200),
            new(0.03000, -0.12200),
            new(0.03000, -0.13200),
            new(0.09000, -0.13200),
            new(0.09000, -0.14200),
            new(0.08000, -0.14200),
            new(0.08000, -0.15200),
            new(0.09000, -0.15200),
            new(0.09000, -0.16200),
            new(0.06000, -0.16200),
            new(0.06000, -0.15200),
            new(0.07000, -0.15200),
            new(0.07000, -0.14200),
            new(0.03000, -0.14200),
            new(0.03000, -0.15200),
            new(0.02000, -0.15200),
            new(0.02000, -0.21300),
            new(0.01000, -0.21300),
        };

        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Geo clean-up....");
        sw.Restart();
        poly = GeoWrangler.clockwiseAndReorderXY(poly);
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Conversion....");
        sw.Restart();
        Paths64 done = new() { GeoWrangler.path64FromPathD(poly) };
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        bool abort = false;

        Console.WriteLine("  Decomposition (vertical)....");
        sw.Restart();
        Paths64 ns = GeoWrangler.rectangular_decomposition(ref abort, done, maxRayLength: rayLength);

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
        Path64 points_1 = TestGeometry.ortho_fractal_1();
        points_1 = GeoWrangler.close(points_1);
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        partFour_do(points_1, "complex_loop");

        Console.WriteLine(" Part 2....");
        Console.WriteLine("  Preparing....");
        sw.Start();
        Path64 points_2 = TestGeometry.ortho_fractal_2();
        points_2 = GeoWrangler.close(points_2);
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");
        partFour_do(points_2, "complex_loop_rot");
    }

    private static void partFour_do(Path64 points, string baseString)
    {
        System.Diagnostics.Stopwatch sw = new();

        bool vertical = true;
        bool abort = false;

        Console.WriteLine("  Keyhole....");
        // Give the keyholder a whirl:
        sw.Restart();
        Path64 toDecomp = GeoWrangler.makeKeyHole(points, reverseEval:false, biDirectionalEval:true)[0];
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Query....");
        sw.Restart();
        Path64 bounds = GeoWrangler.getBounds(toDecomp);
        PointD dist = GeoWrangler.distanceBetweenPoints_point(bounds[0], bounds[1]);
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Decomposition (vertical)....");
        sw.Restart();
        Paths64 decompOut = GeoWrangler.rectangular_decomposition(ref abort, toDecomp,
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

        Paths64 out_decomp = new();
        for (int i = 0; i < polydata.Count; i++)
        {
         PathD points = new (polydata[i]);
         points = GeoWrangler.removeDuplicates(points);
         points = GeoWrangler.stripColinear(points);
         points = GeoWrangler.clockwiseAndReorderXY(points);
         Path64 toKeyHoler = GeoWrangler.path64FromPathD(points);
         
         Path64 toDecomp = GeoWrangler.makeKeyHole(GeoWrangler.sliverGapRemoval(toKeyHoler), reverseEval:false, biDirectionalEval:false)[0];
         Path64  bounds = GeoWrangler.getBounds(toDecomp);
         PointD dist = GeoWrangler.distanceBetweenPoints_point(bounds[0], bounds[1]);

         Paths64 decompOut = GeoWrangler.rectangular_decomposition(ref abort, toDecomp,
          maxRayLength: (long) Math.Max(Math.Abs(dist.x), Math.Abs(dist.y)), vertical: vertical);
         
         out_decomp.AddRange(decompOut);
        }
    }


    private static void writeToLayout(string filename, Path64 orig, Paths64 decomped)
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

        gcell.addPolygon(orig, 1, 0);

        for (int i = 0; i < decomped.Count; i++)
        {
            gcell.addBox(decomped[i], 1, 1);
        }

        g.setDrawing(drawing_);
        g.setValid(true);

        gds.gdsWriter gw = new(g, "../../../../../decomp_out/" + filename + "_partitiontest.gds");
        gw.save();

        oasis.oasWriter ow = new(g, "../../../../../decomp_out/" + filename + "_partitiontest.oas");
        ow.save();
    }
}