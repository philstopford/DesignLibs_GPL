using Clipper2Lib;
using geoWrangler;
using geoCoreLib;
using PartitionTestGeometrySource;
using NUnit.Framework;
using utility;

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

        Console.WriteLine("Part Four (takes less time than part three)");
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
        writeToLayout("l", L, l);
        Assert.AreEqual(l.Count, 2);
        Assert.AreEqual(l[0][0].x, 0);
        Assert.AreEqual(l[0][0].y, 0);
        Assert.AreEqual(l[0][1].x, 10);
        Assert.AreEqual(l[0][1].y, 0);
        Assert.AreEqual(l[0][2].x, 10);
        Assert.AreEqual(l[0][2].y, 50);
        Assert.AreEqual(l[0][3].x, 0);
        Assert.AreEqual(l[0][3].y, 50);
        Assert.AreEqual(l[0][4].x, 0);
        Assert.AreEqual(l[0][4].y, 0);
        Assert.AreEqual(l[1][0].x, 10);
        Assert.AreEqual(l[1][0].y, 0);
        Assert.AreEqual(l[1][1].x, 60);
        Assert.AreEqual(l[1][1].y, 0);
        Assert.AreEqual(l[1][2].x, 60);
        Assert.AreEqual(l[1][2].y, 20);
        Assert.AreEqual(l[1][3].x, 10);
        Assert.AreEqual(l[1][3].y, 20);
        Assert.AreEqual(l[1][4].x, 10);
        Assert.AreEqual(l[1][4].y, 0);
        
        PathsD lccw = GeoWrangler.rectangular_decomposition(ref abort, L_ccw, maxRayLength: rayLength);
        writeToLayout("lccw", L_ccw, lccw);
        Assert.AreEqual(lccw.Count, 2);
        Assert.AreEqual(lccw[0][0].x, 0);
        Assert.AreEqual(lccw[0][0].y, 0);
        Assert.AreEqual(lccw[0][1].x, 10);
        Assert.AreEqual(lccw[0][1].y, 0);
        Assert.AreEqual(lccw[0][2].x, 10);
        Assert.AreEqual(lccw[0][2].y, 50);
        Assert.AreEqual(lccw[0][3].x, 0);
        Assert.AreEqual(lccw[0][3].y, 50);
        Assert.AreEqual(lccw[0][4].x, 0);
        Assert.AreEqual(lccw[0][4].y, 0);
        Assert.AreEqual(lccw[1][0].x, 10);
        Assert.AreEqual(lccw[1][0].y, 0);
        Assert.AreEqual(lccw[1][1].x, 60);
        Assert.AreEqual(lccw[1][1].y, 0);
        Assert.AreEqual(lccw[1][2].x, 60);
        Assert.AreEqual(lccw[1][2].y, 20);
        Assert.AreEqual(lccw[1][3].x, 10);
        Assert.AreEqual(lccw[1][3].y, 20);
        Assert.AreEqual(lccw[1][4].x, 10);
        Assert.AreEqual(lccw[1][4].y, 0);
        
        PathsD rl = GeoWrangler.rectangular_decomposition(ref abort, rL, maxRayLength: rayLength);
        writeToLayout("rl", rL, rl);
        Assert.AreEqual(rl.Count, 2);
        Assert.AreEqual(rl[0][0].x, 10);
        Assert.AreEqual(rl[0][0].y, 0);
        Assert.AreEqual(rl[0][1].x, 60);
        Assert.AreEqual(rl[0][1].y, 0);
        Assert.AreEqual(rl[0][2].x, 60);
        Assert.AreEqual(rl[0][2].y, 50);
        Assert.AreEqual(rl[0][3].x, 10);
        Assert.AreEqual(rl[0][3].y, 50);
        Assert.AreEqual(rl[0][4].x, 10);
        Assert.AreEqual(rl[0][4].y, 0);
        Assert.AreEqual(rl[1][0].x, 0);
        Assert.AreEqual(rl[1][0].y, 0);
        Assert.AreEqual(rl[1][1].x, 10);
        Assert.AreEqual(rl[1][1].y, 0);
        Assert.AreEqual(rl[1][2].x, 10);
        Assert.AreEqual(rl[1][2].y, 20);
        Assert.AreEqual(rl[1][3].x, 0);
        Assert.AreEqual(rl[1][3].y, 20);
        Assert.AreEqual(rl[1][4].x, 0);
        Assert.AreEqual(rl[1][4].y, 0);


        PathsD u = GeoWrangler.rectangular_decomposition(ref abort, U, maxRayLength: rayLength);
        writeToLayout("u", U, u);
        Assert.AreEqual(u.Count, 3);
        Assert.AreEqual(u[0][0].x, 0);
        Assert.AreEqual(u[0][0].y, 0);
        Assert.AreEqual(u[0][1].x, 10);
        Assert.AreEqual(u[0][1].y, 0);
        Assert.AreEqual(u[0][2].x, 10);
        Assert.AreEqual(u[0][2].y, 50);
        Assert.AreEqual(u[0][3].x, 0);
        Assert.AreEqual(u[0][3].y, 50);
        Assert.AreEqual(u[0][4].x, 0);
        Assert.AreEqual(u[0][4].y, 0);
        Assert.AreEqual(u[1][0].x, 60);
        Assert.AreEqual(u[1][0].y, 0);
        Assert.AreEqual(u[1][1].x, 120);
        Assert.AreEqual(u[1][1].y, 0);
        Assert.AreEqual(u[1][2].x, 120);
        Assert.AreEqual(u[1][2].y, 80);
        Assert.AreEqual(u[1][3].x, 60);
        Assert.AreEqual(u[1][3].y, 80);
        Assert.AreEqual(u[1][4].x, 60);
        Assert.AreEqual(u[1][4].y, 0);
        Assert.AreEqual(u[2][0].x, 10);
        Assert.AreEqual(u[2][0].y, 0);
        Assert.AreEqual(u[2][1].x, 60);
        Assert.AreEqual(u[2][1].y, 0);
        Assert.AreEqual(u[2][2].x, 60);
        Assert.AreEqual(u[2][2].y, 20);
        Assert.AreEqual(u[2][3].x, 10);
        Assert.AreEqual(u[2][3].y, 20);
        Assert.AreEqual(u[2][4].x, 10);
        Assert.AreEqual(u[2][4].y, 0);
        
        PathsD t = GeoWrangler.rectangular_decomposition(ref abort, T, maxRayLength: rayLength);
        writeToLayout("t", T, t);
        Assert.AreEqual(t.Count, 3);
        Assert.AreEqual(t[0][0].x, 60);
        Assert.AreEqual(t[0][0].y, 50);
        Assert.AreEqual(t[0][1].x, 80);
        Assert.AreEqual(t[0][1].y, 50);
        Assert.AreEqual(t[0][2].x, 80);
        Assert.AreEqual(t[0][2].y, 80);
        Assert.AreEqual(t[0][3].x, 60);
        Assert.AreEqual(t[0][3].y, 80);
        Assert.AreEqual(t[0][4].x, 60);
        Assert.AreEqual(t[0][4].y, 50);
        Assert.AreEqual(t[1][0].x, 40);
        Assert.AreEqual(t[1][0].y, 0);
        Assert.AreEqual(t[1][1].x, 60);
        Assert.AreEqual(t[1][1].y, 0);
        Assert.AreEqual(t[1][2].x, 60);
        Assert.AreEqual(t[1][2].y, 80);
        Assert.AreEqual(t[1][3].x, 40);
        Assert.AreEqual(t[1][3].y, 80);
        Assert.AreEqual(t[1][4].x, 40);
        Assert.AreEqual(t[1][4].y, 0);
        Assert.AreEqual(t[2][0].x, 0);
        Assert.AreEqual(t[2][0].y, 50);
        Assert.AreEqual(t[2][1].x, 40);
        Assert.AreEqual(t[2][1].y, 50);
        Assert.AreEqual(t[2][2].x, 40);
        Assert.AreEqual(t[2][2].y, 80);
        Assert.AreEqual(t[2][3].x, 0);
        Assert.AreEqual(t[2][3].y, 80);
        Assert.AreEqual(t[2][4].x, 0);
        Assert.AreEqual(t[2][4].y, 50);
        
        PathsD x = GeoWrangler.rectangular_decomposition(ref abort, X, maxRayLength: rayLength);
        writeToLayout("x", X, x);
        Assert.AreEqual(x.Count, 3);
        Assert.AreEqual(x[0][0].x, 0);
        Assert.AreEqual(x[0][0].y, 50);
        Assert.AreEqual(x[0][1].x, 60);
        Assert.AreEqual(x[0][1].y, 50);
        Assert.AreEqual(x[0][2].x, 60);
        Assert.AreEqual(x[0][2].y, 80);
        Assert.AreEqual(x[0][3].x, 0);
        Assert.AreEqual(x[0][3].y, 80);
        Assert.AreEqual(x[0][4].x, 0);
        Assert.AreEqual(x[0][4].y, 50);
        Assert.AreEqual(x[1][0].x, 60);
        Assert.AreEqual(x[1][0].y, 20);
        Assert.AreEqual(x[1][1].x, 80);
        Assert.AreEqual(x[1][1].y, 20);
        Assert.AreEqual(x[1][2].x, 80);
        Assert.AreEqual(x[1][2].y, 100);
        Assert.AreEqual(x[1][3].x, 60);
        Assert.AreEqual(x[1][3].y, 100);
        Assert.AreEqual(x[1][4].x, 60);
        Assert.AreEqual(x[1][4].y, 20);
        Assert.AreEqual(x[2][0].x, 80);
        Assert.AreEqual(x[2][0].y, 50);
        Assert.AreEqual(x[2][1].x, 100);
        Assert.AreEqual(x[2][1].y, 50);
        Assert.AreEqual(x[2][2].x, 100);
        Assert.AreEqual(x[2][2].y, 80);
        Assert.AreEqual(x[2][3].x, 80);
        Assert.AreEqual(x[2][3].y, 80);
        Assert.AreEqual(x[2][4].x, 80);
        Assert.AreEqual(x[2][4].y, 50);
        
        PathsD s = GeoWrangler.rectangular_decomposition(ref abort, S, maxRayLength: rayLength);
        writeToLayout("s", S, s);
        Assert.AreEqual(s.Count, 5);
        Assert.AreEqual(s[0][0].x, 0);
        Assert.AreEqual(s[0][0].y, 50);
        Assert.AreEqual(s[0][1].x, 20);
        Assert.AreEqual(s[0][1].y, 50);
        Assert.AreEqual(s[0][2].x, 20);
        Assert.AreEqual(s[0][2].y, 110);
        Assert.AreEqual(s[0][3].x, 0);
        Assert.AreEqual(s[0][3].y, 110);
        Assert.AreEqual(s[0][4].x, 0);
        Assert.AreEqual(s[0][4].y, 50);
        Assert.AreEqual(s[1][0].x, 0);
        Assert.AreEqual(s[1][0].y, 0);
        Assert.AreEqual(s[1][1].x, 20);
        Assert.AreEqual(s[1][1].y, 0);
        Assert.AreEqual(s[1][2].x, 20);
        Assert.AreEqual(s[1][2].y, 20);
        Assert.AreEqual(s[1][3].x, 0);
        Assert.AreEqual(s[1][3].y, 20);
        Assert.AreEqual(s[1][4].x, 0);
        Assert.AreEqual(s[1][4].y, 0);
        Assert.AreEqual(s[2][0].x, 80);
        Assert.AreEqual(s[2][0].y, 0);
        Assert.AreEqual(s[2][1].x, 100);
        Assert.AreEqual(s[2][1].y, 0);
        Assert.AreEqual(s[2][2].x, 100);
        Assert.AreEqual(s[2][2].y, 60);
        Assert.AreEqual(s[2][3].x, 80);
        Assert.AreEqual(s[2][3].y, 60);
        Assert.AreEqual(s[2][4].x, 80);
        Assert.AreEqual(s[2][4].y, 0);
        Assert.AreEqual(s[3][0].x, 80);
        Assert.AreEqual(s[3][0].y, 80);
        Assert.AreEqual(s[3][1].x, 100);
        Assert.AreEqual(s[3][1].y, 80);
        Assert.AreEqual(s[3][2].x, 100);
        Assert.AreEqual(s[3][2].y, 110);
        Assert.AreEqual(s[3][3].x, 80);
        Assert.AreEqual(s[3][3].y, 110);
        Assert.AreEqual(s[3][4].x, 80);
        Assert.AreEqual(s[3][4].y, 80);
        Assert.AreEqual(s[4][0].x, 20);
        Assert.AreEqual(s[4][0].y, 0);
        Assert.AreEqual(s[4][1].x, 80);
        Assert.AreEqual(s[4][1].y, 0);
        Assert.AreEqual(s[4][2].x, 80);
        Assert.AreEqual(s[4][2].y, 110);
        Assert.AreEqual(s[4][3].x, 20);
        Assert.AreEqual(s[4][3].y, 110);
        Assert.AreEqual(s[4][4].x, 20);
        Assert.AreEqual(s[4][4].y, 0);
        
        PathsD ns = GeoWrangler.rectangular_decomposition(ref abort, nS, maxRayLength: rayLength);
        writeToLayout("ns", nS, ns);
        Assert.AreEqual(ns.Count, 5);
        Assert.AreEqual(ns[0][0].x, 0);
        Assert.AreEqual(ns[0][0].y, -150);
        Assert.AreEqual(ns[0][1].x, 20);
        Assert.AreEqual(ns[0][1].y, -150);
        Assert.AreEqual(ns[0][2].x, 20);
        Assert.AreEqual(ns[0][2].y, -90);
        Assert.AreEqual(ns[0][3].x, 0);
        Assert.AreEqual(ns[0][3].y, -90);
        Assert.AreEqual(ns[0][4].x, 0);
        Assert.AreEqual(ns[0][4].y, -150);
        Assert.AreEqual(ns[1][0].x, 0);
        Assert.AreEqual(ns[1][0].y, -200);
        Assert.AreEqual(ns[1][1].x, 20);
        Assert.AreEqual(ns[1][1].y, -200);
        Assert.AreEqual(ns[1][2].x, 20);
        Assert.AreEqual(ns[1][2].y, -180);
        Assert.AreEqual(ns[1][3].x, 0);
        Assert.AreEqual(ns[1][3].y, -180);
        Assert.AreEqual(ns[1][4].x, 0);
        Assert.AreEqual(ns[1][4].y, -200);
        Assert.AreEqual(ns[2][0].x, 80);
        Assert.AreEqual(ns[2][0].y, -200);
        Assert.AreEqual(ns[2][1].x, 100);
        Assert.AreEqual(ns[2][1].y, -200);
        Assert.AreEqual(ns[2][2].x, 100);
        Assert.AreEqual(ns[2][2].y, -140);
        Assert.AreEqual(ns[2][3].x, 80);
        Assert.AreEqual(ns[2][3].y, -140);
        Assert.AreEqual(ns[2][4].x, 80);
        Assert.AreEqual(ns[2][4].y, -200);
        Assert.AreEqual(ns[3][0].x, 80);
        Assert.AreEqual(ns[3][0].y, -120);
        Assert.AreEqual(ns[3][1].x, 100);
        Assert.AreEqual(ns[3][1].y, -120);
        Assert.AreEqual(ns[3][2].x, 100);
        Assert.AreEqual(ns[3][2].y, -90);
        Assert.AreEqual(ns[3][3].x, 80);
        Assert.AreEqual(ns[3][3].y, -90);
        Assert.AreEqual(ns[3][4].x, 80);
        Assert.AreEqual(ns[3][4].y, -120);
        Assert.AreEqual(ns[4][0].x, 20);
        Assert.AreEqual(ns[4][0].y, -200);
        Assert.AreEqual(ns[4][1].x, 80);
        Assert.AreEqual(ns[4][1].y, -200);
        Assert.AreEqual(ns[4][2].x, 80);
        Assert.AreEqual(ns[4][2].y, -90);
        Assert.AreEqual(ns[4][3].x, 20);
        Assert.AreEqual(ns[4][3].y, -90);
        Assert.AreEqual(ns[4][4].x, 20);
        Assert.AreEqual(ns[4][4].y, -200);
        
        PathsD c1 = GeoWrangler.rectangular_decomposition(ref abort, C1, maxRayLength: rayLength);
        writeToLayout("c1", C1, c1);
        Assert.AreEqual(c1.Count, 17);
        // Use area because result is complex and hash not reliable due to floats.
        Assert.AreEqual(Clipper.Area(c1), 4920);

        PathsD c2 = GeoWrangler.rectangular_decomposition(ref abort, C2, maxRayLength: rayLength);
        writeToLayout("c2", C2, c2);
        Assert.AreEqual(c2.Count, 81);
        // Use area because result is complex and hash not reliable due to floats.
        Assert.AreEqual(Clipper.Area(c2), 24600);
        
        PathsD c3 = GeoWrangler.rectangular_decomposition(ref abort, C3, maxRayLength: rayLength);
        writeToLayout("c3", C3, c3);
        Assert.AreEqual(c3.Count, 13);
        // Use area because result is complex and hash not reliable due to floats.
        Assert.AreEqual(Clipper.Area(c3), 5424);
        
        PathsD c10r15 = GeoWrangler.rectangular_decomposition(ref abort, C10R15, maxRayLength: rayLength);
        writeToLayout("c10r15", C10R15, c10r15);
        Assert.AreEqual(c10r15.Count, 1);
        Assert.LessOrEqual(Math.Abs(Clipper.Area(c10r15))-14742059.915, 0.001);

        PathsD s1 = GeoWrangler.rectangular_decomposition(ref abort, S1, maxRayLength: rayLength);
        writeToLayout("s1", S1, s1);
        Assert.AreEqual(s1.Count, 4);
        Assert.AreEqual(s1[0][0].x, -50);
        Assert.AreEqual(s1[0][0].y, -50);
        Assert.AreEqual(s1[0][1].x, 0);
        Assert.AreEqual(s1[0][1].y, -50);
        Assert.AreEqual(s1[0][2].x, 0);
        Assert.AreEqual(s1[0][2].y, 0);
        Assert.AreEqual(s1[0][3].x, -50);
        Assert.AreEqual(s1[0][3].y, 0);
        Assert.AreEqual(s1[0][4].x, -50);
        Assert.AreEqual(s1[0][4].y, -50);
        Assert.AreEqual(s1[1][0].x, 0);
        Assert.AreEqual(s1[1][0].y, -50);
        Assert.AreEqual(s1[1][1].x, 100);
        Assert.AreEqual(s1[1][1].y, -50);
        Assert.AreEqual(s1[1][2].x, 100);
        Assert.AreEqual(s1[1][2].y, 120);
        Assert.AreEqual(s1[1][3].x, 0);
        Assert.AreEqual(s1[1][3].y, 120);
        Assert.AreEqual(s1[1][4].x, 0);
        Assert.AreEqual(s1[1][4].y, -50);
        Assert.AreEqual(s1[2][0].x, 150);
        Assert.AreEqual(s1[2][0].y, -50);
        Assert.AreEqual(s1[2][1].x, 200);
        Assert.AreEqual(s1[2][1].y, -50);
        Assert.AreEqual(s1[2][2].x, 200);
        Assert.AreEqual(s1[2][2].y, 300);
        Assert.AreEqual(s1[2][3].x, 150);
        Assert.AreEqual(s1[2][3].y, 300);
        Assert.AreEqual(s1[2][4].x, 150);
        Assert.AreEqual(s1[2][4].y, -50);
        Assert.AreEqual(s1[3][0].x, 100);
        Assert.AreEqual(s1[3][0].y, -50);
        Assert.AreEqual(s1[3][1].x, 150);
        Assert.AreEqual(s1[3][1].y, -50);
        Assert.AreEqual(s1[3][2].x, 150);
        Assert.AreEqual(s1[3][2].y, 200);
        Assert.AreEqual(s1[3][3].x, 100);
        Assert.AreEqual(s1[3][3].y, 200);
        Assert.AreEqual(s1[3][4].x, 100);
        Assert.AreEqual(s1[3][4].y, -50);
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
        
        ClipperD c = new(Constants.roundingDecimalPrecision);
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
        done = Clipper.ScalePaths(done, 1000);
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        bool abort = false;

        Console.WriteLine("  Decomposition (vertical)....");
        sw.Restart();
        PathsD ns = GeoWrangler.rectangular_decomposition(ref abort, done, maxRayLength: rayLength);
        
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Writing....");
        writeToLayout("complex", done[0], ns);

        double area = Clipper.Area(ns);
        Assert.AreEqual(ns.Count, 200);
        Assert.AreEqual(Math.Abs(area), 221320);
        
        Console.WriteLine("  Decomposition (horizontal)....");
        sw.Restart();
        ns = GeoWrangler.rectangular_decomposition(ref abort, done, maxRayLength: rayLength, vertical: false);

        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Writing....");
        writeToLayout("complex_horizontal", done[0], ns);

        area = Clipper.Area(ns);
        Assert.AreEqual(ns.Count, 607);
        Assert.AreEqual(Math.Abs(area), 221320);

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

        partFour_do(points_1, "complex_loop", 401, 401, 325380.0);

        Console.WriteLine(" Part 2....");
        Console.WriteLine("  Preparing....");
        sw.Start();
        PathD points_2 = TestGeometry.ortho_fractal_2();
        points_2 = GeoWrangler.close(points_2);
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");
        partFour_do(points_2, "complex_loop_rot", 401, 401, 325380);
    }

    private static void partFour_do(PathD points, string baseString, int expectedCountV, int expectedCountH, double expectedArea)
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
        
        Assert.AreEqual(decompOut.Count, expectedCountV);
        Assert.AreEqual(Clipper.Area(decompOut), expectedArea);

        Console.WriteLine("  Decomposition (horizontal)....");
        sw.Restart();
        decompOut = GeoWrangler.rectangular_decomposition(ref abort, toDecomp,
            maxRayLength: (long)Math.Max(Math.Abs(dist.x), Math.Abs(dist.y)) * 1, vertical: !vertical);
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Writing....");
        writeToLayout(baseString + "_horizontal", points, decompOut);

        Assert.AreEqual(decompOut.Count, expectedCountH);
        Assert.AreEqual(Clipper.Area(decompOut), expectedArea);
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
        
        PathsD out_decomp = new();
        for (int i = 0; i < polydata.Count; i++)
        {
         PathD points = new (polydata[i]);
         points = GeoWrangler.removeDuplicates(points);
         points = GeoWrangler.stripCollinear(points);
         points = GeoWrangler.clockwiseAndReorderXY(points);
         
         PathD toDecomp = GeoWrangler.makeKeyHole(GeoWrangler.sliverGapRemoval(points), reverseEval:false, biDirectionalEval:false)[0];
         PathD  bounds = GeoWrangler.getBounds(toDecomp);
         PointD dist = GeoWrangler.distanceBetweenPoints_point(bounds[0], bounds[1]);

         PathsD decompOut = GeoWrangler.rectangular_decomposition(ref abort, toDecomp,
          maxRayLength: (long) Math.Max(Math.Abs(dist.x), Math.Abs(dist.y)), vertical: vertical);
         
         out_decomp.AddRange(decompOut);
        }
        
        Assert.AreEqual(out_decomp.Count, 19);
        Assert.AreEqual(Clipper.Area(out_decomp), 13843.0);
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