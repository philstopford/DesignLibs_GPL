using System.Diagnostics;
using Clipper2Lib;
using geoCoreLib;
using geoWrangler;
using PartitionTestGeometrySource;

namespace UnitTests;

public class DecompositionTests
{
    private static string root_loc = "/d/development/decomposition_out/";
    // [SetUp]
    public static void DecompositionSetup()
    {
        Console.WriteLine("Part One");
        partOne();

        Console.WriteLine("Part Two");
        partTwo();

        Console.WriteLine("Part Three (complex decomposition tests, lengthy)");
        partThree_moderatelycomplex();
        partThree_morecomplex();
        partThree_evenmorecomplex();
        partThree_canarycomplex();
        partThree_gainingacomplex();
        partThree_extremelycomplex();

        Console.WriteLine("Part Four (takes less time than part three)");
        partFour_1();
        partFour_2();

        Console.WriteLine("Part Five");
        partFive();
    }

    [Test]
    public static void partOne()
    {
        // L
        PathD L = TestGeometry.getL();
        PathD rL = TestGeometry.getRL();

        // Reversed orientation.
        PathD L_ccw = new(L);
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
        // Sign change expected.
        Assert.That(-Clipper.Area(l), Is.EqualTo(Clipper.Area(L)));
        Assert.That(l.Count, Is.EqualTo(2));
        Assert.That(l[0][0].x, Is.EqualTo(0));
        Assert.That(l[0][0].y, Is.EqualTo(0));
        Assert.That(l[0][1].x, Is.EqualTo(10));
        Assert.That(l[0][1].y, Is.EqualTo(0));
        Assert.That(l[0][2].x, Is.EqualTo(10));
        Assert.That(l[0][2].y, Is.EqualTo(50));
        Assert.That(l[0][3].x, Is.EqualTo(0));
        Assert.That(l[0][3].y, Is.EqualTo(50));
        Assert.That(l[0][4].x, Is.EqualTo(0));
        Assert.That(l[0][4].y, Is.EqualTo(0));
        Assert.That(l[1][0].x, Is.EqualTo(10));
        Assert.That(l[1][0].y, Is.EqualTo(0));
        Assert.That(l[1][1].x, Is.EqualTo(60));
        Assert.That(l[1][1].y, Is.EqualTo(0));
        Assert.That(l[1][2].x, Is.EqualTo(60));
        Assert.That(l[1][2].y, Is.EqualTo(20));
        Assert.That(l[1][3].x, Is.EqualTo(10));
        Assert.That(l[1][3].y, Is.EqualTo(20));
        Assert.That(l[1][4].x, Is.EqualTo(10));
        Assert.That(l[1][4].y, Is.EqualTo(0));

        PathsD lccw = GeoWrangler.rectangular_decomposition(ref abort, L_ccw, maxRayLength: rayLength);
        writeToLayout("lccw", L_ccw, lccw);
        Assert.That(Clipper.Area(lccw), Is.EqualTo(Clipper.Area(L_ccw)));
        Assert.That(lccw.Count, Is.EqualTo(2));
        Assert.That(lccw[0][0].x, Is.EqualTo(0));
        Assert.That(lccw[0][0].y, Is.EqualTo(0));
        Assert.That(lccw[0][1].x, Is.EqualTo(10));
        Assert.That(lccw[0][1].y, Is.EqualTo(0));
        Assert.That(lccw[0][2].x, Is.EqualTo(10));
        Assert.That(lccw[0][2].y, Is.EqualTo(50));
        Assert.That(lccw[0][3].x, Is.EqualTo(0));
        Assert.That(lccw[0][3].y, Is.EqualTo(50));
        Assert.That(lccw[0][4].x, Is.EqualTo(0));
        Assert.That(lccw[0][4].y, Is.EqualTo(0));
        Assert.That(lccw[1][0].x, Is.EqualTo(10));
        Assert.That(lccw[1][0].y, Is.EqualTo(0));
        Assert.That(lccw[1][1].x, Is.EqualTo(60));
        Assert.That(lccw[1][1].y, Is.EqualTo(0));
        Assert.That(lccw[1][2].x, Is.EqualTo(60));
        Assert.That(lccw[1][2].y, Is.EqualTo(20));
        Assert.That(lccw[1][3].x, Is.EqualTo(10));
        Assert.That(lccw[1][3].y, Is.EqualTo(20));
        Assert.That(lccw[1][4].x, Is.EqualTo(10));
        Assert.That(lccw[1][4].y, Is.EqualTo(0));

        PathsD rl = GeoWrangler.rectangular_decomposition(ref abort, rL, maxRayLength: rayLength);
        writeToLayout("rl", rL, rl);
        // Sign change expected.
        Assert.That(-Clipper.Area(rl), Is.EqualTo(Clipper.Area(rL)));
        Assert.That(rl.Count, Is.EqualTo(2));
        Assert.That(rl[0][0].x, Is.EqualTo(10));
        Assert.That(rl[0][0].y, Is.EqualTo(0));
        Assert.That(rl[0][1].x, Is.EqualTo(60));
        Assert.That(rl[0][1].y, Is.EqualTo(0));
        Assert.That(rl[0][2].x, Is.EqualTo(60));
        Assert.That(rl[0][2].y, Is.EqualTo(50));
        Assert.That(rl[0][3].x, Is.EqualTo(10));
        Assert.That(rl[0][3].y, Is.EqualTo(50));
        Assert.That(rl[0][4].x, Is.EqualTo(10));
        Assert.That(rl[0][4].y, Is.EqualTo(0));
        Assert.That(rl[1][0].x, Is.EqualTo(0));
        Assert.That(rl[1][0].y, Is.EqualTo(0));
        Assert.That(rl[1][1].x, Is.EqualTo(10));
        Assert.That(rl[1][1].y, Is.EqualTo(0));
        Assert.That(rl[1][2].x, Is.EqualTo(10));
        Assert.That(rl[1][2].y, Is.EqualTo(20));
        Assert.That(rl[1][3].x, Is.EqualTo(0));
        Assert.That(rl[1][3].y, Is.EqualTo(20));
        Assert.That(rl[1][4].x, Is.EqualTo(0));
        Assert.That(rl[1][4].y, Is.EqualTo(0));

        PathsD u = GeoWrangler.rectangular_decomposition(ref abort, U, maxRayLength: rayLength);
        writeToLayout("u", U, u);
        // Sign change expected.
        Assert.That(-Clipper.Area(u), Is.EqualTo(Clipper.Area(U)));
        Assert.That(u.Count, Is.EqualTo(3));
        Assert.That(u[0][0].x, Is.EqualTo(0));
        Assert.That(u[0][0].y, Is.EqualTo(0));
        Assert.That(u[0][1].x, Is.EqualTo(10));
        Assert.That(u[0][1].y, Is.EqualTo(0));
        Assert.That(u[0][2].x, Is.EqualTo(10));
        Assert.That(u[0][2].y, Is.EqualTo(50));
        Assert.That(u[0][3].x, Is.EqualTo(0));
        Assert.That(u[0][3].y, Is.EqualTo(50));
        Assert.That(u[0][4].x, Is.EqualTo(0));
        Assert.That(u[0][4].y, Is.EqualTo(0));
        Assert.That(u[1][0].x, Is.EqualTo(60));
        Assert.That(u[1][0].y, Is.EqualTo(0));
        Assert.That(u[1][1].x, Is.EqualTo(120));
        Assert.That(u[1][1].y, Is.EqualTo(0));
        Assert.That(u[1][2].x, Is.EqualTo(120));
        Assert.That(u[1][2].y, Is.EqualTo(80));
        Assert.That(u[1][3].x, Is.EqualTo(60));
        Assert.That(u[1][3].y, Is.EqualTo(80));
        Assert.That(u[1][4].x, Is.EqualTo(60));
        Assert.That(u[1][4].y, Is.EqualTo(0));
        Assert.That(u[2][0].x, Is.EqualTo(10));
        Assert.That(u[2][0].y, Is.EqualTo(0));
        Assert.That(u[2][1].x, Is.EqualTo(60));
        Assert.That(u[2][1].y, Is.EqualTo(0));
        Assert.That(u[2][2].x, Is.EqualTo(60));
        Assert.That(u[2][2].y, Is.EqualTo(20));
        Assert.That(u[2][3].x, Is.EqualTo(10));
        Assert.That(u[2][3].y, Is.EqualTo(20));
        Assert.That(u[2][4].x, Is.EqualTo(10));
        Assert.That(u[2][4].y, Is.EqualTo(0));

        PathsD t = GeoWrangler.rectangular_decomposition(ref abort, T, maxRayLength: rayLength);
        writeToLayout("t", T, t);
        // Sign change expected.
        Assert.That(-Clipper.Area(t), Is.EqualTo(Clipper.Area(T)));
        Assert.That(t.Count, Is.EqualTo(3));
        Assert.That(t[0][0].x, Is.EqualTo(60));
        Assert.That(t[0][0].y, Is.EqualTo(50));
        Assert.That(t[0][1].x, Is.EqualTo(80));
        Assert.That(t[0][1].y, Is.EqualTo(50));
        Assert.That(t[0][2].x, Is.EqualTo(80));
        Assert.That(t[0][2].y, Is.EqualTo(80));
        Assert.That(t[0][3].x, Is.EqualTo(60));
        Assert.That(t[0][3].y, Is.EqualTo(80));
        Assert.That(t[0][4].x, Is.EqualTo(60));
        Assert.That(t[0][4].y, Is.EqualTo(50));
        Assert.That(t[1][0].x, Is.EqualTo(40));
        Assert.That(t[1][0].y, Is.EqualTo(0));
        Assert.That(t[1][1].x, Is.EqualTo(60));
        Assert.That(t[1][1].y, Is.EqualTo(0));
        Assert.That(t[1][2].x, Is.EqualTo(60));
        Assert.That(t[1][2].y, Is.EqualTo(80));
        Assert.That(t[1][3].x, Is.EqualTo(40));
        Assert.That(t[1][3].y, Is.EqualTo(80));
        Assert.That(t[1][4].x, Is.EqualTo(40));
        Assert.That(t[1][4].y, Is.EqualTo(0));
        Assert.That(t[2][0].x, Is.EqualTo(0));
        Assert.That(t[2][0].y, Is.EqualTo(50));
        Assert.That(t[2][1].x, Is.EqualTo(40));
        Assert.That(t[2][1].y, Is.EqualTo(50));
        Assert.That(t[2][2].x, Is.EqualTo(40));
        Assert.That(t[2][2].y, Is.EqualTo(80));
        Assert.That(t[2][3].x, Is.EqualTo(0));
        Assert.That(t[2][3].y, Is.EqualTo(80));
        Assert.That(t[2][4].x, Is.EqualTo(0));
        Assert.That(t[2][4].y, Is.EqualTo(50));

        PathsD x = GeoWrangler.rectangular_decomposition(ref abort, X, maxRayLength: rayLength);
        writeToLayout("x", X, x);
        // Sign change expected.
        Assert.That(-Clipper.Area(x), Is.EqualTo(Clipper.Area(X)));
        Assert.That(x.Count, Is.EqualTo(3));
        Assert.That(x[0][0].x, Is.EqualTo(0));
        Assert.That(x[0][0].y, Is.EqualTo(50));
        Assert.That(x[0][1].x, Is.EqualTo(60));
        Assert.That(x[0][1].y, Is.EqualTo(50));
        Assert.That(x[0][2].x, Is.EqualTo(60));
        Assert.That(x[0][2].y, Is.EqualTo(80));
        Assert.That(x[0][3].x, Is.EqualTo(0));
        Assert.That(x[0][3].y, Is.EqualTo(80));
        Assert.That(x[0][4].x, Is.EqualTo(0));
        Assert.That(x[0][4].y, Is.EqualTo(50));
        Assert.That(x[1][0].x, Is.EqualTo(60));
        Assert.That(x[1][0].y, Is.EqualTo(20));
        Assert.That(x[1][1].x, Is.EqualTo(80));
        Assert.That(x[1][1].y, Is.EqualTo(20));
        Assert.That(x[1][2].x, Is.EqualTo(80));
        Assert.That(x[1][2].y, Is.EqualTo(100));
        Assert.That(x[1][3].x, Is.EqualTo(60));
        Assert.That(x[1][3].y, Is.EqualTo(100));
        Assert.That(x[1][4].x, Is.EqualTo(60));
        Assert.That(x[1][4].y, Is.EqualTo(20));
        Assert.That(x[2][0].x, Is.EqualTo(80));
        Assert.That(x[2][0].y, Is.EqualTo(50));
        Assert.That(x[2][1].x, Is.EqualTo(100));
        Assert.That(x[2][1].y, Is.EqualTo(50));
        Assert.That(x[2][2].x, Is.EqualTo(100));
        Assert.That(x[2][2].y, Is.EqualTo(80));
        Assert.That(x[2][3].x, Is.EqualTo(80));
        Assert.That(x[2][3].y, Is.EqualTo(80));
        Assert.That(x[2][4].x, Is.EqualTo(80));
        Assert.That(x[2][4].y, Is.EqualTo(50));

        PathsD s = GeoWrangler.rectangular_decomposition(ref abort, S, maxRayLength: rayLength);
        writeToLayout("s", S, s);
        // Sign change expected.
        Assert.That(-Clipper.Area(s), Is.EqualTo(Clipper.Area(S)));
        Assert.That(s.Count, Is.EqualTo(5));
        Assert.That(s[0][0].x, Is.EqualTo(0));
        Assert.That(s[0][0].y, Is.EqualTo(50));
        Assert.That(s[0][1].x, Is.EqualTo(20));
        Assert.That(s[0][1].y, Is.EqualTo(50));
        Assert.That(s[0][2].x, Is.EqualTo(20));
        Assert.That(s[0][2].y, Is.EqualTo(110));
        Assert.That(s[0][3].x, Is.EqualTo(0));
        Assert.That(s[0][3].y, Is.EqualTo(110));
        Assert.That(s[0][4].x, Is.EqualTo(0));
        Assert.That(s[0][4].y, Is.EqualTo(50));
        Assert.That(s[1][0].x, Is.EqualTo(0));
        Assert.That(s[1][0].y, Is.EqualTo(0));
        Assert.That(s[1][1].x, Is.EqualTo(20));
        Assert.That(s[1][1].y, Is.EqualTo(0));
        Assert.That(s[1][2].x, Is.EqualTo(20));
        Assert.That(s[1][2].y, Is.EqualTo(20));
        Assert.That(s[1][3].x, Is.EqualTo(0));
        Assert.That(s[1][3].y, Is.EqualTo(20));
        Assert.That(s[1][4].x, Is.EqualTo(0));
        Assert.That(s[1][4].y, Is.EqualTo(0));
        Assert.That(s[2][0].x, Is.EqualTo(80));
        Assert.That(s[2][0].y, Is.EqualTo(0));
        Assert.That(s[2][1].x, Is.EqualTo(100));
        Assert.That(s[2][1].y, Is.EqualTo(0));
        Assert.That(s[2][2].x, Is.EqualTo(100));
        Assert.That(s[2][2].y, Is.EqualTo(60));
        Assert.That(s[2][3].x, Is.EqualTo(80));
        Assert.That(s[2][3].y, Is.EqualTo(60));
        Assert.That(s[2][4].x, Is.EqualTo(80));
        Assert.That(s[2][4].y, Is.EqualTo(0));
        Assert.That(s[3][0].x, Is.EqualTo(80));
        Assert.That(s[3][0].y, Is.EqualTo(80));
        Assert.That(s[3][1].x, Is.EqualTo(100));
        Assert.That(s[3][1].y, Is.EqualTo(80));
        Assert.That(s[3][2].x, Is.EqualTo(100));
        Assert.That(s[3][2].y, Is.EqualTo(110));
        Assert.That(s[3][3].x, Is.EqualTo(80));
        Assert.That(s[3][3].y, Is.EqualTo(110));
        Assert.That(s[3][4].x, Is.EqualTo(80));
        Assert.That(s[3][4].y, Is.EqualTo(80));
        Assert.That(s[4][0].x, Is.EqualTo(20));
        Assert.That(s[4][0].y, Is.EqualTo(0));
        Assert.That(s[4][1].x, Is.EqualTo(80));
        Assert.That(s[4][1].y, Is.EqualTo(0));
        Assert.That(s[4][2].x, Is.EqualTo(80));
        Assert.That(s[4][2].y, Is.EqualTo(110));
        Assert.That(s[4][3].x, Is.EqualTo(20));
        Assert.That(s[4][3].y, Is.EqualTo(110));
        Assert.That(s[4][4].x, Is.EqualTo(20));
        Assert.That(s[4][4].y, Is.EqualTo(0));

        PathsD ns = GeoWrangler.rectangular_decomposition(ref abort, nS, maxRayLength: rayLength);
        writeToLayout("ns", nS, ns);
        // Sign change expected.
        Assert.That(-Clipper.Area(ns), Is.EqualTo(Clipper.Area(nS)));
        Assert.That(ns.Count, Is.EqualTo(5));
        Assert.That(ns[0][0].x, Is.EqualTo(0));
        Assert.That(ns[0][0].y, Is.EqualTo(-150));
        Assert.That(ns[0][1].x, Is.EqualTo(20));
        Assert.That(ns[0][1].y, Is.EqualTo(-150));
        Assert.That(ns[0][2].x, Is.EqualTo(20));
        Assert.That(ns[0][2].y, Is.EqualTo(-90));
        Assert.That(ns[0][3].x, Is.EqualTo(0));
        Assert.That(ns[0][3].y, Is.EqualTo(-90));
        Assert.That(ns[0][4].x, Is.EqualTo(0));
        Assert.That(ns[0][4].y, Is.EqualTo(-150));
        Assert.That(ns[1][0].x, Is.EqualTo(0));
        Assert.That(ns[1][0].y, Is.EqualTo(-200));
        Assert.That(ns[1][1].x, Is.EqualTo(20));
        Assert.That(ns[1][1].y, Is.EqualTo(-200));
        Assert.That(ns[1][2].x, Is.EqualTo(20));
        Assert.That(ns[1][2].y, Is.EqualTo(-180));
        Assert.That(ns[1][3].x, Is.EqualTo(0));
        Assert.That(ns[1][3].y, Is.EqualTo(-180));
        Assert.That(ns[1][4].x, Is.EqualTo(0));
        Assert.That(ns[1][4].y, Is.EqualTo(-200));
        Assert.That(ns[2][0].x, Is.EqualTo(80));
        Assert.That(ns[2][0].y, Is.EqualTo(-200));
        Assert.That(ns[2][1].x, Is.EqualTo(100));
        Assert.That(ns[2][1].y, Is.EqualTo(-200));
        Assert.That(ns[2][2].x, Is.EqualTo(100));
        Assert.That(ns[2][2].y, Is.EqualTo(-140));
        Assert.That(ns[2][3].x, Is.EqualTo(80));
        Assert.That(ns[2][3].y, Is.EqualTo(-140));
        Assert.That(ns[2][4].x, Is.EqualTo(80));
        Assert.That(ns[2][4].y, Is.EqualTo(-200));
        Assert.That(ns[3][0].x, Is.EqualTo(80));
        Assert.That(ns[3][0].y, Is.EqualTo(-120));
        Assert.That(ns[3][1].x, Is.EqualTo(100));
        Assert.That(ns[3][1].y, Is.EqualTo(-120));
        Assert.That(ns[3][2].x, Is.EqualTo(100));
        Assert.That(ns[3][2].y, Is.EqualTo(-90));
        Assert.That(ns[3][3].x, Is.EqualTo(80));
        Assert.That(ns[3][3].y, Is.EqualTo(-90));
        Assert.That(ns[3][4].x, Is.EqualTo(80));
        Assert.That(ns[3][4].y, Is.EqualTo(-120));
        Assert.That(ns[4][0].x, Is.EqualTo(20));
        Assert.That(ns[4][0].y, Is.EqualTo(-200));
        Assert.That(ns[4][1].x, Is.EqualTo(80));
        Assert.That(ns[4][1].y, Is.EqualTo(-200));
        Assert.That(ns[4][2].x, Is.EqualTo(80));
        Assert.That(ns[4][2].y, Is.EqualTo(-90));
        Assert.That(ns[4][3].x, Is.EqualTo(20));
        Assert.That(ns[4][3].y, Is.EqualTo(-90));
        Assert.That(ns[4][4].x, Is.EqualTo(20));
        Assert.That(ns[4][4].y, Is.EqualTo(-200));

        PathsD c1 = GeoWrangler.rectangular_decomposition(ref abort, C1, maxRayLength: rayLength);
        writeToLayout("c1", C1, c1);
        // Sign change expected.
        Assert.That(-Clipper.Area(c1), Is.EqualTo(Clipper.Area(C1)));
        Assert.That(c1.Count, Is.EqualTo(17));
        // Use area because result is complex and hash not reliable due to floats.
        Assert.That(Clipper.Area(c1), Is.EqualTo(4920));

        PathsD c2 = GeoWrangler.rectangular_decomposition(ref abort, C2, maxRayLength: rayLength);
        writeToLayout("c2", C2, c2);
        // Sign change expected.
        Assert.That(-Clipper.Area(c2), Is.EqualTo(Clipper.Area(C2)));
        Assert.That(c2.Count, Is.EqualTo(81));
        // Use area because result is complex and hash not reliable due to floats.
        Assert.That(Clipper.Area(c2), Is.EqualTo(24600));

        PathsD c3 = GeoWrangler.rectangular_decomposition(ref abort, C3, maxRayLength: rayLength);
        writeToLayout("c3", C3, c3);
        // Sign change expected.
        Assert.That(-Clipper.Area(c3), Is.EqualTo(Clipper.Area(C3)));
        Assert.That(c3.Count, Is.EqualTo(13));
        // Use area because result is complex and hash not reliable due to floats.
        Assert.That(Clipper.Area(c3), Is.EqualTo(5424));

        PathsD c10r15 = GeoWrangler.rectangular_decomposition(ref abort, C10R15, maxRayLength: rayLength);
        writeToLayout("c10r15", C10R15, c10r15);
        // Sign change expected.
        Assert.That(Clipper.Area(c10r15), Is.EqualTo(Clipper.Area(C10R15)));
        Assert.That(c10r15.Count, Is.EqualTo(1));
        Assert.That(Clipper.Area(c10r15), Is.EqualTo(-14742059.915).Within(0.001));

        PathsD s1 = GeoWrangler.rectangular_decomposition(ref abort, S1, maxRayLength: rayLength);
        writeToLayout("s1", S1, s1);
        // Sign change expected.
        Assert.That(-Clipper.Area(s1), Is.EqualTo(Clipper.Area(S1)));
        Assert.That(s1.Count, Is.EqualTo(4));
        Assert.That(s1[0][0].x, Is.EqualTo(-50));
        Assert.That(s1[0][0].y, Is.EqualTo(-50));
        Assert.That(s1[0][1].x, Is.EqualTo(0));
        Assert.That(s1[0][1].y, Is.EqualTo(-50));
        Assert.That(s1[0][2].x, Is.EqualTo(0));
        Assert.That(s1[0][2].y, Is.EqualTo(0));
        Assert.That(s1[0][3].x, Is.EqualTo(-50));
        Assert.That(s1[0][3].y, Is.EqualTo(0));
        Assert.That(s1[0][4].x, Is.EqualTo(-50));
        Assert.That(s1[0][4].y, Is.EqualTo(-50));
        Assert.That(s1[1][0].x, Is.EqualTo(0));
        Assert.That(s1[1][0].y, Is.EqualTo(-50));
        Assert.That(s1[1][1].x, Is.EqualTo(100));
        Assert.That(s1[1][1].y, Is.EqualTo(-50));
        Assert.That(s1[1][2].x, Is.EqualTo(100));
        Assert.That(s1[1][2].y, Is.EqualTo(120));
        Assert.That(s1[1][3].x, Is.EqualTo(0));
        Assert.That(s1[1][3].y, Is.EqualTo(120));
        Assert.That(s1[1][4].x, Is.EqualTo(0));
        Assert.That(s1[1][4].y, Is.EqualTo(-50));
        Assert.That(s1[2][0].x, Is.EqualTo(150));
        Assert.That(s1[2][0].y, Is.EqualTo(-50));
        Assert.That(s1[2][1].x, Is.EqualTo(200));
        Assert.That(s1[2][1].y, Is.EqualTo(-50));
        Assert.That(s1[2][2].x, Is.EqualTo(200));
        Assert.That(s1[2][2].y, Is.EqualTo(300));
        Assert.That(s1[2][3].x, Is.EqualTo(150));
        Assert.That(s1[2][3].y, Is.EqualTo(300));
        Assert.That(s1[2][4].x, Is.EqualTo(150));
        Assert.That(s1[2][4].y, Is.EqualTo(-50));
        Assert.That(s1[3][0].x, Is.EqualTo(100));
        Assert.That(s1[3][0].y, Is.EqualTo(-50));
        Assert.That(s1[3][1].x, Is.EqualTo(150));
        Assert.That(s1[3][1].y, Is.EqualTo(-50));
        Assert.That(s1[3][2].x, Is.EqualTo(150));
        Assert.That(s1[3][2].y, Is.EqualTo(200));
        Assert.That(s1[3][3].x, Is.EqualTo(100));
        Assert.That(s1[3][3].y, Is.EqualTo(200));
        Assert.That(s1[3][4].x, Is.EqualTo(100));
        Assert.That(s1[3][4].y, Is.EqualTo(-50));
    }

    [Test]
    public static void partTwo()
    {
        PathsD incoming = new();
        PathD lPieces = new()
        {
            new(0.00000, 0.00000),
            new(0.00000, 0.05000),
            new(0.01000, 0.05000),
            new(0.01000, 0.00000),
            new(0.00000, 0.00000)
        };


        PathD lPiece2 = new()
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

    [Test]
    public static void partThree_moderatelycomplex()
    {
        Stopwatch sw = new();
        int rayLength = 1000;

        Console.WriteLine(" ModeratelyComplex");
        Console.WriteLine("  Preparing....");
        sw.Start();

        PathD poly = TestGeometry.moderatelycomplex();

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
        writeToLayout("moderatelycomplex", done[0], ns);

        Assert.That(ns.Count, Is.EqualTo(81));
        // Sign change expected.
        Assert.That(-Clipper.Area(ns), Is.EqualTo(Clipper.Area(done)));

        Console.WriteLine("  Decomposition (horizontal)....");
        sw.Restart();
        ns = GeoWrangler.rectangular_decomposition(ref abort, done, maxRayLength: rayLength, vertical: false);

        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Writing....");
        writeToLayout("moderatelycomplex_horizontal", done[0], ns);

        Assert.That(ns.Count, Is.EqualTo(81));
        // Sign change expected.
        Assert.That(-Clipper.Area(ns), Is.EqualTo(Clipper.Area(done)));

        Console.WriteLine("  Done.");
    }

    [Test]
    public static void partThree_morecomplex()
    {
        Stopwatch sw = new();
        int rayLength = 1000;

        Console.WriteLine(" MoreComplex");
        Console.WriteLine("  Preparing....");
        sw.Start();

        PathD poly = TestGeometry.morecomplex();

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
        writeToLayout("morecomplex", done[0], ns);

        Assert.That(ns.Count, Is.EqualTo(161));
        // Sign change expected.
        Assert.That(-Clipper.Area(ns), Is.EqualTo(Clipper.Area(done)));

        Console.WriteLine("  Decomposition (horizontal)....");
        sw.Restart();
        ns = GeoWrangler.rectangular_decomposition(ref abort, done, maxRayLength: rayLength, vertical: false);

        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Writing....");
        writeToLayout("morecomplex_horizontal", done[0], ns);

        Assert.That(ns.Count, Is.EqualTo(161));
        // Sign change expected.
        Assert.That(-Clipper.Area(ns), Is.EqualTo(Clipper.Area(done)));

        Console.WriteLine("  Done.");
    }
    [Test]
    public static void partThree_evenmorecomplex()
    {
        Stopwatch sw = new();
        int rayLength = 1000;

        Console.WriteLine(" EvenMoreComplex");
        Console.WriteLine("  Preparing....");
        sw.Start();

        PathD poly = TestGeometry.evenmorecomplex();

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
        writeToLayout("evenmorecomplex", done[0], ns);

        Assert.That(ns.Count, Is.EqualTo(241));
        // Sign change expected.
        Assert.That(-Clipper.Area(ns), Is.EqualTo(Clipper.Area(done)));

        Console.WriteLine("  Decomposition (horizontal)....");
        sw.Restart();
        ns = GeoWrangler.rectangular_decomposition(ref abort, done, maxRayLength: rayLength, vertical: false);

        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Writing....");
        writeToLayout("evenmorecomplex_horizontal", done[0], ns);

        Assert.That(ns.Count, Is.EqualTo(241));
        // Sign change expected.
        Assert.That(-Clipper.Area(ns), Is.EqualTo(Clipper.Area(done)));

        Console.WriteLine("  Done.");
    }

    // This test corresponds to a detected failure from coincident edges and was the 'minimal' case where it started.
    [Test]
    public static void partThree_canarycomplex()
    {
        Stopwatch sw = new();
        int rayLength = 1000;

        Console.WriteLine(" CanaryComplex");
        Console.WriteLine("  Preparing....");
        sw.Start();

        PathD canary = TestGeometry.canarycomplex();
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        bool abort = false;

        Console.WriteLine("  Decomposition (vertical)....");
        sw.Restart();
        PathsD test = GeoWrangler.rectangular_decomposition(ref abort, canary, maxRayLength: rayLength);
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Writing....");
        writeToLayout("canarycomplex", canary, test);
        Assert.That(test.Count, Is.EqualTo(381));
        Assert.That(-Clipper.Area(test), Is.EqualTo(Clipper.Area(canary)));

        Console.WriteLine("  Decomposition (horizontal)....");
        sw.Restart();
        test = GeoWrangler.rectangular_decomposition(ref abort, canary, maxRayLength: rayLength, vertical: false);
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Writing....");
        writeToLayout("canarycomplex_horizontal", canary, test);
        Assert.That(test.Count, Is.EqualTo(386));
        Assert.That(-Clipper.Area(test), Is.EqualTo(Clipper.Area(canary)));

        Console.WriteLine("  Done.");
    }

    // This test case decomposes properly, but seems to trigger an issue in GDS output where a box ends up zero height
    // in the horizontal decomposition case. OASIS is fine.
    [Test]
    public static void partThree_canarycomplex2()
    {
        Stopwatch sw = new();
        int rayLength = 1000;

        Console.WriteLine(" CanaryComplex2");
        Console.WriteLine("  Preparing....");
        sw.Start();

        bool abort = false;

        PathD canary = Clipper.MakePath(new double[]
        {
            /*
            0.06, 0.44,
            0.06, 0.45,
            0.07, 0.45,
            0.07, 0.46,
            0, 0.46,
            0, 0.47,
            0.079, 0.47,
            0.079, 0.459,
            0.08, 0.459,
            0.08, 0.45,
            0.09, 0.45,
            0.09, 0.44,
            0.06, 0.44
            */
            0.07,   0.45,
            0.07,   0.46,
            0.06,   0.46,
            0.06,   0.47,
            0.079,  0.47,
            0.079,  0.459,
            0.08,   0.459,
            0.08,   0.45,
        });

        canary = Clipper.ScalePath(canary, 1000);

        Console.WriteLine("  Decomposition (horizontal)....");
        sw.Restart();
        PathsD test = GeoWrangler.rectangular_decomposition(ref abort, canary, maxRayLength: rayLength, vertical: false);
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Writing....");
        writeToLayout("canarycomplex2_horizontal", canary, test);
        Console.WriteLine("  Done.");
    }

    [Test]
    public static void partThree_gainingacomplex()
    {
        Stopwatch sw = new();
        int rayLength = 1000;

        Console.WriteLine(" GainingAComplex");
        Console.WriteLine("  Preparing....");
        sw.Start();

        PathD poly = TestGeometry.gainingacomplex();

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
        writeToLayout("gainingacomplex", done[0], ns);

        Assert.That(ns.Count, Is.EqualTo(401));
        // Sign change expected.
        Assert.That(-Clipper.Area(ns), Is.EqualTo(Clipper.Area(done)));

        Console.WriteLine("  Decomposition (horizontal)....");
        sw.Restart();
        ns = GeoWrangler.rectangular_decomposition(ref abort, done, maxRayLength: rayLength, vertical: false);

        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Writing....");
        writeToLayout("gainingacomplex_horizontal", done[0], ns);

        Assert.That(ns.Count, Is.EqualTo(401));
        // Sign change expected.
        Assert.That(-Clipper.Area(ns), Is.EqualTo(Clipper.Area(done)));

        Console.WriteLine("  Done.");
    }
    [Test]
    public static void partThree_extremelycomplex()
    {
        Stopwatch sw = new();
        int rayLength = 1000;

        Console.WriteLine(" ExtremelyComplex");
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
        writeToLayout("extremelycomplex", done[0], ns);

        Assert.That(ns.Count, Is.EqualTo(721));
        Assert.That(-Clipper.Area(ns), Is.EqualTo(Clipper.Area(done)));

        Console.WriteLine("  Decomposition (horizontal)....");
        sw.Restart();
        ns = GeoWrangler.rectangular_decomposition(ref abort, done, maxRayLength: rayLength, vertical: false);

        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Writing....");
        writeToLayout("extremelycomplex_horizontal", done[0], ns);

        Assert.That(ns.Count, Is.EqualTo(721));
        Assert.That(-Clipper.Area(ns), Is.EqualTo(Clipper.Area(done)));

        Console.WriteLine("  Done.");
    }

    [Test]
    public static void partFour_1()
    {
        Stopwatch sw = new();
        Console.WriteLine(" Part 1....");
        Console.WriteLine("  Preparing....");
        sw.Start();
        PathD points_1 = TestGeometry.ortho_fractal_1();
        points_1 = GeoWrangler.close(points_1);
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        partFour_do(points_1, "complex_loop");
    }

    [Test]
    public static void partFour_2()
    {
        Stopwatch sw = new();
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
        Stopwatch sw = new();

        bool vertical = true;
        bool abort = false;

        Console.WriteLine("  Keyhole....");
        // Give the keyholder a whirl:
        sw.Restart();
        PathD toDecomp = GeoWrangler.makeKeyHole(points, reverseEval: false, biDirectionalEval: true, customSizing: GeoWrangler.decomp_keyhole_sizing)[0];
        writeToLayout(baseString + "_kh", points, new PathsD { toDecomp });
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

        Assert.That(decompOut.Count, Is.EqualTo(481));
        Assert.That(Clipper.Area(decompOut), Is.EqualTo(147399));

        Console.WriteLine("  Decomposition (horizontal)....");
        sw.Restart();
        decompOut = GeoWrangler.rectangular_decomposition(ref abort, toDecomp,
            maxRayLength: (long)Math.Max(Math.Abs(dist.x), Math.Abs(dist.y)) * 1, vertical: !vertical);
        sw.Stop();
        Console.WriteLine("     done in " + sw.Elapsed.TotalSeconds + ".");

        Console.WriteLine("  Writing....");
        writeToLayout(baseString + "_horizontal", points, decompOut);

        Assert.That(decompOut.Count, Is.EqualTo(482));
        Assert.That(Clipper.Area(decompOut), Is.EqualTo(147399));
        Console.WriteLine("  Done.");

        sw.Stop();
    }

    [Test]
    public static void partFive()
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
            PathD points = new(polydata[i]);
            points = GeoWrangler.removeDuplicates(points);
            points = GeoWrangler.stripCollinear(points);
            points = GeoWrangler.clockwiseAndReorderXY(points);

            PathD toDecomp = GeoWrangler.makeKeyHole(GeoWrangler.sliverGapRemoval(points), reverseEval: false, biDirectionalEval: false)[0];
            PathD bounds = GeoWrangler.getBounds(toDecomp);
            PointD dist = GeoWrangler.distanceBetweenPoints_point(bounds[0], bounds[1]);

            PathsD decompOut = GeoWrangler.rectangular_decomposition(ref abort, toDecomp,
             maxRayLength: (long)Math.Max(Math.Abs(dist.x), Math.Abs(dist.y)), vertical: vertical);

            out_decomp.AddRange(decompOut);
        }

        Assert.That(out_decomp.Count, Is.EqualTo(19));
        Assert.That(Clipper.Area(out_decomp), Is.EqualTo(13843.0));
    }

    private static void writeToLayout(string filename, PathD orig, PathsD decomped)
    {
        // Can the system define geometry and write it correctly to Oasis and GDS files.
        GeoCore g = new();
        g.reset();
        int scale = 100; // for 0.01 nm resolution.

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
            databaseunits = 1000 * scale,
            userunits = 0.001 / scale,
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

        gds.gdsWriter gw = new(g, root_loc + filename + "_partitiontest.gds");
        gw.save();

        oasis.oasWriter ow = new(g, root_loc + filename + "_partitiontest.oas");
        ow.save();
    }
}