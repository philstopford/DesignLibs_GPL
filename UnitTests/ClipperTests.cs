using geoWrangler;

namespace UnitTests;

using PolyTree64 = Clipper2Lib.PolyTree64;
using Point64 = Clipper2Lib.Point64;
using PathD = Clipper2Lib.PathD;
using PathsD = Clipper2Lib.PathsD;
using Path64 = Clipper2Lib.Path64;
using Paths64 = Clipper2Lib.Paths64;
using IntPoint = ClipperLib1.IntPoint;
using Path = List<ClipperLib1.IntPoint>;
using Paths = List<List<ClipperLib1.IntPoint>>;
using SvgWriter = Clipper2Lib.SvgWriter;
using SvgUtils = Clipper2Lib.SvgUtils;

public class ClipperTests
{
    [SetUp]
    public static void ClipperTest()
    {
        // Clipper1 tests
        clipper1_collinearTest();
        clipper1_collinearOffsetTest();
        clipper1_unionTest();
        clipper1_notTest();
        clipper1_edgeOffsetTest();
        clipper1_zFillCallbackTest();
        clipper1_coincident_openPathTest();
        clipper1_keyHole_test1();
        clipper1_keyHole_test2();
        clipper1_openPath_clipTest1();
        clipper1_openPath_clipTest2();
        clipper1_offsetTest();
        clipper1_leftChordTest();
        clipper1_rightChordTest();

        // Clipper 2 tests
        clipper2_openpath_test();
        clipper2_openpath_complextest();
        clipper2_openpath_parallelLines();
        lineClipTest1();
        lineClipTest2();
        clipper2_H_subtract_test();
        clipper2_chordTest();
        clipper2_collinearTest();
        clipper2_collinearOffsetTest();
        clipper2_unionTest();
        clipper2_notTest();
        clipper2_edgeOffsetTest();
        clipper2_zFillCallbackTest();
        clipper2_coincident_openPathTest();
        clipper2_keyHole_test2();
        clipper2_openPath_clipTest1();
        clipper2_openPath_clipTest2();
        clipper2_keyHole_test1();
        clipper2_offsetTest();
        clipper2_leftChordTest();
        clipper2_rightChordTest();
        
        // Offset tests
        OpenPathOffsetTest1();
        OpenPathOffsetTest2();
        ClosedPathOffsetTest();

        clipper1_clipper2_compare_S();
        clipper1_clipper2_compare_X();
        clipper1_clipper2_compare_subtraction();
    }
    
    const double clipper1_keyhole_sizing = 500;

    [Test]
    public static void clipper1_collinearTest()
    {
        Path collinear = new()
        {
            new(-10, -10),
            new(-10, 0),
            new(-10, 10),
            new(0, 10),
            new(10, 10),
            new(10, 0),
            new(10, -10),
            new(0, -10),
            new(-10, -10),
        };

        ClipperLib1.Clipper c = new() {PreserveCollinear = true};
        c.AddPath(collinear, ClipperLib1.PolyType.ptSubject, true);
        c.AddPath(collinear, ClipperLib1.PolyType.ptClip, true);
        Paths output = new();
        c.Execute(ClipperLib1.ClipType.ctUnion, output);
        double area = output.Sum(t => ClipperLib1.Clipper.Area(t));
        Assert.AreEqual(400, area);
        Assert.AreEqual(1, output.Count);
        Assert.AreEqual(8, output[0].Count);
    }

    [Test]
    public static void clipper1_collinearOffsetTest()
    {
        Path collinear = new()
        {
            new(-10, -10),
            new(-10, 0),
            new(-10, 10),
            new(0, 10),
            new(10, 10),
            new(10, 0),
            new(10, -10),
            new(0, -10),
            new(-10, -10),
        };

        ClipperLib1.ClipperOffset co = new() {PreserveCollinear = true};
        co.AddPath(collinear, ClipperLib1.JoinType.jtMiter, ClipperLib1.EndType.etClosedPolygon);
        Paths temp = new();
        co.Execute(ref temp, 1.0);

        double area = temp.Sum(t => ClipperLib1.Clipper.Area(t));
        Assert.AreEqual(484, area);
        Assert.AreEqual(1, temp.Count);
        Assert.AreEqual(8, temp[0].Count);
    }
    
    [Test]
    public static void clipper1_unionTest()
    {
        Path simpleFirstPath = new()
        {
            new(100000, -300000),
            new(100000, 0),
            new(200000, 0),
            new(200000, -300000),
            new(100000, -300000),
        };
        
        Path simpleSecondPath = new()
        {
            new (0,-200000),
            new (0,-100000),
            new (300000,-100000),
            new (300000,-200000),
            new (0,-200000),
        };
        
        ClipperLib1.Clipper cs = new();
        
        cs.AddPath(simpleFirstPath, ClipperLib1.PolyType.ptSubject, true);
        cs.AddPath(simpleSecondPath, ClipperLib1.PolyType.ptClip, true);
        
        Paths simpleOutputPoints = new();
        cs.Execute(ClipperLib1.ClipType.ctUnion, simpleOutputPoints);

        double area1 = simpleOutputPoints.Sum(t => ClipperLib1.Clipper.Area(t));
        Assert.AreEqual(50000000000, area1);
        Assert.AreEqual(1, simpleOutputPoints.Count);
        Assert.AreEqual(12, simpleOutputPoints[0].Count);

        Path firstPath = new() {
        new (100000,-300000),
        new (100000,-290000),
        new (100000,-280000),
        new (100000,-270000),
        new (100000,-260000),
        new (100000,-250000),
        new (100000,-240000),
        new (100000,-230000),
        new (100000,-220000),
        new (100000,-210000),
        new (100000,-200000),
        new (100000,-190000),
        new (100000,-180000),
        new (100000,-170000),
        new (100000,-160000),
        new (100000,-150000),
        new (100000,-140000),
        new (100000,-130000),
        new (100000,-120000),
        new (100000,-110000),
        new (100000,-100000),
        new (100000,-90000),
        new (100000,-80000),
        new (100000,-70000),
        new (100000,-60000),
        new (100000,-50000),
        new (100000,-40000),
        new (100000,-30000),
        new (100000,-20000),
        new (100000,-10000),
        new (100000,0),
        new (110000,0),
        new (120000,0),
        new (130000,0),
        new (140000,0),
        new (150000,0),
        new (160000,0),
        new (170000,0),
        new (180000,0),
        new (190000,0),
        new (200000,0),
        new (200000,-10000),
        new (200000,-20000),
        new (200000,-30000),
        new (200000,-40000),
        new (200000,-50000),
        new (200000,-60000),
        new (200000,-70000),
        new (200000,-80000),
        new (200000,-90000),
        new (200000,-100000),
        new (200000,-110000),
        new (200000,-120000),
        new (200000,-130000),
        new (200000,-140000),
        new (200000,-150000),
        new (200000,-160000),
        new (200000,-170000),
        new (200000,-180000),
        new (200000,-190000),
        new (200000,-200000),
        new (200000,-210000),
        new (200000,-220000),
        new (200000,-230000),
        new (200000,-240000),
        new (200000,-250000),
        new (200000,-260000),
        new (200000,-270000),
        new (200000,-280000),
        new (200000,-290000),
        new (200000,-300000),
        new (190000,-300000),
        new (180000,-300000),
        new (170000,-300000),
        new (160000,-300000),
        new (150000,-300000),
        new (140000,-300000),
        new (130000,-300000),
        new (120000,-300000),
        new (110000,-300000),
        new (100000,-300000),
        };

        Path secondPath = new() {
        new (0,-200000),
        new (0,-190000),
        new (0,-180000),
        new (0,-170000),
        new (0,-160000),
        new (0,-150000),
        new (0,-140000),
        new (0,-130000),
        new (0,-120000),
        new (0,-110000),
        new (0,-100000),
        new (10000,-100000),
        new (20000,-100000),
        new (30000,-100000),
        new (40000,-100000),
        new (50000,-100000),
        new (60000,-100000),
        new (70000,-100000),
        new (80000,-100000),
        new (90000,-100000),
        new (100000,-100000),
        new (110000,-100000),
        new (120000,-100000),
        new (130000,-100000),
        new (140000,-100000),
        new (150000,-100000),
        new (160000,-100000),
        new (170000,-100000),
        new (180000,-100000),
        new (190000,-100000),
        new (200000,-100000),
        new (210000,-100000),
        new (220000,-100000),
        new (230000,-100000),
        new (240000,-100000),
        new (250000,-100000),
        new (260000,-100000),
        new (270000,-100000),
        new (280000,-100000),
        new (290000,-100000),
        new (300000,-100000),
        new (300000,-110000),
        new (300000,-120000),
        new (300000,-130000),
        new (300000,-140000),
        new (300000,-150000),
        new (300000,-160000),
        new (300000,-170000),
        new (300000,-180000),
        new (300000,-190000),
        new (300000,-200000),
        new (290000,-200000),
        new (280000,-200000),
        new (270000,-200000),
        new (260000,-200000),
        new (250000,-200000),
        new (240000,-200000),
        new (230000,-200000),
        new (220000,-200000),
        new (210000,-200000),
        new (200000,-200000),
        new (190000,-200000),
        new (180000,-200000),
        new (170000,-200000),
        new (160000,-200000),
        new (150000,-200000),
        new (140000,-200000),
        new (130000,-200000),
        new (120000,-200000),
        new (110000,-200000),
        new (100000,-200000),
        new (90000,-200000),
        new (80000,-200000),
        new (70000,-200000),
        new (60000,-200000),
        new (50000,-200000),
        new (40000,-200000),
        new (30000,-200000),
        new (20000,-200000),
        new (10000,-200000),
        new (0,-200000),
        };

        ClipperLib1.Clipper c = new();

        c.AddPath(firstPath, ClipperLib1.PolyType.ptSubject, true);
        c.AddPath(secondPath, ClipperLib1.PolyType.ptClip, true);

        Paths outputPoints = new();
        c.Execute(ClipperLib1.ClipType.ctUnion,  outputPoints);
        double area = outputPoints.Sum(t => ClipperLib1.Clipper.Area(t));
        Assert.AreEqual(50000000000, area);
        Assert.AreEqual(1, outputPoints.Count);
        Assert.AreEqual(12, outputPoints[0].Count);
    }

    [Test]
    public static void clipper1_notTest()
    {
        Path firstPath = new() {
        new (-250000,-250000),
        new (-250000,250000),
        new (250000,250000),
        new (250000,-250000),
        new (-250000,-250000),
        };
        
        Path secondPath = new() {
        new (-150000,-150000),
        new (-150000,150000),
        new (150000,150000),
        new (150000,-150000),
        new (-150000,-150000),
        new (-2147483647,-2147483647),
        new (-2147483647,2147483647),
        new (2147483647,2147483647),
        new (2147483647,-2147483647),
        new (-2147483647,-2147483647),
        new (-150000,-150000),
        };
        
        
        ClipperLib1.Clipper c = new();
        
        c.AddPath(firstPath, ClipperLib1.PolyType.ptSubject, true);
        c.AddPath(secondPath,ClipperLib1.PolyType.ptClip, true);
        
        Paths outputPoints = new();
        
        c.Execute(ClipperLib1.ClipType.ctIntersection, outputPoints);

        double area = outputPoints.Sum(t => ClipperLib1.Clipper.Area(t));
        Assert.AreEqual(160000000000, area);
        Assert.AreEqual(2, outputPoints.Count);
        Assert.AreEqual(4, outputPoints[0].Count);
        Assert.AreEqual(4, outputPoints[1].Count);
    }

    [Test]
    public static void clipper1_edgeOffsetTest()
    {
        Path edge = new()
        {
            new(-100000, 99500),
            new(-100000, 200500)
        };

        ClipperLib1.ClipperOffset co = new();
        co.AddPath(edge, ClipperLib1.JoinType.jtMiter, ClipperLib1.EndType.etOpenSquare);
        Paths p = new();
        co.Execute(ref p, 500);
        double area = p.Sum(t => ClipperLib1.Clipper.Area(t));
        Assert.AreEqual(102000000, area);
        Assert.AreEqual(1, p.Count);
        Assert.AreEqual(4, p[0].Count);
    }
    
    private static void clipper1_zFillTest(IntPoint bot1, IntPoint top1, IntPoint bot2, IntPoint top2, ref IntPoint pt)
    {
        pt.Z = -1;
    }

    [Test]
    public static void clipper1_zFillCallbackTest()
    {
        Path outer = new()
        {
            new(-1000, -1000),
            new(-1000, 1000),
            new(1000, 1000),
            new(1000, -1000)
        };
        
        Path cutter = new()
        {
            new(-100, -1100),
            new(-100, -900),
            new(100, -900),
            new(100, -1100),
        };

        ClipperLib1.Clipper c = new();

        c.ZFillFunction = clipper1_zFillTest;

        c.AddPath(outer, ClipperLib1.PolyType.ptSubject, true);
        c.AddPath(cutter, ClipperLib1.PolyType.ptClip, true);

        Paths solution = new();
        c.Execute(ClipperLib1.ClipType.ctIntersection, solution);
        double area = solution.Sum(t => ClipperLib1.Clipper.Area(t));
        Assert.AreEqual(20000, area);
        Assert.AreEqual(1, solution.Count);
        Assert.AreEqual(4, solution[0].Count);
    }
    
    [Test]
    public static void clipper1_coincident_openPathTest()
    {
        Path lPoly = new()
        {
            new(-1000, -1000),
            new(-1000, 1000),
            new(1000, 1000),
            new(1000, -1000)
        };

        Path t = new()
        {
            new (-1000, -1100),
            new (-1000, 500),
        };

        ClipperLib1.Clipper c = new()
        {
            PreserveCollinear = true,
            StrictlySimple = false
        };
        c.AddPath(lPoly, ClipperLib1.PolyType.ptClip, true);
        c.AddPath(t, ClipperLib1.PolyType.ptSubject, false);

        ClipperLib1.PolyTree pt = new();
        c.Execute(ClipperLib1.ClipType.ctIntersection, pt);
        Paths solution = ClipperLib1.Clipper.OpenPathsFromPolyTree(pt);
        
        Assert.AreEqual(1, solution.Count);
        Assert.AreEqual(2, solution[0].Count);

        Path t2 = new()
        {
            new (-1000, -1100),
            new (-1000, 500),
            new (-900, 500),
        };

        ClipperLib1.Clipper c2 = new()
        {
            PreserveCollinear = true,
            StrictlySimple = false
        };
        
        c2.AddPath(lPoly, ClipperLib1.PolyType.ptClip, true);
        c2.AddPath(t2, ClipperLib1.PolyType.ptSubject, false);

        ClipperLib1.PolyTree pt2 = new();
        c2.Execute(ClipperLib1.ClipType.ctIntersection, pt2);
        Paths solution2 = ClipperLib1.Clipper.OpenPathsFromPolyTree(pt2);

        Assert.AreEqual(1, solution2.Count);
        Assert.AreEqual(3, solution2[0].Count);

        Path t2b = new()
        {
            new (-900, 500),
            new (-1000, 500),
            new (-1000, -1100),
        };

        ClipperLib1.Clipper c2b = new()
        {
            PreserveCollinear = true,
            StrictlySimple = false
        };
        
        c2b.AddPath(lPoly, ClipperLib1.PolyType.ptClip, true);
        c2b.AddPath(t2b, ClipperLib1.PolyType.ptSubject, false);

        ClipperLib1.PolyTree pt2b = new();
        c2b.Execute(ClipperLib1.ClipType.ctIntersection, pt2b);
        Paths solution2b = ClipperLib1.Clipper.OpenPathsFromPolyTree(pt2b);

        Assert.AreEqual(1, solution2b.Count);
        Assert.AreEqual(3, solution2b[0].Count);

        Path t3 = new();
        int x = 0;
        int y = -1100;
        while (y < 1200)
        {
            t3.Add(new()
            {
                X = x,
                Y = y
            });
            y += 100;
        }
        ClipperLib1.Clipper c3 = new()
        {
            PreserveCollinear = true,
            StrictlySimple = false
        };
        c3.ZFillFunction = clipper1_zFillTest;
        c3.AddPath(lPoly, ClipperLib1.PolyType.ptClip, true);
        c3.AddPath(t3, ClipperLib1.PolyType.ptSubject, false);

        ClipperLib1.PolyTree pt3 = new();
        c3.Execute(ClipperLib1.ClipType.ctIntersection, pt3);
        Paths solution3 = ClipperLib1.Clipper.OpenPathsFromPolyTree(pt3);

        Assert.AreEqual(1, solution3.Count);
        Assert.AreEqual(21, solution3[0].Count);
    }

    [Test]
    public static void clipper1_keyHole_test2()
    {
        Path lPoly = new()
        {
            new IntPoint(200000, 0),
            new IntPoint(200000, 1100000),
            new IntPoint(1000000, 1100000),
            new IntPoint(1000000, 800000),
            new IntPoint(800000, 800000),
            new IntPoint(800000, 0),
            new IntPoint(200000, 0)
        };

        Paths t = new()
        {
            new Path()
            {
                new(800000, 800000),
                new(800000, 1100000)
            }
        };
        
        // Turn the new edges into cutters and slice. Not terribly elegant and we're relying on rounding to squash notches later.
        ClipperLib1.ClipperOffset co = new();
        co.AddPaths(t, ClipperLib1.JoinType.jtMiter, ClipperLib1.EndType.etOpenSquare);
        ClipperLib1.PolyTree tp = new();
        co.Execute(ref tp, 1.0);

        Paths cutters = ClipperLib1.Clipper.ClosedPathsFromPolyTree(tp);

        double area = cutters.Sum(t => ClipperLib1.Clipper.Area(t));
        Assert.AreEqual(600004, area);
        Assert.AreEqual(1, cutters.Count);
        Assert.AreEqual(4, cutters[0].Count);

        ClipperLib1.Clipper c = new();

        c.AddPath(lPoly, ClipperLib1.PolyType.ptSubject, true);

        // Take first cutter only - we only cut once, no matter how many potential cutters we have.
        c.AddPath(cutters[0], ClipperLib1.PolyType.ptClip, true);
        Paths f = new();
        c.Execute(ClipperLib1.ClipType.ctDifference, f, ClipperLib1.PolyFillType.pftEvenOdd, ClipperLib1.PolyFillType.pftEvenOdd);
        double area2 = f.Sum(t => ClipperLib1.Clipper.Area(t));
        Assert.AreEqual(719999399999, area2);
        Assert.AreEqual(2, f.Count);
        Assert.AreEqual(6, f[0].Count);
        Assert.AreEqual(4, f[1].Count);
    }
    
    [Test]
    public static void clipper1_openPath_clipTest1()
    {
        Path lPoly = new() {
        new IntPoint(0,0),
        new IntPoint(0,200000),
        new IntPoint(200000,200000),
        new IntPoint(200000,500000),
        new IntPoint(0,500000),
        new IntPoint(0,1100000),
        new IntPoint(1000000,1100000),
        new IntPoint(1000000,800000),
        new IntPoint(800000,800000),
        new IntPoint(800000,600000),
        new IntPoint(1000000,600000),
        new IntPoint(1000000,0),
        new IntPoint(0,0)
        };

        Path t = new () {
        new IntPoint(0,200000),
        new IntPoint(0,-9800000)
        };

        ClipperLib1.Clipper c = new();

        c.AddPath(t, ClipperLib1.PolyType.ptSubject, false);
        c.AddPath(lPoly, ClipperLib1.PolyType.ptClip, true);

        ClipperLib1.PolyTree pt = new();

        c.Execute(ClipperLib1.ClipType.ctIntersection, pt);
        Paths solution = ClipperLib1.Clipper.OpenPathsFromPolyTree(pt);
 
        Assert.AreEqual(1, solution.Count);
    }

    [Test]
    public static void clipper1_keyHole_test1()
    {
        Console.WriteLine("Clipper1 Test1");
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
        Paths kHSource = new()
        {
            outer,
            inner1
        };
        
        ClipperLib1.ClipperOffset co = new();
        co.AddPaths(kHSource, ClipperLib1.JoinType.jtMiter, ClipperLib1.EndType.etClosedPolygon);
        Paths out_ = new ();
        co.Execute(ref out_, clipper1_keyhole_sizing);

        double area = out_.Sum(t => ClipperLib1.Clipper.Area(t));
        Assert.AreEqual(121200000000, area);
        Assert.AreEqual(2, out_.Count);
        Assert.AreEqual(4, out_[0].Count);
        Assert.AreEqual(4, out_[1].Count);
    }

    [Test]
    public static void clipper1_openPath_clipTest2()
    {
        Paths rays = new()
        {
            new Path() {new(100000, 200000), new(100000, -9800000)}
        };

        Paths collisionPaths = new()
        {
            new Path()
            {
                new(0, 0),
                new(0, 500000),
                new(100000, 500000),
                new(100000, 200000),
                new(600000, 200000),
                new(600000, 800000),
                new(1200000, 800000),
                new(1200000, 0),
                new(0, 0)

            }
        };

        ClipperLib1.Clipper c = new ClipperLib1.Clipper();
        c.AddPaths(rays, ClipperLib1.PolyType.ptSubject, false);
        c.AddPaths(collisionPaths, ClipperLib1.PolyType.ptClip, true);
        ClipperLib1.PolyTree pt = new ClipperLib1.PolyTree();
        c.Execute(ClipperLib1.ClipType.ctIntersection, pt);
        Paths solution = ClipperLib1.Clipper.OpenPathsFromPolyTree(pt);
        Assert.AreEqual(1, solution.Count);
        Assert.AreEqual(2, solution[0].Count);
    }

    [Test]
    public static void clipper1_offsetTest()
    {
        Path lPoly = new ()
        {
            new IntPoint(0, 0),
            new IntPoint(0, 500000),
            new IntPoint(100000, 500000),
            new IntPoint(100000, 200000),
            new IntPoint(600000, 200000),
            new IntPoint(600000, 0),
            new IntPoint(0, 0)
        };

        Path newEdge = new()
        {
            new IntPoint(100000, 200000),
            new IntPoint(100000, 0)
        };

        Paths newEdges = new()
        {
            newEdge
        };

        ClipperLib1.ClipperOffset co = new ClipperLib1.ClipperOffset();
        co.AddPaths(newEdges, ClipperLib1.JoinType.jtMiter, ClipperLib1.EndType.etOpenSquare);
        ClipperLib1.PolyTree tp = new();
        co.Execute(ref tp, 1.0);

        Paths cutters = ClipperLib1.Clipper.ClosedPathsFromPolyTree(tp);

        ClipperLib1.Clipper c = new ClipperLib1.Clipper();
        c.AddPath(lPoly, ClipperLib1.PolyType.ptSubject, true);
        c.AddPaths(cutters, ClipperLib1.PolyType.ptClip, true);
        Paths solution = new();
        c.Execute(ClipperLib1.ClipType.ctDifference, solution);

        double area = solution.Sum(t => ClipperLib1.Clipper.Area(t));
        Assert.AreEqual(149999599999, area);
        Assert.AreEqual(2, solution.Count);
        Assert.AreEqual(6, solution[0].Count);
        Assert.AreEqual(4, solution[1].Count);
    }
    
    [Test]
    public static void clipper1_leftChordTest()
    {
        Path testPath = new List<IntPoint>() {
         new(-200000,0),
         new(-300000,0),
         new(-310453,-548),
         new(-320791,-2185),
         new(-330902,-4894),
         new(-340674,-8645),
         new(-350000,-13397),
         new(-358779,-19098),
         new(-366913,-25686),
         new(-374314,-33087),
         new(-380902,-41221),
         new(-386603,-50000),
         new(-391355,-59326),
         new(-395106,-69098),
         new(-397815,-79209),
         new(-399452,-89547),
         new(-400000,-100000),
         new(-400000,-700000),
         new(-400000,-700000),
         new(-398907,-710396),
         new(-397815,-720791),
         new(-395106,-730902),
         new(-391355,-740674),
         new(-386603,-750000),
         new(-380902,-758779),
         new(-374314,-766913),
         new(-366913,-774314),
         new(-358779,-780902),
         new(-350000,-786603),
         new(-340674,-791355),
         new(-330902,-795106),
         new(-320791,-797815),
         new(-310453,-799452),
         new(-300000,-800000),
         };

        Path a = new List<IntPoint>() {
          new(-460000,-390000),
          new(-459955,-395235),
          new(-459594,-405679),
          new(-458873,-416047),
          new(-457796,-426288),
          new(-456369,-436353),
          new(-454102,-448610),
          new(-451315,-460421),
          new(-448030,-471696),
          new(-444272,-482349),
          new(-440068,-492300),
          new(-435451,-501472),
          new(-429416,-511353),
          new(-422903,-519904),
          new(-414796,-528076),
          new(-406257,-534189),
          new(-396132,-538540),
          new(-385806,-540000),
          new(-381067,-540036),
          new(-369255,-540441),
          new(-357570,-541292),
          new(-346100,-542583),
          new(-334932,-544305),
          new(-324151,-546443),
          new(-313840,-548983),
          new(-304076,-551904),
          new(-293187,-555882),
          new(-283312,-560333),
          new(-273218,-566059),
          new(-264802,-572279),
          new(-257399,-579871),
          new(-252063,-588852),
          new(-250000,-599117),
          new(-242658,-614841),
          new(-234819,-621423),
          new(-225801,-626830),
          new(-216572,-631139),
          new(-206066,-635097),
          new(-196418,-638096),
          new(-186036,-640798),
          new(-175000,-643183),
          new(-163393,-645233),
          new(-151303,-646931),
          new(-141346,-648029),
          new(-131187,-648888),
          new(-120876,-649505),
          new(-110463,-649876),
          new(-100000,-650000),
          new(-94765,-649967),
          new(-84321,-649701),
          new(-73953,-649169),
          new(-63712,-648376),
          new(-53647,-647324),
          new(-41390,-645654),
          new(-29579,-643601),
          new(-18304,-641180),
          new(-7651,-638411),
          new(2300,-635314),
          new(13206,-631197),
          new(22873,-626688),
          new(32442,-620997),
          new(40954,-614029),
          new(47244,-605763),
          new(50000,-595332),
          new(50438,-591472),
          new(55347,-581946),
          new(63107,-574604),
          new(72568,-568506),
          new(82553,-563595),
          new(92112,-559765),
          new(102721,-556206),
          new(114298,-552945),
          new(124199,-550567),
          new(134615,-548408),
          new(145495,-546477),
          new(156787,-544784),
          new(168436,-543337),
          new(180385,-542143),
          new(192576,-541209),
          new(204949,-540538),
          new(217444,-540135),
          new(230000,-540000),
          new(237542,-540000),
          new(241117,-539909),
          new(251801,-538540),
          new(262329,-535544),
          new(272585,-530954),
          new(282456,-524819),
          new(290312,-518575),
          new(297765,-511353),
          new(304760,-503206),
          new(311244,-494199),
          new(317167,-484398),
          new(322483,-473879),
          new(327154,-462721),
          new(331142,-451010),
          new(333821,-441303),
          new(336031,-431346),
          new(337761,-421187),
          new(339003,-410876),
          new(339750,-400463),
          new(340000,-390000),
          new(339909,-384765),
          new(339178,-374321),
          new(337721,-363953),
          new(335544,-353712),
          new(332658,-343647),
          new(329078,-333809),
          new(324819,-324244),
          new(319904,-315000),
          new(314356,-306121),
          new(308202,-297651),
          new(301472,-289630),
          new(294199,-282099),
          new(286418,-275093),
          new(278168,-268647),
          new(269488,-262793),
          new(260421,-257558),
          new(251010,-252968),
          new(241303,-249046),
          new(231346,-245811),
          new(221187,-243278),
          new(210876,-241460),
          new(200463,-240365),
          new(190000,-240000),
          new(189927,-240000),
          new(180166,-239562),
          new(168038,-237784),
          new(156076,-234653),
          new(146687,-231190),
          new(137509,-226893),
          new(128587,-221783),
          new(119964,-215885),
          new(111681,-209227),
          new(103779,-201842),
          new(96298,-193766),
          new(89272,-185039),
          new(82737,-175702),
          new(76724,-165801),
          new(71262,-155385),
          new(66379,-144505),
          new(62097,-133213),
          new(58439,-121564),
          new(55421,-109615),
          new(53058,-97424),
          new(51362,-85051),
          new(50341,-72556),
          new(50000,-60000),
          new(50000,-30000),
          new(49909,-24765),
          new(49178,-14321),
          new(47721,-3953),
          new(45544,6288),
          new(42658,16353),
          new(39078,26191),
          new(34819,35756),
          new(29904,45000),
          new(24356,53879),
          new(18202,62349),
          new(11472,70370),
          new(4199,77901),
          new(-3582,84907),
          new(-11832,91353),
          new(-20512,97207),
          new(-29579,102442),
          new(-38990,107032),
          new(-48697,110954),
          new(-58654,114189),
          new(-68813,116722),
          new(-79124,118540),
          new(-89537,119635),
          new(-100000,120000),
          new(-105235,119909),
          new(-115679,119178),
          new(-126047,117721),
          new(-136288,115544),
          new(-146353,112658),
          new(-156191,109078),
          new(-165756,104819),
          new(-175000,99904),
          new(-183879,94356),
          new(-192349,88202),
          new(-200370,81472),
          new(-207901,74199),
          new(-214907,66418),
          new(-221353,58168),
          new(-227207,49488),
          new(-232442,40421),
          new(-237032,31010),
          new(-240954,21303),
          new(-244189,11346),
          new(-246722,1187),
          new(-248540,-9124),
          new(-249635,-19537),
          new(-250000,-30000),
          new(-250000,-60000),
          new(-250062,-66282),
          new(-250555,-78815),
          new(-251539,-91257),
          new(-253010,-103546),
          new(-254959,-115623),
          new(-257378,-127429),
          new(-260255,-138907),
          new(-263575,-150000),
          new(-267323,-160655),
          new(-271480,-170819),
          new(-276026,-180444),
          new(-280939,-189481),
          new(-287560,-199886),
          new(-294665,-209227),
          new(-302202,-217432),
          new(-310113,-224438),
          new(-320015,-231190),
          new(-330260,-236067),
          new(-340735,-239014),
          new(-351326,-240000),
          new(-357014,-240206),
          new(-368327,-241847),
          new(-379453,-245111),
          new(-390272,-249963),
          new(-398966,-255181),
          new(-407297,-261425),
          new(-415203,-268647),
          new(-422623,-276794),
          new(-429500,-285801),
          new(-435782,-295602),
          new(-441421,-306121),
          new(-446374,-317279),
          new(-450605,-328990),
          new(-453446,-338697),
          new(-455790,-348654),
          new(-457625,-358813),
          new(-458942,-369124),
          new(-459735,-379537),
          new(-460000,-390000),
          };

        // Let's see if we can track the origin of the chords.
        for (int p = 0; p < a.Count; p++)
        {
            a[p] = new()
            {
                X = a[p].X,
                Y = a[p].Y,
                Z = 1
            };
        }

        for (int p = 0; p < testPath.Count; p++)
        {
            testPath[p] = new()
            {
                X = testPath[p].X,
                Y = testPath[p].Y,
                Z = 2
            };
        }
        
        ClipperLib1.Clipper c = new();
        c.ZFillFunction = clipper1_zFillTest;
        c.AddPath(testPath, ClipperLib1.PolyType.ptSubject, false);
        c.AddPath(a, ClipperLib1.PolyType.ptClip, true);

        ClipperLib1.PolyTree pt = new();
        c.Execute(ClipperLib1.ClipType.ctIntersection, pt);
        Paths open = ClipperLib1.Clipper.OpenPathsFromPolyTree(pt);
        
        Assert.AreEqual(2, open.Count);
        Assert.AreEqual(2, open[0].Count);
        Assert.AreEqual(2, open[1].Count);
    }
    
    [Test]
    public static void clipper1_rightChordTest()
    {
        Path testPath = new ()
        {
            new(-400000, -700000),
            new(-398907, -710396),
            new(-397815, -720791),
            new(-395106, -730902),
            new(-391355, -740674),
            new(-386603, -750000),
            new(-380902, -758779),
            new(-374314, -766913),
            new(-366913, -774314),
            new(-358779, -780902),
            new(-350000, -786603),
            new(-340674, -791355),
            new(-330902, -795106),
            new(-320791, -797815),
            new(-310453, -799452),
            new(-300000, -800000),
            new(-200000, -800000),
            new(-189547, -799452),
            new(-179209, -797815),
            new(-169098, -795106),
            new(-159326, -791355),
            new(-150000, -786603),
            new(-141221, -780902),
            new(-133087, -774314),
            new(-125686, -766913),
            new(-119098, -758779),
            new(-113397, -750000),
            new(-108645, -740674),
            new(-104894, -730902),
            new(-102185, -720791),
            new(-100548, -710453),
            new(-100000, -700000),
            new(-100000, -100000),
            new(-100548, -89547),
            new(-102185, -79209),
            new(-104894, -69098),
            new(-108645, -59326),
            new(-113397, -50000),
            new(-119098, -41221),
            new(-125686, -33087),
            new(-133087, -25686),
            new(-141221, -19098),
            new(-150000, -13397),
            new(-159326, -8645),
            new(-169098, -4894),
            new(-179209, -2185),
            new(-189547, -548),
            new(-200000, 0)
        };

        Path a = new ()
        {
            new(-460000, -390000),
            new(-459955, -395235),
            new(-459594, -405679),
            new(-458873, -416047),
            new(-457796, -426288),
            new(-456369, -436353),
            new(-454102, -448610),
            new(-451315, -460421),
            new(-448030, -471696),
            new(-444272, -482349),
            new(-440068, -492300),
            new(-435451, -501472),
            new(-429416, -511353),
            new(-422903, -519904),
            new(-414796, -528076),
            new(-406257, -534189),
            new(-396132, -538540),
            new(-385806, -540000),
            new(-381067, -540036),
            new(-369255, -540441),
            new(-357570, -541292),
            new(-346100, -542583),
            new(-334932, -544305),
            new(-324151, -546443),
            new(-313840, -548983),
            new(-304076, -551904),
            new(-293187, -555882),
            new(-283312, -560333),
            new(-273218, -566059),
            new(-264802, -572279),
            new(-257399, -579871),
            new(-252063, -588852),
            new(-250000, -599117),
            new(-242658, -614841),
            new(-234819, -621423),
            new(-225801, -626830),
            new(-216572, -631139),
            new(-206066, -635097),
            new(-196418, -638096),
            new(-186036, -640798),
            new(-175000, -643183),
            new(-163393, -645233),
            new(-151303, -646931),
            new(-141346, -648029),
            new(-131187, -648888),
            new(-120876, -649505),
            new(-110463, -649876),
            new(-100000, -650000),
            new(-94765, -649967),
            new(-84321, -649701),
            new(-73953, -649169),
            new(-63712, -648376),
            new(-53647, -647324),
            new(-41390, -645654),
            new(-29579, -643601),
            new(-18304, -641180),
            new(-7651, -638411),
            new(2300, -635314),
            new(13206, -631197),
            new(22873, -626688),
            new(32442, -620997),
            new(40954, -614029),
            new(47244, -605763),
            new(50000, -595332),
            new(50438, -591472),
            new(55347, -581946),
            new(63107, -574604),
            new(72568, -568506),
            new(82553, -563595),
            new(92112, -559765),
            new(102721, -556206),
            new(114298, -552945),
            new(124199, -550567),
            new(134615, -548408),
            new(145495, -546477),
            new(156787, -544784),
            new(168436, -543337),
            new(180385, -542143),
            new(192576, -541209),
            new(204949, -540538),
            new(217444, -540135),
            new(230000, -540000),
            new(237542, -540000),
            new(241117, -539909),
            new(251801, -538540),
            new(262329, -535544),
            new(272585, -530954),
            new(282456, -524819),
            new(290312, -518575),
            new(297765, -511353),
            new(304760, -503206),
            new(311244, -494199),
            new(317167, -484398),
            new(322483, -473879),
            new(327154, -462721),
            new(331142, -451010),
            new(333821, -441303),
            new(336031, -431346),
            new(337761, -421187),
            new(339003, -410876),
            new(339750, -400463),
            new(340000, -390000),
            new(339909, -384765),
            new(339178, -374321),
            new(337721, -363953),
            new(335544, -353712),
            new(332658, -343647),
            new(329078, -333809),
            new(324819, -324244),
            new(319904, -315000),
            new(314356, -306121),
            new(308202, -297651),
            new(301472, -289630),
            new(294199, -282099),
            new(286418, -275093),
            new(278168, -268647),
            new(269488, -262793),
            new(260421, -257558),
            new(251010, -252968),
            new(241303, -249046),
            new(231346, -245811),
            new(221187, -243278),
            new(210876, -241460),
            new(200463, -240365),
            new(190000, -240000),
            new(189927, -240000),
            new(180166, -239562),
            new(168038, -237784),
            new(156076, -234653),
            new(146687, -231190),
            new(137509, -226893),
            new(128587, -221783),
            new(119964, -215885),
            new(111681, -209227),
            new(103779, -201842),
            new(96298, -193766),
            new(89272, -185039),
            new(82737, -175702),
            new(76724, -165801),
            new(71262, -155385),
            new(66379, -144505),
            new(62097, -133213),
            new(58439, -121564),
            new(55421, -109615),
            new(53058, -97424),
            new(51362, -85051),
            new(50341, -72556),
            new(50000, -60000),
            new(50000, -30000),
            new(49909, -24765),
            new(49178, -14321),
            new(47721, -3953),
            new(45544, 6288),
            new(42658, 16353),
            new(39078, 26191),
            new(34819, 35756),
            new(29904, 45000),
            new(24356, 53879),
            new(18202, 62349),
            new(11472, 70370),
            new(4199, 77901),
            new(-3582, 84907),
            new(-11832, 91353),
            new(-20512, 97207),
            new(-29579, 102442),
            new(-38990, 107032),
            new(-48697, 110954),
            new(-58654, 114189),
            new(-68813, 116722),
            new(-79124, 118540),
            new(-89537, 119635),
            new(-100000, 120000),
            new(-105235, 119909),
            new(-115679, 119178),
            new(-126047, 117721),
            new(-136288, 115544),
            new(-146353, 112658),
            new(-156191, 109078),
            new(-165756, 104819),
            new(-175000, 99904),
            new(-183879, 94356),
            new(-192349, 88202),
            new(-200370, 81472),
            new(-207901, 74199),
            new(-214907, 66418),
            new(-221353, 58168),
            new(-227207, 49488),
            new(-232442, 40421),
            new(-237032, 31010),
            new(-240954, 21303),
            new(-244189, 11346),
            new(-246722, 1187),
            new(-248540, -9124),
            new(-249635, -19537),
            new(-250000, -30000),
            new(-250000, -60000),
            new(-250062, -66282),
            new(-250555, -78815),
            new(-251539, -91257),
            new(-253010, -103546),
            new(-254959, -115623),
            new(-257378, -127429),
            new(-260255, -138907),
            new(-263575, -150000),
            new(-267323, -160655),
            new(-271480, -170819),
            new(-276026, -180444),
            new(-280939, -189481),
            new(-287560, -199886),
            new(-294665, -209227),
            new(-302202, -217432),
            new(-310113, -224438),
            new(-320015, -231190),
            new(-330260, -236067),
            new(-340735, -239014),
            new(-351326, -240000),
            new(-357014, -240206),
            new(-368327, -241847),
            new(-379453, -245111),
            new(-390272, -249963),
            new(-398966, -255181),
            new(-407297, -261425),
            new(-415203, -268647),
            new(-422623, -276794),
            new(-429500, -285801),
            new(-435782, -295602),
            new(-441421, -306121),
            new(-446374, -317279),
            new(-450605, -328990),
            new(-453446, -338697),
            new(-455790, -348654),
            new(-457625, -358813),
            new(-458942, -369124),
            new(-459735, -379537),
            new(-460000, -390000),
        };

        // Let's see if we can track the origin of the chords.
        for (int p = 0; p < a.Count; p++)
        {
            a[p] = new()
            {
                X = a[p].X,
                Y = a[p].Y,
                Z = 1
            };
        }

        for (int p = 0; p < testPath.Count; p++)
        {
            testPath[p] = new()
            {
                X = testPath[p].X,
                Y = testPath[p].Y,
                Z = 2
            };
        }

        ClipperLib1.Clipper c = new();
        c.ZFillFunction = clipper1_zFillTest;
        c.AddPath(testPath, ClipperLib1.PolyType.ptSubject, false);
        c.AddPath(a, ClipperLib1.PolyType.ptClip, true);

        ClipperLib1.PolyTree pt = new();
        c.Execute(ClipperLib1.ClipType.ctIntersection, pt);
        Paths open = ClipperLib1.Clipper.OpenPathsFromPolyTree(pt);
        Assert.AreEqual(1, open.Count);
        Assert.AreEqual(17, open[0].Count);
    }
    
    private static readonly Paths64 lineClip_collisionGeometry = new () {
        new () {
            new (0,-250000),
            new (0,-230000),
            new (0,-210000),
            new (0,-190000),
            new (0,-170000),
            new (0,-150000),
            new (0,-130000),
            new (0,-110000),
            new (0,-90000),
            new (0,-70000),
            new (0,-50000),
            new (4323,-29663),
            new (16543,-12843),
            new (34549,-2447),
            new (50000,0),
            new (70337,-4323),
            new (87157,-16543),
            new (97553,-34549),
            new (100000,-50000),
            new (100000,-70000),
            new (100000,-90000),
            new (100000,-110000),
            new (100000,-130000),
            new (100000,-150000),
            new (100000,-170000),
            new (100000,-190000),
            new (100000,-210000),
            new (100000,-230000),
            new (100000,-250000),
            new (95677,-270337),
            new (83457,-287157),
            new (65451,-297553),
            new (50000,-300000),
            new (29663,-295677),
            new (12843,-283457),
            new (0,-250000),
        },
        new () {
            new (200000,-250000),
            new (200000,-230000),
            new (200000,-210000),
            new (200000,-190000),
            new (200000,-170000),
            new (200000,-150000),
            new (200000,-130000),
            new (200000,-110000),
            new (200000,-90000),
            new (200000,-70000),
            new (200000,-50000),
            new (204323,-29663),
            new (216543,-12843),
            new (234549,-2447),
            new (250000,0),
            new (270337,-4323),
            new (287157,-16543),
            new (297553,-34549),
            new (300000,-50000),
            new (300000,-70000),
            new (300000,-90000),
            new (300000,-110000),
            new (300000,-130000),
            new (300000,-150000),
            new (300000,-170000),
            new (300000,-190000),
            new (300000,-210000),
            new (300000,-230000),
            new (300000,-250000),
            new (295677,-270337),
            new (283457,-287157),
            new (265451,-297553),
            new (250000,-300000),
            new (229663,-295677),
            new (212843,-283457),
            new (200000,-250000),
        }
    };

    [Test]
    public static void lineClipTest1()
    {
        Paths64 castLines = new() {
        new() {
        new (0,-250000),
        new (-291700,-320076),
        },
        new() {
        new (0,-230000),
        new (-300000,-230000),
        },
        new() {
        new (0,-210000),
        new (-300000,-210000),
        },
        new() {
        new (0,-190000),
        new (-300000,-190000),
        },
        new() {
        new (0,-170000),
        new (-300000,-170000),
        },
        new() {
        new (0,-150000),
        new (-300000,-150000),
        },
        new() {
        new (0,-130000),
        new (-300000,-130000),
        },
        new() {
        new (0,-110000),
        new (-300000,-110000),
        },
        new() {
        new (0,-90000),
        new (-300000,-90000),
        },
        new() {
        new (0,-70000),
        new (-300000,-70000),
        },
        new() {
        new (0,-50000),
        new (-298292,-18037),
        },
        new() {
        new (4323,-29663),
        new (-269743,92353),
        },
        new() {
        new (16543,-12843),
        new (-184197,210100),
        },
        new() {
        new (34549,-2447),
        new (-72957,277629),
        },
        new() {
        new (50000,0),
        new (65704,299588),
        },
        new() {
        new (70337,-4323),
        new (192352,269743),
        },
        new() {
        new (87157,-16543),
        new (310099,184197),
        },
        new() {
        new (97553,-34549),
        new (377628,72957),
        },
        new() {
        new (100000,-50000),
        new (399288,-29350),
        },
        new() {
        new (100000,-70000),
        new (400000,-70000),
        },
        new() {
        new (100000,-90000),
        new (400000,-90000),
        },
        new() {
        new (100000,-110000),
        new (400000,-110000),
        },
        new() {
        new (100000,-130000),
        new (400000,-130000),
        },
        new() {
        new (100000,-150000),
        new (400000,-150000),
        },
        new() {
        new (100000,-170000),
        new (400000,-170000),
        },
        new() {
        new (100000,-190000),
        new (400000,-190000),
        },
        new() {
        new (100000,-210000),
        new (400000,-210000),
        },
        new() {
        new (100000,-230000),
        new (400000,-230000),
        },
        new() {
        new (100000,-250000),
        new (398292,-281963),
        },
        new() {
        new (95677,-270337),
        new (369743,-392353),
        },
        new() {
        new (83457,-287157),
        new (284197,-510100),
        },
        new() {
        new (65451,-297553),
        new (172957,-577629),
        },
        new() {
        new (50000,-300000),
        new (34296,-599588),
        },
        new() {
        new (29663,-295677),
        new (-92352,-569743),
        },
        new() {
        new (12843,-283457),
        new (-238759,-446847),
        },
        new() {
        new (0,-250000),
        new (-291700,-320076),}
        };

        // Create our Clipper1 collision geometry from the Clipper2 definition
        Paths collsionGeometry1 = new();
        foreach (List<Clipper2Lib.Point64> t in lineClip_collisionGeometry)
        {
            List<ClipperLib1.IntPoint> a = new List<ClipperLib1.IntPoint>();
            for (int j = 0; j < t.Count; j++)
            {
                a.Add(new ClipperLib1.IntPoint(t[j].X, t[j].Y));
            }
            collsionGeometry1.Add(a);
        }

        foreach (Path64 t in castLines)
        {
            Clipper2Lib.Clipper64 c2 = new();
            c2.AddOpenSubject(t);
            c2.AddClip(lineClip_collisionGeometry);
            Paths64 unused2 = new();
            Paths64 clipped2 = new();
            c2.Execute(Clipper2Lib.ClipType.Difference, Clipper2Lib.FillRule.EvenOdd, unused2, clipped2);

            ClipperLib1.Clipper c1 = new();
            List<ClipperLib1.IntPoint> c1Ray = new()
                {
                    new (t[0].X, t[0].Y),
                    new (t[1].X, t[1].Y),
                }
                ;
            c1.AddPath(c1Ray, ClipperLib1.PolyType.ptSubject, false);
            c1.AddPaths(collsionGeometry1, ClipperLib1.PolyType.ptClip, true);
            ClipperLib1.PolyTree pt = new();
            c1.Execute(ClipperLib1.ClipType.ctDifference, pt);
            Paths clipped1 = ClipperLib1.Clipper.OpenPathsFromPolyTree(pt);
            
            // Results ought to be ~consistent. ClipperLib2 reverses the direction compared to ClipperLib1
            Assert.True(((Math.Abs(clipped1[0][0].X - clipped2[0][0].X) <= 1) && (Math.Abs(clipped1[0][0].Y - clipped2[0][0].Y) <= 1)) ||
                        ((Math.Abs(clipped1[0][0].X - clipped2[0][1].X) <= 1) && (Math.Abs(clipped1[0][0].Y - clipped2[0][1].Y) <= 1)));
        }
    }

    [Test]
    public static void lineClipTest2()
    {
        Paths64 castLines = new() {
            new() {
                new (200000,-250000),
                new (-91700,-320076),
            },
            new() {
                new (200000,-230000),
                new (-100000,-230000),
            },
            new() {
                new (200000,-210000),
                new (-100000,-210000),
            },
            new() {
                new (200000,-190000),
                new (-100000,-190000),
            },
            new() {
                new (200000,-170000),
                new (-100000,-170000),
            },
            new() {
                new (200000,-150000),
                new (-100000,-150000),
            },
            new() {
                new (200000,-130000),
                new (-100000,-130000),
            },
            new() {
                new (200000,-110000),
                new (-100000,-110000),
            },
            new() {
                new (200000,-90000),
                new (-100000,-90000),
            },
            new() {
                new (200000,-70000),
                new (-100000,-70000),
            },
            new() {
                new (200000,-50000),
                new (-98292,-18037),
            },
            new() {
                new (204323,-29663),
                new (-69743,92353),
            },
            new() {
                new (216543,-12843),
                new (15803,210100),
            },
            new() {
                new (234549,-2447),
                new (127043,277629),
            },
            new() {
                new (250000,0),
                new (265704,299588),
            },
            new() {
                new (270337,-4323),
                new (392352,269743),
            },
            new() {
                new (287157,-16543),
                new (510099,184197),
            },
            new() {
                new (297553,-34549),
                new (577628,72957),
            },
            new() {
                new (300000,-50000),
                new (599288,-29350),
            },
            new() {
                new (300000,-70000),
                new (600000,-70000),
            },
            new() {
                new (300000,-90000),
                new (600000,-90000),
            },
            new() {
                new (300000,-110000),
                new (600000,-110000),
            },
            new() {
                new (300000,-130000),
                new (600000,-130000),
            },
            new() {
                new (300000,-150000),
                new (600000,-150000),
            },
            new() {
                new (300000,-170000),
                new (600000,-170000),
            },
            new() {
                new (300000,-190000),
                new (600000,-190000),
            },
            new() {
                new (300000,-210000),
                new (600000,-210000),
            },
            new() {
                new (300000,-230000),
                new (600000,-230000),
            },
            new() {
                new (300000,-250000),
                new (598292,-281963),
            },
            new() {
                new (295677,-270337),
                new (569743,-392353),
            },
            new() {
                new (283457,-287157),
                new (484197,-510100),
            },
            new() {
                new (265451,-297553),
                new (372957,-577629),
            },
            new() {
                new (250000,-300000),
                new (234296,-599588),
            },
            new() {
                new (229663,-295677),
                new (107648,-569743),
            },
            new() {
                new (212843,-283457),
                new (-38759,-446847),
            },
            new() {
                new (200000,-250000),
                new (-91700,-320076)
            }
        };        
        
        // Create our Clipper1 collision geometry from the Clipper2 definition
        Paths collsionGeometry1 = new();
        foreach (List<Clipper2Lib.Point64> t in lineClip_collisionGeometry)
        {
            List<ClipperLib1.IntPoint> a = new List<ClipperLib1.IntPoint>();
            for (int j = 0; j < t.Count; j++)
            {
                a.Add(new ClipperLib1.IntPoint(t[j].X, t[j].Y));
            }
            collsionGeometry1.Add(a);
        }

        foreach (Path64 t in castLines)
        {
            Clipper2Lib.Clipper64 c2 = new();
            c2.AddOpenSubject(t);
            c2.AddClip(lineClip_collisionGeometry);
            Paths64 unused2 = new();
            Paths64 clipped2 = new();
            c2.Execute(Clipper2Lib.ClipType.Difference, Clipper2Lib.FillRule.EvenOdd, unused2, clipped2);

            ClipperLib1.Clipper c1 = new();
            List<ClipperLib1.IntPoint> c1Ray = new()
                {
                    new (t[0].X, t[0].Y),
                    new (t[1].X, t[1].Y),
                }
                ;
            c1.AddPath(c1Ray, ClipperLib1.PolyType.ptSubject, false);
            c1.AddPaths(collsionGeometry1, ClipperLib1.PolyType.ptClip, true);
            ClipperLib1.PolyTree pt = new();
            c1.Execute(ClipperLib1.ClipType.ctDifference, pt);
            Paths clipped1 = ClipperLib1.Clipper.OpenPathsFromPolyTree(pt);
            
            // Results ought to be ~consistent. ClipperLib2 reverses the direction compared to ClipperLib1
            Assert.True(((Math.Abs(clipped1[0][0].X - clipped2[0][0].X) <= 1) && (Math.Abs(clipped1[0][0].Y - clipped2[0][0].Y) <= 1)) ||
                        ((Math.Abs(clipped1[0][0].X - clipped2[0][1].X) <= 1) && (Math.Abs(clipped1[0][0].Y - clipped2[0][1].Y) <= 1)));
        }
    }
    
    [Test]
    public static void clipper2_openpath_parallelLines()
    {
        Path64 a1 = Clipper2Lib.Clipper.MakePath(new [] {10, 0, 20, 0});
        Path64 a2 = Clipper2Lib.Clipper.MakePath(new [] {30, 0, 40, 0});
        Path64 a3 = Clipper2Lib.Clipper.MakePath(new [] {50, 0, 60, 0});
        Path64 b1 = Clipper2Lib.Clipper.MakePath(new [] {20, 30, 30, 30});
        Path64 b2 = Clipper2Lib.Clipper.MakePath(new [] {40, 30, 50, 30});
        Path64 b3 = Clipper2Lib.Clipper.MakePath(new [] {60, 30, 70, 30});

        Clipper2Lib.Clipper64 c = new();
        c.AddOpenSubject(a1);
        c.AddOpenSubject(a2);
        c.AddOpenSubject(a3);
        c.AddOpenSubject(b1);
        c.AddOpenSubject(b2);
        c.AddOpenSubject(b3);

        Path64 bounds = Clipper2Lib.Clipper.MakePath(new[] { -100, -100, 100, 100 });
        c.AddClip(bounds);
        
        Paths64 o = new();
        c.Execute(Clipper2Lib.ClipType.Intersection, Clipper2Lib.FillRule.EvenOdd, new Paths64(), o);

        Assert.AreEqual(0, o.Count);
    }
    
    [Test]
    public static void clipper2_openpath_complextest()
    {
        Clipper2Lib.PathD clippingPath = Clipper2Lib.Clipper.MakePath(new[]
        {
            15.00000, -35.00000,
            13.60000, -34.95000,
            12.22000, -34.81000,
            10.84000, -34.56000,
            9.49000, -34.23000,
            8.16000, -33.79000,
            6.87000, -33.27000,
            5.61000, -32.66000,
            4.40000, -31.96000,
            3.24000, -31.18000,
            2.14000, -30.32000,
            1.11000, -29.39000,
            0.14000, -28.38000,
            -0.76000, -27.31000,
            -1.58000, -26.18000,
            -2.32000, -25.00000,
            -2.98000, -23.77000,
            -3.54000, -22.49000,
            -4.02000, -21.18000,
            -4.41000, -19.84000,
            -4.70000, -18.47000,
            -4.89000, -17.09000,
            -5.00000, -15.00000,
            -5.00000, -15.00000,
            -4.98000, -13.60000,
            -4.93000, -12.22000,
            -4.84000, -10.84000,
            -4.71000, -9.49000,
            -4.55000, -8.16000,
            -4.35000, -6.87000,
            -4.13000, -5.61000,
            -3.86000, -4.40000,
            -3.57000, -3.24000,
            -3.25000, -2.14000,
            -2.81000, -0.86000,
            -2.33000, 0.32000,
            -1.81000, 1.38000,
            -1.15000, 2.49000,
            -0.45000, 3.41000,
            0.41000, 4.23000,
            1.43000, 4.81000,
            2.48000, 5.00000,
            3.64000, 4.97000,
            4.78000, 4.88000,
            5.98000, 4.70000,
            7.08000, 4.46000,
            8.20000, 4.10000,
            9.19000, 3.58000,
            9.91000, 2.74000,
            9.95000, 2.48000,
            9.95000, 1.37000,
            9.96000, 0.27000,
            9.98000, -0.05000,
            9.99000, -1.16000,
            10.00000, -2.32000,
            10.00000, -3.67500,
            10.00000, -5.03000,
            9.87000, -6.14000,
            9.49000, -7.21000,
            8.89000, -8.16000,
            8.08000, -8.95000,
            7.11000, -9.53000,
            6.04000, -9.89000,
            5.00000, -10.00000,
            3.88000, -10.13000,
            2.81000, -10.51000,
            1.85000, -11.11000,
            1.06000, -11.92000,
            0.47000, -12.89000,
            0.11000, -13.96000,
            0.00000, -15.00000,
            0.13000, -16.12000,
            0.51000, -17.19000,
            1.11000, -18.15000,
            1.92000, -18.94000,
            2.89000, -19.53000,
            3.96000, -19.89000,
            5.00000, -20.00000,
            6.12000, -20.13000,
            7.19000, -20.51000,
            8.15000, -21.11000,
            8.94000, -21.92000,
            9.53000, -22.89000,
            9.89000, -23.96000,
            10.00000, -25.00000,
            10.13000, -26.12000,
            10.51000, -27.19000,
            11.11000, -28.15000,
            11.92000, -28.94000,
            12.89000, -29.53000,
            13.96000, -29.89000,
            15.00000, -30.00000,
            16.12000, -29.87000,
            17.19000, -29.49000,
            18.15000, -28.89000,
            18.94000, -28.08000,
            19.53000, -27.11000,
            19.89000, -26.04000,
            20.00000, -25.00000,
            20.13000, -23.88000,
            20.51000, -22.81000,
            21.11000, -21.85000,
            21.92000, -21.06000,
            22.89000, -20.47000,
            23.96000, -20.11000,
            25.00000, -20.00000,
            26.12000, -19.87000,
            27.19000, -19.49000,
            28.15000, -18.89000,
            28.94000, -18.08000,
            29.53000, -17.11000,
            29.89000, -16.04000,
            30.00000, -15.00000,
            29.87000, -13.88000,
            29.49000, -12.81000,
            28.89000, -11.85000,
            28.08000, -11.06000,
            27.11000, -10.47000,
            26.04000, -10.11000,
            25.00000, -10.00000,
            23.88000, -9.87000,
            22.81000, -9.49000,
            21.85000, -8.89000,
            21.06000, -8.08000,
            20.47000, -7.11000,
            20.11000, -6.04000,
            20.00000, -5.00000,
            19.87000, -3.88000,
            19.50000, -2.81000,
            18.89000, -1.85000,
            18.09000, -1.06000,
            17.13000, -0.47000,
            16.06000, -0.11000,
            15.02000, 0.00000,
            13.91000, 0.06000,
            12.77000, 0.27000,
            11.70000, 0.64000,
            10.72000, 1.25000,
            10.09000, 2.20000,
            10.05000, 2.50000,
            10.80000, 3.36000,
            11.83000, 3.79000,
            12.97000, 4.11000,
            14.18000, 4.36000,
            15.37000, 4.55000,
            16.48000, 4.69000,
            17.65000, 4.80000,
            18.88000, 4.89000,
            20.14000, 4.95000,
            21.44000, 4.99000,
            22.52000, 5.00000,
            23.83000, 4.89000,
            24.91000, 4.63000,
            25.96000, 4.23000,
            27.00000, 3.67000,
            27.99000, 2.98000,
            28.95000, 2.14000,
            29.86000, 1.18000,
            30.71000, 0.09000,
            31.35000, -0.86000,
            31.94000, -1.88000,
            32.49000, -2.96000,
            32.99000, -4.11000,
            33.44000, -5.30000,
            33.83000, -6.55000,
            34.17000, -7.83000,
            34.45000, -9.15000,
            34.68000, -10.50000,
            34.85000, -11.87000,
            34.95000, -13.26000,
            35.00000, -14.65000,
            35.00000, -15.00000,
            34.95000, -16.40000,
            34.81000, -17.78000,
            34.56000, -19.16000,
            34.23000, -20.51000,
            33.79000, -21.84000,
            33.27000, -23.13000,
            32.66000, -24.39000,
            31.96000, -25.60000,
            31.18000, -26.76000,
            30.32000, -27.86000,
            29.39000, -28.89000,
            28.38000, -29.86000,
            27.31000, -30.76000,
            26.18000, -31.58000,
            25.00000, -32.32000,
            23.77000, -32.98000,
            22.49000, -33.54000,
            21.18000, -34.02000,
            19.84000, -34.41000,
            18.47000, -34.70000,
            17.09000, -34.89000,
            15.70000, -34.99000,
        });

        PathD cutMe = Clipper2Lib.Clipper.MakePath(new[]
        {
            10.00000, -3.67500,
            40.00000, -3.67500,
        });

        Clipper2Lib.ClipperD d = new(Constants.roundingDecimalPrecision);
        d.AddOpenSubject(cutMe);
        d.AddClip(clippingPath);
        PathsD out_ = new();
        PathsD unused = new();
        d.Execute(Clipper2Lib.ClipType.Difference, Clipper2Lib.FillRule.EvenOdd, unused, out_);
        
        Assert.AreEqual(2, out_.Count);
        Assert.AreEqual(2, out_[0].Count);
        Assert.AreEqual(2, out_[1].Count);
    }
    
    [Test]
    public static void clipper2_openpath_test()
    {
        Path64 testPath = Clipper2Lib.Clipper.MakePath(new[]
        {
            -50000, -550000,
            -50000, -150000,
            650000, -150000
        });
        
        Path64 b = Clipper2Lib.Clipper.MakePath(new[]
        {
            300000,-800000,
            300000,0,
            500000,0,
            500000,-800000
        });
        
        Clipper2Lib.Clipper64 c = new() {PreserveCollinear = true};
        c.AddOpenSubject(testPath);
        c.AddClip(b);
        Paths64 unused = new();
        Paths64 topChords = new();
        c.Execute(Clipper2Lib.ClipType.Intersection, Clipper2Lib.FillRule.EvenOdd, unused, topChords);

        Path64 testPath2 = Clipper2Lib.Clipper.MakePath(new[]
        {
            650000,-150000,
            650000,-550000,
            -50000,-550000
        });
        
        c.Clear();
        c.AddOpenSubject(testPath2);
        c.AddClip(b);
        Paths64 bottomChords = new();
        c.Execute(Clipper2Lib.ClipType.Intersection, Clipper2Lib.FillRule.EvenOdd, unused, bottomChords);
        
        Path64 testPath3 = Clipper2Lib.Clipper.MakePath(new[]
        {
            300000,-800000,
            300000,0
        });

        Path64 a = Clipper2Lib.Clipper.MakePath(new[]
        {
            -50000, -550000,
            -50000, -150000,
            650000, -150000,
            650000, -550000
        });

        c.Clear();
        c.AddOpenSubject(testPath3);
        c.AddClip(a);
        Paths64 leftChords = new();
        c.Execute(Clipper2Lib.ClipType.Intersection, Clipper2Lib.FillRule.EvenOdd, unused, leftChords);

        Path64 testPath4 = Clipper2Lib.Clipper.MakePath(new[]
        {
            300000,0,
            500000,0,
            500000,-800000,
            300000,-800000
        });

        c.Clear();
        c.AddOpenSubject(testPath4);
        c.AddClip(a);
        Paths64 rightChords = new();
        c.Execute(Clipper2Lib.ClipType.Intersection, Clipper2Lib.FillRule.EvenOdd, unused, rightChords);

        Assert.AreEqual(1, leftChords.Count);
        Assert.AreEqual(2, leftChords[0].Count);
        Assert.AreEqual(1, rightChords.Count);
        Assert.AreEqual(2, rightChords[0].Count);
        Assert.AreEqual(1, bottomChords.Count);
        Assert.AreEqual(2, bottomChords[0].Count);
        Assert.AreEqual(1, topChords.Count);
        Assert.AreEqual(2, topChords[0].Count);
    }

    [Test]
    public static void OpenPathOffsetTest1()
    {
        int[] pointData = new[]
        {
            -360,-100,
            120,260,
            300,-140,
            -340,160
        };

        Path64 sourcePath = Clipper2Lib.Clipper.MakePath(pointData);
        
        Clipper2Lib.ClipperOffset co = new(miterLimit:2, arcTolerance:0.25);
        co.AddPath(sourcePath, Clipper2Lib.JoinType.Miter, Clipper2Lib.EndType.Square);
        Paths64 resizedPolyData = new();
        co.Execute(50, resizedPolyData);

        double area = resizedPolyData.Sum(t => Clipper2Lib.Clipper.Area(t));
        Assert.AreEqual(86501.5, area);
        Assert.AreEqual(2, resizedPolyData.Count);
        Assert.AreEqual(11, resizedPolyData[0].Count);
        Assert.AreEqual(3, resizedPolyData[1].Count);
    }

    [Test]
    public static void OpenPathOffsetTest2()
    {
        int[] pointData = new[]
        {
            0, 0,
            0, 10000,
            0, 20000,
            0, 30000,
            0, 40000,
            0, 50000,
            0, 60000,
            0, 70000,
            0, 80000,
            0, 90000,
            0, 100000,
            0, 110000,
            0, 120000,
            0, 130000,
            0, 140000,
            0, 150000,
            0, 160000,
            0, 170000,
            0, 180000,
            0, 190000,
            0, 200000,
            0, 210000,
            0, 220000,
            0, 230000,
            0, 240000,
            0, 250000,
            0, 260000,
            0, 270000,
            0, 280000,
            0, 290000,
            0, 300000,
            10000, 300000,
            20000, 300000,
            30000, 300000,
            40000, 300000,
            50000, 300000,
            60000, 300000,
            70000, 300000,
            80000, 300000,
            90000, 300000,
            100000, 300000,
            110000, 300000,
            120000, 300000,
            130000, 300000,
            140000, 300000,
            150000, 300000,
            160000, 300000,
            170000, 300000,
            180000, 300000,
            190000, 300000,
            200000, 300000,
            210000, 300000,
            220000, 300000,
            230000, 300000,
            240000, 300000,
            250000, 300000,
            260000, 300000,
            270000, 300000,
            280000, 300000,
            290000, 300000,
            300000, 300000,
            300000, 290000,
            300000, 280000,
            300000, 270000,
            300000, 260000,
            300000, 250000,
            300000, 240000,
            300000, 230000,
            300000, 220000,
            300000, 210000,
            300000, 200000,
            300000, 190000,
            300000, 180000,
            300000, 170000,
            300000, 160000,
            300000, 150000,
            300000, 140000,
            300000, 130000,
            300000, 120000,
            300000, 110000,
            300000, 100000,
            300000, 90000,
            300000, 80000,
            300000, 70000,
            300000, 60000,
            300000, 50000,
            300000, 40000,
            300000, 30000,
            300000, 20000,
            300000, 10000,
            300000, 0,
            290000, 0,
            280000, 0,
            270000, 0,
            260000, 0,
            250000, 0,
            240000, 0,
            230000, 0,
            220000, 0,
            210000, 0,
            200000, 0,
            190000, 0,
            180000, 0,
            170000, 0,
            160000, 0,
            150000, 0,
            140000, 0,
            130000, 0,
            120000, 0,
            110000, 0,
            100000, 0,
            90000, 0,
            80000, 0,
            70000, 0,
            60000, 0,
            50000, 0,
            40000, 0,
            30000, 0,
            20000, 0,
            10000, 0,
            0, 0
        };

        Path64 sourcePath = Clipper2Lib.Clipper.MakePath(pointData);
        
        Clipper2Lib.ClipperOffset co = new() {PreserveCollinear = true, ReverseSolution = true};
        co.AddPath(sourcePath, Clipper2Lib.JoinType.Miter, Clipper2Lib.EndType.Polygon);
        Paths64 resizedPolyData = new();
        co.Execute(Convert.ToDouble(6 * 10000), resizedPolyData);
        
        double area = resizedPolyData.Sum(t => Clipper2Lib.Clipper.Area(t));
        Assert.AreEqual(176400000000, area);
        Assert.AreEqual(1, resizedPolyData.Count);
        Assert.AreEqual(120, resizedPolyData[0].Count);
    }
    
    [Test]
    public static void ClosedPathOffsetTest()
    {
        Path64 closedPath_forOffsetTest = new()
        {
            new (0, 0),
            new (0, 8),
            new (4, 8),
            new (4, 4),
            new (16, 4),
            new (16, 16),
            new (4, 16),
            new (4, 10),
            new (0, 10),
            new (0,20),
            new (20,20),
            new (20,0)

        };
        Path test1 = new();
        for (int pt = 0; pt < closedPath_forOffsetTest.Count; pt++)
        {
            test1.Add(new ClipperLib1.IntPoint(closedPath_forOffsetTest[pt].X, closedPath_forOffsetTest[pt].Y));
        }
        ClipperLib1.ClipperOffset co1 = new();
        co1.AddPath(test1, ClipperLib1.JoinType.jtMiter, ClipperLib1.EndType.etClosedPolygon);
        Paths c1up = new();
        co1.Execute(ref c1up, 2.0);
        
        double area = c1up.Sum(t => ClipperLib1.Clipper.Area(t));
        Assert.AreEqual(512, area);
        Assert.AreEqual(2, c1up.Count);
        Assert.AreEqual(4, c1up[0].Count);
        Assert.AreEqual(4, c1up[1].Count);
        
        co1.Clear();
        co1.AddPaths(c1up, ClipperLib1.JoinType.jtMiter, ClipperLib1.EndType.etClosedPolygon);
        Paths c1down = new();
        co1.Execute(ref c1down, -2.0);

        double area1 = c1down.Sum(t => ClipperLib1.Clipper.Area(t));
        Assert.AreEqual(256, area1);
        Assert.AreEqual(2, c1down.Count);
        Assert.AreEqual(4, c1down[0].Count);
        Assert.AreEqual(4, c1down[1].Count);
        
        Clipper2Lib.ClipperOffset co2 = new() {PreserveCollinear = true, ReverseSolution = true};
        co2.AddPath(closedPath_forOffsetTest, Clipper2Lib.JoinType.Miter, Clipper2Lib.EndType.Polygon);
        Paths64 c2up = new();
        co2.Execute(2.0, c2up);
        
        double area2 = c2up.Sum(t => Clipper2Lib.Clipper.Area(t));
        Assert.AreEqual(512, area2);
        Assert.AreEqual(2, c2up.Count);
        Assert.AreEqual(5, c2up[0].Count);
        Assert.AreEqual(5, c2up[1].Count);
        
        co2.Clear();
        co2.AddPaths(c2up, Clipper2Lib.JoinType.Miter, Clipper2Lib.EndType.Polygon);
        Paths64 c2down = new();
        co2.Execute(-2.0, c2down);
        
        double area3 = c2down.Sum(t => Clipper2Lib.Clipper.Area(t));
        Assert.AreEqual(-256, area3);
        Assert.AreEqual(2, c2down.Count);
        Assert.AreEqual(5, c2down[0].Count);
        Assert.AreEqual(5, c2down[1].Count);
    }
    
    const double clipper2_keyhole_sizing = 500;

    [Test]
    public static void clipper2_H_subtract_test()
    {
        double[] subj0 = new double[] { 0, 0, 9, 0, 9, 12, 0, 12 };
        //Array.Reverse(subj0);
        double[] subj1 = new double[] { 0, 2, 9, 2, 9, 10, 0, 10 };
        //Array.Reverse(subj1);
        PathsD subj = new () { Clipper2Lib.Clipper.MakePath(subj0),
            Clipper2Lib.Clipper.MakePath(subj1)};

        double[] clip0 = new double[] { 0, 3, 9, 3, 9, 4, 0, 4 };
        //Array.Reverse(clip0);
        double[] clip1 = new double[] { 0, 8, 9, 8, 9, 9, 0, 9 };
        //Array.Reverse(clip1);
        double[] clip2 = new double[] { 2, 3, 3, 3, 3, 9, 2, 9 };
        //Array.Reverse(clip2);
        double[] clip3 = new double[] { 6, 3, 7, 3, 7, 9, 6, 9 };
        //Array.Reverse(clip3);
        
        PathsD clip = new() { Clipper2Lib.Clipper.MakePath(clip0),
            Clipper2Lib.Clipper.MakePath(clip1),
            Clipper2Lib.Clipper.MakePath(clip2),
            Clipper2Lib.Clipper.MakePath(clip3) };

        subj = Clipper2Lib.Clipper.Union(subj, Clipper2Lib.FillRule.NonZero);
        clip = Clipper2Lib.Clipper.Union(clip, Clipper2Lib.FillRule.NonZero);

        PathsD sol = Clipper2Lib.Clipper.Difference(subj, clip, Clipper2Lib.FillRule.NonZero);

        SvgWriter svgSrc,svgDst;
        svgSrc = new();
        svgDst = new();
        SvgUtils.AddSolution(svgSrc, subj, true);
        SvgUtils.AddSolution(svgSrc, clip, true);
        SvgUtils.SaveToFile(svgSrc, "svgSrc.svg", Clipper2Lib.FillRule.NonZero, 800, 600, 10);

        SvgUtils.AddSolution(svgDst, sol, true);
        SvgUtils.SaveToFile(svgDst, "svgDst.svg", Clipper2Lib.FillRule.NonZero, 800, 600, 10);
        
        double area = sol.Sum(t => Clipper2Lib.Clipper.Area(t));
        Assert.AreEqual(82, area);
        Assert.AreEqual(5, sol.Count);
        Assert.AreEqual(8, sol[0].Count);
        Assert.AreEqual(4, sol[1].Count);
        Assert.AreEqual(4, sol[2].Count);
        Assert.AreEqual(4, sol[3].Count);
        Assert.AreEqual(4, sol[4].Count);
    }
    
    [Test]
    public static void clipper2_chordTest()
    {
        Path64 testPath = Clipper2Lib.Clipper.MakePath(new[]
        {
            -50000, -550000,
            -50000, -150000,
            650000, -150000
        });
        
        Path64 b = Clipper2Lib.Clipper.MakePath(new[]
        {
            300000,-800000,
            300000,0,
            500000,0,
            500000,-800000
        });
        
        Clipper2Lib.Clipper64 c = new() {PreserveCollinear = true};
        c.AddOpenSubject(testPath);
        c.AddClip(b);
        Paths64 unused = new();
        Paths64 topChords = new();
        c.Execute(Clipper2Lib.ClipType.Intersection, Clipper2Lib.FillRule.EvenOdd, unused, topChords);

        Path64 testPath2 = Clipper2Lib.Clipper.MakePath(new[]
        {
            650000,-150000,
            650000,-550000,
            -50000,-550000
        });
        
        c.Clear();
        c.AddOpenSubject(testPath2);
        c.AddClip(b);
        Paths64 bottomChords = new();
        c.Execute(Clipper2Lib.ClipType.Intersection, Clipper2Lib.FillRule.EvenOdd, unused, bottomChords);
        
        Path64 testPath3 = Clipper2Lib.Clipper.MakePath(new[]
        {
            300000,-800000,
            300000,0
        });

        Path64 a = Clipper2Lib.Clipper.MakePath(new[]
        {
            -50000, -550000,
            -50000, -150000,
            650000, -150000,
            650000, -550000
        });

        c.Clear();
        c.AddOpenSubject(testPath3);
        c.AddClip(a);
        Paths64 leftChords = new();
        c.Execute(Clipper2Lib.ClipType.Intersection, Clipper2Lib.FillRule.EvenOdd, unused, leftChords);

        Path64 testPath4 = Clipper2Lib.Clipper.MakePath(new[]
        {
            300000,0,
            500000,0,
            500000,-800000,
            300000,-800000
        });

        c.Clear();
        c.AddOpenSubject(testPath4);
        c.AddClip(a);
        Paths64 rightChords = new();
        c.Execute(Clipper2Lib.ClipType.Intersection, Clipper2Lib.FillRule.EvenOdd, unused, rightChords);

        Assert.AreEqual(1, leftChords.Count);
        Assert.AreEqual(2, leftChords[0].Count);
        Assert.AreEqual(1, rightChords.Count);
        Assert.AreEqual(2, rightChords[0].Count);
        Assert.AreEqual(1, bottomChords.Count);
        Assert.AreEqual(2, bottomChords[0].Count);
        Assert.AreEqual(1, topChords.Count);
        Assert.AreEqual(2, topChords[0].Count);
    }
    
    [Test]
    public static void clipper2_collinearTest()
    {
        Path64 collinear = new()
        {
            new(-10, -10),
            new(-10, 0),
            new(-10, 10),
            new(0, 10),
            new(10, 10),
            new(10, 0),
            new(10, -10),
            new(0, -10),
            new(-10, -10),
        };

        Clipper2Lib.Clipper64 c = new() {PreserveCollinear = true};
        c.AddSubject(collinear);
        c.AddClip(collinear);
        Paths64 output = new();
        c.Execute(Clipper2Lib.ClipType.Union, Clipper2Lib.FillRule.EvenOdd, output);

        double area = output.Sum(t => Clipper2Lib.Clipper.Area(t));
        Assert.AreEqual(400, area);
        Assert.AreEqual(1, output.Count);
        Assert.AreEqual(8, output[0].Count);
    }

    [Test]
    public static void clipper2_collinearOffsetTest()
    {
        Path64 collinear = Clipper2Lib.Clipper.MakePath(new []
        {
            -10, -10,
            -10, 0,
            -10, 10,
            0, 10,
            10, 10,
            10, 0,
            10, -10,
            0, -10,
            -10, -10,
        });

        Clipper2Lib.ClipperOffset co = new() {PreserveCollinear = true, ReverseSolution = true};
        co.AddPath(collinear, Clipper2Lib.JoinType.Miter, Clipper2Lib.EndType.Polygon);
        Paths64 temp = new();
        co.Execute(1.0, temp);
        
        double area = temp.Sum(t => Clipper2Lib.Clipper.Area(t));
        Assert.AreEqual(484, area);
        Assert.AreEqual(1, temp.Count);
        Assert.AreEqual(8, temp[0].Count);
    }

    [Test]
    public static void clipper2_unionTest()
    {
        Path64 simpleFirstPath = new()
        {
            new(100000, -300000),
            new(100000, 0),
            new(200000, 0),
            new(200000, -300000),
            new(100000, -300000),
        };
        
        Path64 simpleSecondPath = new()
        {
            new (0,-200000),
            new (0,-100000),
            new (300000,-100000),
            new (300000,-200000),
            new (0,-200000),
        };
        
        Clipper2Lib.Clipper64 cs = new();
        
        cs.AddSubject(simpleFirstPath);
        cs.AddClip(simpleSecondPath);
        
        Paths64 simpleOutputPoints = new();
        cs.Execute(Clipper2Lib.ClipType.Union, Clipper2Lib.FillRule.EvenOdd, simpleOutputPoints);
        
        double area = simpleOutputPoints.Sum(t => Clipper2Lib.Clipper.Area(t));
        Assert.AreEqual(50000000000, area);
        Assert.AreEqual(1, simpleOutputPoints.Count);
        Assert.AreEqual(12, simpleOutputPoints[0].Count);

        Path64 firstPath = new() {
        new (100000,-300000),
        new (100000,-290000),
        new (100000,-280000),
        new (100000,-270000),
        new (100000,-260000),
        new (100000,-250000),
        new (100000,-240000),
        new (100000,-230000),
        new (100000,-220000),
        new (100000,-210000),
        new (100000,-200000),
        new (100000,-190000),
        new (100000,-180000),
        new (100000,-170000),
        new (100000,-160000),
        new (100000,-150000),
        new (100000,-140000),
        new (100000,-130000),
        new (100000,-120000),
        new (100000,-110000),
        new (100000,-100000),
        new (100000,-90000),
        new (100000,-80000),
        new (100000,-70000),
        new (100000,-60000),
        new (100000,-50000),
        new (100000,-40000),
        new (100000,-30000),
        new (100000,-20000),
        new (100000,-10000),
        new (100000,0),
        new (110000,0),
        new (120000,0),
        new (130000,0),
        new (140000,0),
        new (150000,0),
        new (160000,0),
        new (170000,0),
        new (180000,0),
        new (190000,0),
        new (200000,0),
        new (200000,-10000),
        new (200000,-20000),
        new (200000,-30000),
        new (200000,-40000),
        new (200000,-50000),
        new (200000,-60000),
        new (200000,-70000),
        new (200000,-80000),
        new (200000,-90000),
        new (200000,-100000),
        new (200000,-110000),
        new (200000,-120000),
        new (200000,-130000),
        new (200000,-140000),
        new (200000,-150000),
        new (200000,-160000),
        new (200000,-170000),
        new (200000,-180000),
        new (200000,-190000),
        new (200000,-200000),
        new (200000,-210000),
        new (200000,-220000),
        new (200000,-230000),
        new (200000,-240000),
        new (200000,-250000),
        new (200000,-260000),
        new (200000,-270000),
        new (200000,-280000),
        new (200000,-290000),
        new (200000,-300000),
        new (190000,-300000),
        new (180000,-300000),
        new (170000,-300000),
        new (160000,-300000),
        new (150000,-300000),
        new (140000,-300000),
        new (130000,-300000),
        new (120000,-300000),
        new (110000,-300000),
        new (100000,-300000),
        };

        Path64 secondPath = new() {
        new (0,-200000),
        new (0,-190000),
        new (0,-180000),
        new (0,-170000),
        new (0,-160000),
        new (0,-150000),
        new (0,-140000),
        new (0,-130000),
        new (0,-120000),
        new (0,-110000),
        new (0,-100000),
        new (10000,-100000),
        new (20000,-100000),
        new (30000,-100000),
        new (40000,-100000),
        new (50000,-100000),
        new (60000,-100000),
        new (70000,-100000),
        new (80000,-100000),
        new (90000,-100000),
        new (100000,-100000),
        new (110000,-100000),
        new (120000,-100000),
        new (130000,-100000),
        new (140000,-100000),
        new (150000,-100000),
        new (160000,-100000),
        new (170000,-100000),
        new (180000,-100000),
        new (190000,-100000),
        new (200000,-100000),
        new (210000,-100000),
        new (220000,-100000),
        new (230000,-100000),
        new (240000,-100000),
        new (250000,-100000),
        new (260000,-100000),
        new (270000,-100000),
        new (280000,-100000),
        new (290000,-100000),
        new (300000,-100000),
        new (300000,-110000),
        new (300000,-120000),
        new (300000,-130000),
        new (300000,-140000),
        new (300000,-150000),
        new (300000,-160000),
        new (300000,-170000),
        new (300000,-180000),
        new (300000,-190000),
        new (300000,-200000),
        new (290000,-200000),
        new (280000,-200000),
        new (270000,-200000),
        new (260000,-200000),
        new (250000,-200000),
        new (240000,-200000),
        new (230000,-200000),
        new (220000,-200000),
        new (210000,-200000),
        new (200000,-200000),
        new (190000,-200000),
        new (180000,-200000),
        new (170000,-200000),
        new (160000,-200000),
        new (150000,-200000),
        new (140000,-200000),
        new (130000,-200000),
        new (120000,-200000),
        new (110000,-200000),
        new (100000,-200000),
        new (90000,-200000),
        new (80000,-200000),
        new (70000,-200000),
        new (60000,-200000),
        new (50000,-200000),
        new (40000,-200000),
        new (30000,-200000),
        new (20000,-200000),
        new (10000,-200000),
        new (0,-200000),
        };

        Clipper2Lib.Clipper64 c = new();

        c.AddSubject(firstPath);
        c.AddClip(secondPath);

        Paths64 outputPoints = new();
        c.Execute(Clipper2Lib.ClipType.Union, Clipper2Lib.FillRule.EvenOdd, outputPoints);        
        double area1 = outputPoints.Sum(t => Clipper2Lib.Clipper.Area(t));
        Assert.AreEqual(50000000000, area1);
        Assert.AreEqual(1, outputPoints.Count);
        Assert.AreEqual(120, outputPoints[0].Count);
    }

    [Test]
    public static void clipper2_notTest()
    {
        Path64 firstPath = new() {
            new (-250000,-250000),
            new (-250000,250000),
            new (250000,250000),
            new (250000,-250000),
            new (-250000,-250000),
        };
        
        Path64 secondPath = new() {
        new (-150000,-150000),
        new (-150000,150000),
        new (150000,150000),
        new (150000,-150000),
        new (-150000,-150000),
        };
        
        Clipper2Lib.Clipper64 c = new();
        
        c.AddSubject(firstPath);
        c.AddClip(secondPath);
        
        Paths64 outputPoints = new();
        
        c.Execute(Clipper2Lib.ClipType.Difference, Clipper2Lib.FillRule.EvenOdd, outputPoints);
        
        double area = outputPoints.Sum(t => Clipper2Lib.Clipper.Area(t));
        Assert.AreEqual(160000000000, area);
        Assert.AreEqual(2, outputPoints.Count);
        Assert.AreEqual(4, outputPoints[0].Count);
        Assert.AreEqual(4, outputPoints[1].Count);
    }

    [Test]
    public static void clipper2_edgeOffsetTest()
    {
        Path64 edge = new()
        {
            new(-100000, 99500),
            new(-100000, 200500)
        };

        Clipper2Lib.ClipperOffset co = new() {PreserveCollinear = true, ReverseSolution = true};
        co.AddPath(edge, Clipper2Lib.JoinType.Miter, Clipper2Lib.EndType.Square);
        Paths64 p = new();
        co.Execute(500, p);
        
        double area = p.Sum(t => Clipper2Lib.Clipper.Area(t));
        Assert.AreEqual(-50750000, area);
        Assert.AreEqual(1, p.Count);
        Assert.AreEqual(5, p[0].Count);
    }

    private static void clipper2_zFillTest(Point64 bot1, Point64 top1, Point64 bot2, Point64 top2, ref Point64 pt)
    {
        pt = new(pt.X, pt.Y, -1);
    }
    
    [Test]
    public static void clipper2_zFillCallbackTest()
    {
        Path64 outer = new()
        {
            new(-1000, -1000),
            new(-1000, 1000),
            new(1000, 1000),
            new(1000, -1000),
        };
        
        Path64 cutter = new()
        {
            new(-100, -1100),
            new(-100, -900),
            new(100, -900),
            new(100, -1100),
        };

        Clipper2Lib.Clipper64 c = new();

        c.ZCallback = clipper2_zFillTest;

        c.AddSubject(outer);
        c.AddClip(cutter);

        Paths64 solution = new();
        c.Execute(Clipper2Lib.ClipType.Intersection, Clipper2Lib.FillRule.EvenOdd, solution);
        
        double area = solution.Sum(t => Clipper2Lib.Clipper.Area(t));
        Assert.AreEqual(20000, area);
        Assert.AreEqual(1, solution.Count);
        Assert.AreEqual(4, solution[0].Count);
    }
    
    [Test]
    public static void clipper2_coincident_openPathTest()
    {
        // All solutionN paths should be zero count due to the open path clipping.
        Path64 lPoly = new()
        {
            new(-1000, -1000),
            new(-1000, 1000),
            new(1000, 1000),
            new(1000, -1000)
        };

        Path64 t = new()
        {
            new (-1000, -1100),
            new (-1000, 500),
        };

        Clipper2Lib.Clipper64 c = new();
        c.AddClip(lPoly);
        c.AddOpenSubject(t);

        Paths64 open = new();
        Paths64 solution = new();
        c.Execute(Clipper2Lib.ClipType.Intersection, Clipper2Lib.FillRule.EvenOdd, solution, open);
        
        double area = solution.Sum(t => Clipper2Lib.Clipper.Area(t));
        Assert.AreEqual(0, area);
        Assert.AreEqual(0, solution.Count);
        
        Path64 t2 = new()
        {
            new (-1000, -1100),
            new (-1000, 500),
            new (-900, 500),
        };

        Clipper2Lib.Clipper64 c2 = new();
        c2.AddClip(lPoly);
        c2.AddOpenSubject(t2);

        Paths64 open2 = new();
        Paths64 solution2 = new();
        c2.Execute(Clipper2Lib.ClipType.Intersection, Clipper2Lib.FillRule.EvenOdd, solution2, open2);

        double area1 = solution2.Sum(t => Clipper2Lib.Clipper.Area(t));
        Assert.AreEqual(0, area1);
        Assert.AreEqual(0, solution2.Count);

        Path64 t2b = new()
        {
            new (-900, 500),
            new (-1000, 500),
            new (-1000, -1100),
        };

        Clipper2Lib.Clipper64 c2b = new();
        c2b.AddClip(lPoly);
        c2b.AddOpenSubject(t2b);

        Paths64 open2b = new();
        Paths64 solution2b = new();
        c2b.Execute(Clipper2Lib.ClipType.Intersection, Clipper2Lib.FillRule.EvenOdd, solution2b, open2b);

        double area2 = solution2b.Sum(t => Clipper2Lib.Clipper.Area(t));
        Assert.AreEqual(0, area2);
        Assert.AreEqual(0, solution2b.Count);

        Path64 t3 = new();
        int x = 0;
        int y = -1100;
        while (y < 1200)
        {
            t3.Add(new(x, y));
            y += 100;
        }
        Clipper2Lib.Clipper64 c3 = new();
        c3.ZCallback = clipper2_zFillTest;
        
        c3.AddClip(lPoly);
        c3.AddOpenSubject(t3);

        Paths64 open3 = new();
        Paths64 solution3 = new();
        c3.Execute(Clipper2Lib.ClipType.Intersection, Clipper2Lib.FillRule.EvenOdd, solution3, open3 );
        
        double area3 = solution3.Sum(t => Clipper2Lib.Clipper.Area(t));
        Assert.AreEqual(0, area3);
        Assert.AreEqual(0, solution3.Count);
    }
    
    [Test]
    public static void clipper2_keyHole_test2()
    {
        Path64 lPoly = new()
        {
            new Point64(200000, 0),
            new Point64(200000, 1100000),
            new Point64(1000000, 1100000),
            new Point64(1000000, 800000),
            new Point64(800000, 800000),
            new Point64(800000, 0),
            new Point64(200000, 0)
        };

        Paths64 t = new()
        {
            new Path64()
            {
                new(800000, 800000),
                new(800000, 1100000)
            }
        };
        
        // Turn the new edges into cutters and slice. Not terribly elegant and we're relying on rounding to squash notches later.
        Clipper2Lib.ClipperOffset co = new() {PreserveCollinear = true, ReverseSolution = true};
        co.AddPaths(t, Clipper2Lib.JoinType.Miter, Clipper2Lib.EndType.Square);

        Paths64 cutters = new();
        co.Execute(2.0, cutters);

        double area = cutters.Sum(t => Clipper2Lib.Clipper.Area(t));
        Assert.AreEqual(-600004, area);
        Assert.AreEqual(1, cutters.Count);
        Assert.AreEqual(5, cutters[0].Count);
        
        Clipper2Lib.Clipper64 c = new();

        c.AddSubject(lPoly);

        // Take first cutter only - we only cut once, no matter how many potential cutters we have.
        c.AddClip(cutters[0]);
        Paths64 f = new();
        c.Execute(Clipper2Lib.ClipType.Difference, Clipper2Lib.FillRule.EvenOdd, f);
        double area1 = f.Sum(t => Clipper2Lib.Clipper.Area(t));
        Assert.AreEqual(719999399999, area1);
        Assert.AreEqual(2, f.Count);
        Assert.AreEqual(4, f[0].Count);
        Assert.AreEqual(6, f[1].Count);
    }
    
    [Test]
    public static void clipper2_openPath_clipTest1()
    {
        Path64 lPoly = new () {
        new(0,0),
        new(0,200000),
        new(200000,200000),
        new(200000,500000),
        new(0,500000),
        new(0,1100000),
        new(1000000,1100000),
        new(1000000,800000),
        new(800000,800000),
        new(800000,600000),
        new(1000000,600000),
        new(1000000,0),
        new(0,0)
        };

        Path64 t = new () {
        new Point64(0,200000),
        new Point64(0,-9800000)
        };

        Clipper2Lib.Clipper64 c = new();

        c.AddOpenSubject(t);
        c.AddClip(lPoly);

        Clipper2Lib.PolyTree64 pt = new();
        Paths64 p = new();

        c.Execute(Clipper2Lib.ClipType.Intersection, Clipper2Lib.FillRule.EvenOdd, pt, p);
        
        Assert.AreEqual(1, p.Count);
        Assert.AreEqual(2, p[0].Count);
    }
    
    [Test]
    public static void clipper2_openPath_clipTest2()
    {
        Paths64 rays2 = new()
        {
            new Path64() {new(100000, 200000), new(100000, -9800000)}
        };

        Paths64 rays = new()
        {
            new Path64() {new(100000, 500000), new(10100000, 500000)}
        };

        Paths64 collisionPaths = new()
        {
            new Path64()
            {
                new(0, 0),
                new(0, 500000),
                new(100000, 500000),
                new(100000, 200000),
                new(600000, 200000),
                new(600000, 800000),
                new(1200000, 800000),
                new(1200000, 0),
                new(0, 0)

            }
        };

        Clipper2Lib.Clipper64 c = new ();
        c.AddOpenSubject(rays);
        c.AddClip(collisionPaths);
        Clipper2Lib.PolyTree64 pt = new ();
        Paths64 solution = new();
        c.Execute(Clipper2Lib.ClipType.Intersection, Clipper2Lib.FillRule.EvenOdd, pt, solution);
        
        Assert.AreEqual(1, solution.Count);
        Assert.AreEqual(2, solution[0].Count);
    }
    
    [Test]
    public static void clipper2_keyHole_test1()
    {
        Console.WriteLine("Clipper2 Test1");
        Path64 outer = new()
        {
            new Point64(-200000, -200000),
            new Point64(200000, -200000),
            new Point64(200000, 200000),
            new Point64(-200000, 200000),
            new Point64(-200000, -200000)
        };

        Path64 inner1 = new()
        {
            new Point64(-100000, -100000),
            new Point64(-100000, 100000),
            new Point64(100000, 100000),
            new Point64(100000, -100000),
            new Point64(-100000, -100000)
        };

        Paths64 kHSource = new()
        {
            outer,
            inner1
        };
        
        Clipper2Lib.ClipperOffset co = new() {PreserveCollinear = true, ReverseSolution = true};
        co.AddPaths(kHSource, Clipper2Lib.JoinType.Miter, Clipper2Lib.EndType.Polygon);
        // ClipperLib2 specifies full width in offset for open path, unlike version 1
        Paths64 out_ = new();
        co.Execute(2*clipper2_keyhole_sizing, out_);
        
        double area = out_.Sum(t => Clipper2Lib.Clipper.Area(t));
        Assert.AreEqual(-122400000000, area);
        Assert.AreEqual(2, out_.Count);
        Assert.AreEqual(4, out_[0].Count);
        Assert.AreEqual(4, out_[1].Count);
    }
    
    [Test]
    public static void clipper2_offsetTest()
    {
        Path64 lPoly = new ()
        {
            new Point64(0, 0),
            new Point64(0, 500000),
            new Point64(100000, 500000),
            new Point64(100000, 200000),
            new Point64(600000, 200000),
            new Point64(600000, 0),
            new Point64(0, 0)
        };

        Path64 newEdge = new()
        {
            new Point64(100000, 200000),
            new Point64(100000, 0)
        };

        Paths64 newEdges = new()
        {
            newEdge
        };

        Clipper2Lib.ClipperOffset co = new() {PreserveCollinear = true, ReverseSolution = true};
        co.AddPaths(newEdges, Clipper2Lib.JoinType.Miter, Clipper2Lib.EndType.Square);
        // ClipperLib2 specifies full width in offset for open path, unlike version 1
        Paths64 cutters = new();
        co.Execute( 2.0, cutters);
        
        double area = cutters.Sum(t => Clipper2Lib.Clipper.Area(t));
        Assert.AreEqual(-400004, area);
        Assert.AreEqual(1, cutters.Count);
        Assert.AreEqual(5, cutters[0].Count);
        
        Paths64 solution = new();
        Clipper2Lib.Clipper64 c = new ();
        c.AddSubject(lPoly);
        c.AddClip(cutters);
        c.Execute(Clipper2Lib.ClipType.Difference, Clipper2Lib.FillRule.EvenOdd, solution);
        
        double area1 = solution.Sum(t => Clipper2Lib.Clipper.Area(t));
        Assert.AreEqual(149999599999, area1);
        Assert.AreEqual(2, solution.Count);
        Assert.AreEqual(6, solution[0].Count);
        Assert.AreEqual(4, solution[1].Count);
    }
    
    [Test]
    public static void clipper2_leftChordTest()
    {
        Path64 testPath = new () {
         new(-200000,0),
         new(-300000,0),
         new(-310453,-548),
         new(-320791,-2185),
         new(-330902,-4894),
         new(-340674,-8645),
         new(-350000,-13397),
         new(-358779,-19098),
         new(-366913,-25686),
         new(-374314,-33087),
         new(-380902,-41221),
         new(-386603,-50000),
         new(-391355,-59326),
         new(-395106,-69098),
         new(-397815,-79209),
         new(-399452,-89547),
         new(-400000,-100000),
         new(-400000,-700000),
         new(-400000,-700000),
         new(-398907,-710396),
         new(-397815,-720791),
         new(-395106,-730902),
         new(-391355,-740674),
         new(-386603,-750000),
         new(-380902,-758779),
         new(-374314,-766913),
         new(-366913,-774314),
         new(-358779,-780902),
         new(-350000,-786603),
         new(-340674,-791355),
         new(-330902,-795106),
         new(-320791,-797815),
         new(-310453,-799452),
         new(-300000,-800000),
         };

        Path64 a = new () {
          new(-460000,-390000),
          new(-459955,-395235),
          new(-459594,-405679),
          new(-458873,-416047),
          new(-457796,-426288),
          new(-456369,-436353),
          new(-454102,-448610),
          new(-451315,-460421),
          new(-448030,-471696),
          new(-444272,-482349),
          new(-440068,-492300),
          new(-435451,-501472),
          new(-429416,-511353),
          new(-422903,-519904),
          new(-414796,-528076),
          new(-406257,-534189),
          new(-396132,-538540),
          new(-385806,-540000),
          new(-381067,-540036),
          new(-369255,-540441),
          new(-357570,-541292),
          new(-346100,-542583),
          new(-334932,-544305),
          new(-324151,-546443),
          new(-313840,-548983),
          new(-304076,-551904),
          new(-293187,-555882),
          new(-283312,-560333),
          new(-273218,-566059),
          new(-264802,-572279),
          new(-257399,-579871),
          new(-252063,-588852),
          new(-250000,-599117),
          new(-242658,-614841),
          new(-234819,-621423),
          new(-225801,-626830),
          new(-216572,-631139),
          new(-206066,-635097),
          new(-196418,-638096),
          new(-186036,-640798),
          new(-175000,-643183),
          new(-163393,-645233),
          new(-151303,-646931),
          new(-141346,-648029),
          new(-131187,-648888),
          new(-120876,-649505),
          new(-110463,-649876),
          new(-100000,-650000),
          new(-94765,-649967),
          new(-84321,-649701),
          new(-73953,-649169),
          new(-63712,-648376),
          new(-53647,-647324),
          new(-41390,-645654),
          new(-29579,-643601),
          new(-18304,-641180),
          new(-7651,-638411),
          new(2300,-635314),
          new(13206,-631197),
          new(22873,-626688),
          new(32442,-620997),
          new(40954,-614029),
          new(47244,-605763),
          new(50000,-595332),
          new(50438,-591472),
          new(55347,-581946),
          new(63107,-574604),
          new(72568,-568506),
          new(82553,-563595),
          new(92112,-559765),
          new(102721,-556206),
          new(114298,-552945),
          new(124199,-550567),
          new(134615,-548408),
          new(145495,-546477),
          new(156787,-544784),
          new(168436,-543337),
          new(180385,-542143),
          new(192576,-541209),
          new(204949,-540538),
          new(217444,-540135),
          new(230000,-540000),
          new(237542,-540000),
          new(241117,-539909),
          new(251801,-538540),
          new(262329,-535544),
          new(272585,-530954),
          new(282456,-524819),
          new(290312,-518575),
          new(297765,-511353),
          new(304760,-503206),
          new(311244,-494199),
          new(317167,-484398),
          new(322483,-473879),
          new(327154,-462721),
          new(331142,-451010),
          new(333821,-441303),
          new(336031,-431346),
          new(337761,-421187),
          new(339003,-410876),
          new(339750,-400463),
          new(340000,-390000),
          new(339909,-384765),
          new(339178,-374321),
          new(337721,-363953),
          new(335544,-353712),
          new(332658,-343647),
          new(329078,-333809),
          new(324819,-324244),
          new(319904,-315000),
          new(314356,-306121),
          new(308202,-297651),
          new(301472,-289630),
          new(294199,-282099),
          new(286418,-275093),
          new(278168,-268647),
          new(269488,-262793),
          new(260421,-257558),
          new(251010,-252968),
          new(241303,-249046),
          new(231346,-245811),
          new(221187,-243278),
          new(210876,-241460),
          new(200463,-240365),
          new(190000,-240000),
          new(189927,-240000),
          new(180166,-239562),
          new(168038,-237784),
          new(156076,-234653),
          new(146687,-231190),
          new(137509,-226893),
          new(128587,-221783),
          new(119964,-215885),
          new(111681,-209227),
          new(103779,-201842),
          new(96298,-193766),
          new(89272,-185039),
          new(82737,-175702),
          new(76724,-165801),
          new(71262,-155385),
          new(66379,-144505),
          new(62097,-133213),
          new(58439,-121564),
          new(55421,-109615),
          new(53058,-97424),
          new(51362,-85051),
          new(50341,-72556),
          new(50000,-60000),
          new(50000,-30000),
          new(49909,-24765),
          new(49178,-14321),
          new(47721,-3953),
          new(45544,6288),
          new(42658,16353),
          new(39078,26191),
          new(34819,35756),
          new(29904,45000),
          new(24356,53879),
          new(18202,62349),
          new(11472,70370),
          new(4199,77901),
          new(-3582,84907),
          new(-11832,91353),
          new(-20512,97207),
          new(-29579,102442),
          new(-38990,107032),
          new(-48697,110954),
          new(-58654,114189),
          new(-68813,116722),
          new(-79124,118540),
          new(-89537,119635),
          new(-100000,120000),
          new(-105235,119909),
          new(-115679,119178),
          new(-126047,117721),
          new(-136288,115544),
          new(-146353,112658),
          new(-156191,109078),
          new(-165756,104819),
          new(-175000,99904),
          new(-183879,94356),
          new(-192349,88202),
          new(-200370,81472),
          new(-207901,74199),
          new(-214907,66418),
          new(-221353,58168),
          new(-227207,49488),
          new(-232442,40421),
          new(-237032,31010),
          new(-240954,21303),
          new(-244189,11346),
          new(-246722,1187),
          new(-248540,-9124),
          new(-249635,-19537),
          new(-250000,-30000),
          new(-250000,-60000),
          new(-250062,-66282),
          new(-250555,-78815),
          new(-251539,-91257),
          new(-253010,-103546),
          new(-254959,-115623),
          new(-257378,-127429),
          new(-260255,-138907),
          new(-263575,-150000),
          new(-267323,-160655),
          new(-271480,-170819),
          new(-276026,-180444),
          new(-280939,-189481),
          new(-287560,-199886),
          new(-294665,-209227),
          new(-302202,-217432),
          new(-310113,-224438),
          new(-320015,-231190),
          new(-330260,-236067),
          new(-340735,-239014),
          new(-351326,-240000),
          new(-357014,-240206),
          new(-368327,-241847),
          new(-379453,-245111),
          new(-390272,-249963),
          new(-398966,-255181),
          new(-407297,-261425),
          new(-415203,-268647),
          new(-422623,-276794),
          new(-429500,-285801),
          new(-435782,-295602),
          new(-441421,-306121),
          new(-446374,-317279),
          new(-450605,-328990),
          new(-453446,-338697),
          new(-455790,-348654),
          new(-457625,-358813),
          new(-458942,-369124),
          new(-459735,-379537),
          new(-460000,-390000),
          };

        // Let's see if we can track the origin of the chords.
        for (int p = 0; p < a.Count; p++)
        {
            a[p] = new(a[p].X, a[p].Y, 1);
        }

        for (int p = 0; p < testPath.Count; p++)
        {
            testPath[p] = new(testPath[p].X, testPath[p].Y, 2);
        }

        Clipper2Lib.Clipper64 c = new();
        c.ZCallback = clipper2_zFillTest;
        c.AddOpenSubject(testPath);
        c.AddClip(a);

        Paths64 unused = new();
        Paths64 open = new();
        c.Execute(Clipper2Lib.ClipType.Intersection, Clipper2Lib.FillRule.EvenOdd, unused, open);

        Assert.AreEqual(2, open.Count);
        Assert.AreEqual(2, open[0].Count);
        Assert.AreEqual(2, open[1].Count);
    }

    [Test]
    public static void clipper2_rightChordTest()
    {
        Path64 testPath = new ()
        {
            new(-400000, -700000),
            new(-398907, -710396),
            new(-397815, -720791),
            new(-395106, -730902),
            new(-391355, -740674),
            new(-386603, -750000),
            new(-380902, -758779),
            new(-374314, -766913),
            new(-366913, -774314),
            new(-358779, -780902),
            new(-350000, -786603),
            new(-340674, -791355),
            new(-330902, -795106),
            new(-320791, -797815),
            new(-310453, -799452),
            new(-300000, -800000),
            new(-200000, -800000),
            new(-189547, -799452),
            new(-179209, -797815),
            new(-169098, -795106),
            new(-159326, -791355),
            new(-150000, -786603),
            new(-141221, -780902),
            new(-133087, -774314),
            new(-125686, -766913),
            new(-119098, -758779),
            new(-113397, -750000),
            new(-108645, -740674),
            new(-104894, -730902),
            new(-102185, -720791),
            new(-100548, -710453),
            new(-100000, -700000),
            new(-100000, -100000),
            new(-100548, -89547),
            new(-102185, -79209),
            new(-104894, -69098),
            new(-108645, -59326),
            new(-113397, -50000),
            new(-119098, -41221),
            new(-125686, -33087),
            new(-133087, -25686),
            new(-141221, -19098),
            new(-150000, -13397),
            new(-159326, -8645),
            new(-169098, -4894),
            new(-179209, -2185),
            new(-189547, -548),
            new(-200000, 0)
        };

        Path64 a = new ()
        {
            new(-460000, -390000),
            new(-459955, -395235),
            new(-459594, -405679),
            new(-458873, -416047),
            new(-457796, -426288),
            new(-456369, -436353),
            new(-454102, -448610),
            new(-451315, -460421),
            new(-448030, -471696),
            new(-444272, -482349),
            new(-440068, -492300),
            new(-435451, -501472),
            new(-429416, -511353),
            new(-422903, -519904),
            new(-414796, -528076),
            new(-406257, -534189),
            new(-396132, -538540),
            new(-385806, -540000),
            new(-381067, -540036),
            new(-369255, -540441),
            new(-357570, -541292),
            new(-346100, -542583),
            new(-334932, -544305),
            new(-324151, -546443),
            new(-313840, -548983),
            new(-304076, -551904),
            new(-293187, -555882),
            new(-283312, -560333),
            new(-273218, -566059),
            new(-264802, -572279),
            new(-257399, -579871),
            new(-252063, -588852),
            new(-250000, -599117),
            new(-242658, -614841),
            new(-234819, -621423),
            new(-225801, -626830),
            new(-216572, -631139),
            new(-206066, -635097),
            new(-196418, -638096),
            new(-186036, -640798),
            new(-175000, -643183),
            new(-163393, -645233),
            new(-151303, -646931),
            new(-141346, -648029),
            new(-131187, -648888),
            new(-120876, -649505),
            new(-110463, -649876),
            new(-100000, -650000),
            new(-94765, -649967),
            new(-84321, -649701),
            new(-73953, -649169),
            new(-63712, -648376),
            new(-53647, -647324),
            new(-41390, -645654),
            new(-29579, -643601),
            new(-18304, -641180),
            new(-7651, -638411),
            new(2300, -635314),
            new(13206, -631197),
            new(22873, -626688),
            new(32442, -620997),
            new(40954, -614029),
            new(47244, -605763),
            new(50000, -595332),
            new(50438, -591472),
            new(55347, -581946),
            new(63107, -574604),
            new(72568, -568506),
            new(82553, -563595),
            new(92112, -559765),
            new(102721, -556206),
            new(114298, -552945),
            new(124199, -550567),
            new(134615, -548408),
            new(145495, -546477),
            new(156787, -544784),
            new(168436, -543337),
            new(180385, -542143),
            new(192576, -541209),
            new(204949, -540538),
            new(217444, -540135),
            new(230000, -540000),
            new(237542, -540000),
            new(241117, -539909),
            new(251801, -538540),
            new(262329, -535544),
            new(272585, -530954),
            new(282456, -524819),
            new(290312, -518575),
            new(297765, -511353),
            new(304760, -503206),
            new(311244, -494199),
            new(317167, -484398),
            new(322483, -473879),
            new(327154, -462721),
            new(331142, -451010),
            new(333821, -441303),
            new(336031, -431346),
            new(337761, -421187),
            new(339003, -410876),
            new(339750, -400463),
            new(340000, -390000),
            new(339909, -384765),
            new(339178, -374321),
            new(337721, -363953),
            new(335544, -353712),
            new(332658, -343647),
            new(329078, -333809),
            new(324819, -324244),
            new(319904, -315000),
            new(314356, -306121),
            new(308202, -297651),
            new(301472, -289630),
            new(294199, -282099),
            new(286418, -275093),
            new(278168, -268647),
            new(269488, -262793),
            new(260421, -257558),
            new(251010, -252968),
            new(241303, -249046),
            new(231346, -245811),
            new(221187, -243278),
            new(210876, -241460),
            new(200463, -240365),
            new(190000, -240000),
            new(189927, -240000),
            new(180166, -239562),
            new(168038, -237784),
            new(156076, -234653),
            new(146687, -231190),
            new(137509, -226893),
            new(128587, -221783),
            new(119964, -215885),
            new(111681, -209227),
            new(103779, -201842),
            new(96298, -193766),
            new(89272, -185039),
            new(82737, -175702),
            new(76724, -165801),
            new(71262, -155385),
            new(66379, -144505),
            new(62097, -133213),
            new(58439, -121564),
            new(55421, -109615),
            new(53058, -97424),
            new(51362, -85051),
            new(50341, -72556),
            new(50000, -60000),
            new(50000, -30000),
            new(49909, -24765),
            new(49178, -14321),
            new(47721, -3953),
            new(45544, 6288),
            new(42658, 16353),
            new(39078, 26191),
            new(34819, 35756),
            new(29904, 45000),
            new(24356, 53879),
            new(18202, 62349),
            new(11472, 70370),
            new(4199, 77901),
            new(-3582, 84907),
            new(-11832, 91353),
            new(-20512, 97207),
            new(-29579, 102442),
            new(-38990, 107032),
            new(-48697, 110954),
            new(-58654, 114189),
            new(-68813, 116722),
            new(-79124, 118540),
            new(-89537, 119635),
            new(-100000, 120000),
            new(-105235, 119909),
            new(-115679, 119178),
            new(-126047, 117721),
            new(-136288, 115544),
            new(-146353, 112658),
            new(-156191, 109078),
            new(-165756, 104819),
            new(-175000, 99904),
            new(-183879, 94356),
            new(-192349, 88202),
            new(-200370, 81472),
            new(-207901, 74199),
            new(-214907, 66418),
            new(-221353, 58168),
            new(-227207, 49488),
            new(-232442, 40421),
            new(-237032, 31010),
            new(-240954, 21303),
            new(-244189, 11346),
            new(-246722, 1187),
            new(-248540, -9124),
            new(-249635, -19537),
            new(-250000, -30000),
            new(-250000, -60000),
            new(-250062, -66282),
            new(-250555, -78815),
            new(-251539, -91257),
            new(-253010, -103546),
            new(-254959, -115623),
            new(-257378, -127429),
            new(-260255, -138907),
            new(-263575, -150000),
            new(-267323, -160655),
            new(-271480, -170819),
            new(-276026, -180444),
            new(-280939, -189481),
            new(-287560, -199886),
            new(-294665, -209227),
            new(-302202, -217432),
            new(-310113, -224438),
            new(-320015, -231190),
            new(-330260, -236067),
            new(-340735, -239014),
            new(-351326, -240000),
            new(-357014, -240206),
            new(-368327, -241847),
            new(-379453, -245111),
            new(-390272, -249963),
            new(-398966, -255181),
            new(-407297, -261425),
            new(-415203, -268647),
            new(-422623, -276794),
            new(-429500, -285801),
            new(-435782, -295602),
            new(-441421, -306121),
            new(-446374, -317279),
            new(-450605, -328990),
            new(-453446, -338697),
            new(-455790, -348654),
            new(-457625, -358813),
            new(-458942, -369124),
            new(-459735, -379537),
            new(-460000, -390000),
        };

        // Let's see if we can track the origin of the chords.
        for (int p = 0; p < a.Count; p++)
        {
            a[p] = new(a[p].X, a[p].Y, 1);
        }

        for (int p = 0; p < testPath.Count; p++)
        {
            testPath[p] = new(testPath[p].X, testPath[p].Y, 2);
        }

        Clipper2Lib.Clipper64 c = new();
        c.ZCallback = clipper2_zFillTest;
        c.AddOpenSubject(testPath);
        c.AddClip(a);

        Paths64 unused = new();
        Paths64 open = new();
        c.Execute(Clipper2Lib.ClipType.Intersection, Clipper2Lib.FillRule.EvenOdd, unused, open);
        
        Assert.AreEqual(1, open.Count);
        Assert.AreEqual(17, open[0].Count);
    }
    
    [Test]
    public static void clipper1_clipper2_compare_S()
    {
        Clipper2Lib.Path64 BP = new()
        {
            new(1000, 27000),
            new(1000, 2000),
            new(27000, 2000),
            new(27000, 27000),
            new(1000, 27000)
        };

        Paths64 iPoly = new()
        {
            new()
            {
                new(1000, 2000),
                new(1000, 7000),
                new(6000, 7000),
                new(6000, 14000),
                new(1000, 14000),
                new(1000, 27000),
                new(27000, 27000),
                new(27000, 23000),
                new(16000, 23000),
                new(16000, 16000),
                new(27000, 16000),
                new(27000, 2000),
            }
        };

        List<ClipperLib1.IntPoint> BP1 = new();
        for (int pt = 0; pt < BP.Count; pt++)
        {
            BP1.Add(new (BP[pt].X, BP[pt].Y));
        }

        /*
        ClipperLib2.ClipperOffset co = new();
        co.AddPath(BP, ClipperLib2.JoinType.Miter, ClipperLib2.EndType.Polygon);
        BP = ClipperLib2.ClipperFunc.Paths(co.Execute(0.999))[0];
        */
        
        Paths iPoly1 = new();
        for (int p = 0; p < iPoly.Count; p++)
        {
            List<ClipperLib1.IntPoint> t = new();
            for (int pt = 0; pt < iPoly[p].Count; pt++)
            {
                t.Add(new (iPoly[p][pt].X, iPoly[p][pt].Y));
            }
            iPoly1.Add(t);
        }

        Clipper2Lib.ClipperOffset co = new() {PreserveCollinear = true, ReverseSolution = true};
        co.AddPaths(iPoly, Clipper2Lib.JoinType.Miter, Clipper2Lib.EndType.Polygon);
        co.Execute(1.0001, iPoly);
        
        double area = iPoly.Sum(t => Clipper2Lib.Clipper.Area(t));
        Assert.AreEqual(538134004, area);
        Assert.AreEqual(1, iPoly.Count);
        Assert.AreEqual(12, iPoly[0].Count);

        ClipperLib1.Clipper c1 = new() {PreserveCollinear = false};
        c1.AddPath(BP1, ClipperLib1.PolyType.ptSubject, true);
        c1.AddPaths(iPoly1, ClipperLib1.PolyType.ptClip, true);

        Paths o1 = new();
        c1.Execute(ClipperLib1.ClipType.ctDifference, o1);
        
        Clipper2Lib.Clipper64 c2 = new() {PreserveCollinear = false};
        c2.AddSubject(BP);
        c2.AddClip(iPoly);

        Paths64 o2 = new();
        c2.Execute(Clipper2Lib.ClipType.Difference, Clipper2Lib.FillRule.EvenOdd, o2);
        
        double areac1 = o1.Sum(t => ClipperLib1.Clipper.Area(t));
        Assert.AreEqual(112000000, areac1);
        Assert.AreEqual(2, o1.Count);
        Assert.AreEqual(4, o1[0].Count);
        Assert.AreEqual(4, o1[1].Count);
        
        double areac2 = o2.Sum(t => Clipper2Lib.Clipper.Area(t));
        Assert.AreEqual(111954004, areac2);
        Assert.AreEqual(2, o2.Count);
        Assert.AreEqual(4, o2[0].Count);
        Assert.AreEqual(4, o2[1].Count);
    }
    
    [Test]
    public static void clipper1_clipper2_compare_X()
    {
        Clipper2Lib.Paths64 xShape = new()
        {
            new()
            {
                new(0, 13000),
                new(0, 21000),
                new(13000, 21000),
                new(13000, 30000),
                new(24000, 30000),
                new(24000, 21000),
                new(48000, 21000),
                new(48000, 13000),
                new(24000, 13000),
                new(24000, 0),
                new(13000, 0),
                new(13000, 13000)
            }
        };

        Clipper2Lib.Path64 bounds = new()
        {
            new (0,30000),
            new (0,0),
            new (48000,0),
            new (48000,30000),
            new (0,30000)
        };

        Clipper2Lib.Clipper64 c2 = new();
        c2.AddSubject( bounds );
        c2.AddClip( xShape);
        Clipper2Lib.Paths64 c2out = new();
        c2.Execute(Clipper2Lib.ClipType.Difference, Clipper2Lib.FillRule.EvenOdd, c2out);
        
        double area = c2out.Sum(t => Clipper2Lib.Clipper.Area(t));
        Assert.AreEqual(814000000, area);
        Assert.AreEqual(4, c2out.Count);
        Assert.AreEqual(4, c2out[0].Count);
        Assert.AreEqual(4, c2out[1].Count);
        Assert.AreEqual(4, c2out[2].Count);
        Assert.AreEqual(4, c2out[3].Count);

        List<List<ClipperLib1.IntPoint>> xShape1 = new();
        for (int p = 0; p < xShape.Count; p++)
        {
            List<ClipperLib1.IntPoint> t = new();
            for (int pt = 0; pt < xShape[p].Count; pt++)
            {
                t.Add(new (xShape[p][pt].X, xShape[p][pt].Y));
            }
            xShape1.Add(t);
        }

        List<ClipperLib1.IntPoint> bounds1 = new();
        for (int pt = 0; pt < bounds.Count; pt++)
        {
            bounds1.Add(new (bounds[pt].X, bounds[pt].Y));
        }

        ClipperLib1.Clipper c1 = new() {PreserveCollinear = false};
        c1.AddPath(bounds1, ClipperLib1.PolyType.ptSubject, true);
        c1.AddPaths(xShape1, ClipperLib1.PolyType.ptClip, true);

        List<List<ClipperLib1.IntPoint>> o1 = new();
        c1.Execute(ClipperLib1.ClipType.ctDifference, o1);

        double area2 = o1.Sum(t => ClipperLib1.Clipper.Area(t));
        Assert.AreEqual(814000000, area2);
        Assert.AreEqual(4, o1.Count);
        Assert.AreEqual(4, o1[0].Count);
        Assert.AreEqual(4, o1[1].Count);
        Assert.AreEqual(4, o1[2].Count);
        Assert.AreEqual(4, o1[3].Count);
    }
    
    [Test]
    public static void clipper1_clipper2_compare_subtraction()
    {
        Paths64 layerAPaths = new()
        {
            new()
            {
                new(0, 0),
                new(0, 900000),
                new(500000, 900000),
                new(500000, 0),
                new(0, 0),
            }
        };


        Paths64 layerBPaths = new()
        {
            new()
            {
                new(0, 0),
                new(0, 11900),
                new(0, 23800),
                new(0, 35700),
                new(0, 47600),
                new(0, 59500),
                new(12500, 59500),
                new(25000, 59500),
                new(37500, 59500),
                new(50000, 59500),
                new(50000, 50000),
                new(350000, 50000),
                new(350000, 350000),
                new(50000, 350000),
                new(50000, 339661),
                new(50000, 329321),
                new(50000, 318982),
                new(50000, 308643),
                new(50000, 298304),
                new(50000, 287964),
                new(50000, 277625),
                new(50000, 267286),
                new(50000, 256946),
                new(50000, 246607),
                new(50000, 236268),
                new(50000, 225929),
                new(50000, 215589),
                new(50000, 205250),
                new(50000, 194911),
                new(50000, 184571),
                new(50000, 174232),
                new(50000, 163893),
                new(50000, 153554),
                new(50000, 143214),
                new(50000, 132875),
                new(50000, 122536),
                new(50000, 112196),
                new(50000, 101857),
                new(50000, 91518),
                new(50000, 81179),
                new(50000, 70839),
                new(50000, 60500),
                new(37500, 60500),
                new(25000, 60500),
                new(12500, 60500),
                new(0, 60500),
                new(0, 70731),
                new(0, 80962),
                new(0, 91192),
                new(0, 101423),
                new(0, 111654),
                new(0, 121885),
                new(0, 132115),
                new(0, 142346),
                new(0, 152577),
                new(0, 162808),
                new(0, 173038),
                new(0, 183269),
                new(0, 193500),
                new(0, 203731),
                new(0, 213962),
                new(0, 224192),
                new(0, 234423),
                new(0, 244654),
                new(0, 254885),
                new(0, 265115),
                new(0, 275346),
                new(0, 285577),
                new(0, 295808),
                new(0, 306038),
                new(0, 316269),
                new(0, 326500),
                new(0, 336731),
                new(0, 346962),
                new(0, 357192),
                new(0, 367423),
                new(0, 377654),
                new(0, 387885),
                new(0, 398115),
                new(0, 408346),
                new(0, 418577),
                new(0, 428808),
                new(0, 439038),
                new(0, 449269),
                new(0, 459500),
                new(12500, 459500),
                new(25000, 459500),
                new(37500, 459500),
                new(50000, 459500),
                new(50000, 450000),
                new(350000, 450000),
                new(350000, 750000),
                new(50000, 750000),
                new(50000, 739661),
                new(50000, 729321),
                new(50000, 718982),
                new(50000, 708643),
                new(50000, 698304),
                new(50000, 687964),
                new(50000, 677625),
                new(50000, 667286),
                new(50000, 656946),
                new(50000, 646607),
                new(50000, 636268),
                new(50000, 625929),
                new(50000, 615589),
                new(50000, 605250),
                new(50000, 594911),
                new(50000, 584571),
                new(50000, 574232),
                new(50000, 563893),
                new(50000, 553554),
                new(50000, 543214),
                new(50000, 532875),
                new(50000, 522536),
                new(50000, 512196),
                new(50000, 501857),
                new(50000, 491518),
                new(50000, 481179),
                new(50000, 470839),
                new(50000, 460500),
                new(37500, 460500),
                new(25000, 460500),
                new(12500, 460500),
                new(0, 460500),
                new(0, 470721),
                new(0, 480942),
                new(0, 491163),
                new(0, 501384),
                new(0, 511605),
                new(0, 521826),
                new(0, 532047),
                new(0, 542267),
                new(0, 552488),
                new(0, 562709),
                new(0, 572930),
                new(0, 583151),
                new(0, 593372),
                new(0, 603593),
                new(0, 613814),
                new(0, 624035),
                new(0, 634256),
                new(0, 644477),
                new(0, 654698),
                new(0, 664919),
                new(0, 675140),
                new(0, 685360),
                new(0, 695581),
                new(0, 705802),
                new(0, 716023),
                new(0, 726244),
                new(0, 736465),
                new(0, 746686),
                new(0, 756907),
                new(0, 767128),
                new(0, 777349),
                new(0, 787570),
                new(0, 797791),
                new(0, 808012),
                new(0, 818233),
                new(0, 828453),
                new(0, 838674),
                new(0, 848895),
                new(0, 859116),
                new(0, 869337),
                new(0, 879558),
                new(0, 889779),
                new(0, 900000),
                new(500000, 900000),
                new(500000, 0),
                new(0, 0),
                new(-2147483647, -2147483647),
                new(-2147483647, 2147483647),
                new(2147483647, 2147483647),
                new(2147483647, -2147483647),
                new(-2147483647, -2147483647),
                new(0, 0),
            }
        };

        Clipper2Lib.Clipper64 c = new()
        {
            PreserveCollinear = true
        };

        c.AddSubject(layerBPaths);
        c.AddClip(layerAPaths);

        Paths64 solution_cl = new();

        c.Execute(Clipper2Lib.ClipType.Difference, Clipper2Lib.FillRule.EvenOdd, solution_cl);

        double area = solution_cl.Sum(t => Clipper2Lib.Clipper.Area(t));
        Assert.AreEqual(1.8446743606529683E+19d, area);
        Assert.AreEqual(2, solution_cl.Count);
        Assert.AreEqual(4, solution_cl[0].Count);
        Assert.AreEqual(4, solution_cl[1].Count);

        c.PreserveCollinear = false;

        Paths64 solution_ncl = new();

        c.Execute(Clipper2Lib.ClipType.Difference, Clipper2Lib.FillRule.EvenOdd, solution_ncl);

        double area2 = solution_ncl.Sum(t => Clipper2Lib.Clipper.Area(t));
        Assert.AreEqual(1.8446743606529683E+19, area2);
        Assert.AreEqual(2, solution_ncl.Count);
        Assert.AreEqual(4, solution_ncl[0].Count);
        Assert.AreEqual(4, solution_ncl[1].Count);

        ClipperLib1.Clipper c1 = new()
        {
            PreserveCollinear = true
        };

        Paths layerBPaths1 = new();

        foreach (List<Clipper2Lib.Point64> t in layerBPaths)
        {
            List<ClipperLib1.IntPoint> r = new();
            for (int pt = 0; pt < t.Count; pt++)
            {
                r.Add(new(t[pt].X, t[pt].Y));
            }

            layerBPaths1.Add(r);
        }

        Paths layerAPaths1 = new();
        foreach (List<Clipper2Lib.Point64> t in layerAPaths)
        {
            List<ClipperLib1.IntPoint> r = new();
            for (int pt = 0; pt < t.Count; pt++)
            {
                r.Add(new(t[pt].X, t[pt].Y));
            }

            layerAPaths1.Add(r);
        }

        c1.AddPaths(layerBPaths1, ClipperLib1.PolyType.ptSubject, true);
        c1.AddPaths(layerAPaths1, ClipperLib1.PolyType.ptClip, true);

        Paths solution_cl1 = new();

        c1.Execute(ClipperLib1.ClipType.ctDifference, solution_cl1);

        double area3 = solution_cl1.Sum(t => ClipperLib1.Clipper.Area(t));
        Assert.AreEqual(1.8446743606529683E+19d, area3);
        Assert.AreEqual(2, solution_cl1.Count);
        Assert.AreEqual(4, solution_cl1[0].Count);
        Assert.AreEqual(4, solution_cl1[1].Count);

        c1.PreserveCollinear = false;

        Paths solution_ncl1 = new();

        c1.Execute(ClipperLib1.ClipType.ctDifference, solution_ncl1);

        double area4 = solution_ncl1.Sum(t => ClipperLib1.Clipper.Area(t));
        Assert.AreEqual(1.8446743606529683E+19d, area4);
        Assert.AreEqual(2, solution_ncl1.Count);
        Assert.AreEqual(4, solution_ncl1[0].Count);
        Assert.AreEqual(4, solution_ncl1[1].Count);
    }
}