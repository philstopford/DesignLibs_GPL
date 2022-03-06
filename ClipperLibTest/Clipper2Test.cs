
using ClipperLib2;

namespace ClipperLib2Test;

using Path64 = List<Point64>;
using Paths64 = List<List<Point64>>;

public static class Clipper2Test
{
    const double keyhole_sizing = 500;

    private static void zFillTest(Point64 bot1, Point64 top1, Point64 bot2, Point64 top2, ref Point64 pt)
    {
        pt.Z = -1;
    }
    public static void zFillCallbackTest()
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
            new(-100, -11000),
            new(-100, 100),
            new(100, 100),
            new(100, -1100),
        };

        Clipper c = new();

        c.ZFill = zFillTest;

        c.AddSubject(outer);
        c.AddClip(cutter);

        Paths64 solution = new();
        c.Execute(ClipType.Intersection, FillRule.EvenOdd, solution);
    }
    public static void coincident_openPathTest()
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

        Path64 t = new()
        {
            new Point64(200001, 1),
            new Point64(200000, 0),
            new Point64(200000, 2000000),
            new Point64(200001, 1100000-1),
        };

        Clipper c = new();
        c.AddClip(lPoly);
        c.AddOpenSubject(t);

        Paths64 open = new();
        Paths64 solution = new();
        c.Execute(ClipType.Intersection, FillRule.EvenOdd, solution, open);
        
        Path64 t2 = new()
        {
            new Point64(200000, 2000000),
            new Point64(200000, 0),
            new Point64(200001, 0),
        };

        Clipper c2 = new();
        c2.AddClip(lPoly);
        c2.AddOpenSubject(t2);

        Paths64 open2 = new();
        Paths64 solution2 = new();
        c2.Execute(ClipType.Intersection, FillRule.EvenOdd, solution2, open2);
        
        Path64 t3 = new()
        {
            new Point64(200000, 0),
            new Point64(900000, 0),
        };

        Clipper c3 = new();
        c3.AddClip(lPoly);
        c3.AddOpenSubject(t3);

        Paths64 open3 = new();
        Paths64 solution3 = new();
        c3.Execute(ClipType.Intersection, FillRule.EvenOdd, solution3, open3);
    }
    
    public static void keyHole_test2()
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
                new Point64(800000, 800000),
                new Point64(800000, 1100000)
            }
        };
        
        // Turn the new edges into cutters and slice. Not terribly elegant and we're relying on rounding to squash notches later.
        ClipperOffset co = new();
        co.AddPaths(t, JoinType.Miter, EndType.Square);

        Paths64 cutters = ClipperFunc.Paths(co.Execute(2.0));

        Clipper c = new();

        c.AddSubject(lPoly);

        // Take first cutter only - we only cut once, no matter how many potential cutters we have.
        c.AddClip(cutters[0]);
        Paths64 f = new();
        c.Execute(ClipType.Difference, FillRule.EvenOdd, f);
    }
    public static void openPath_clipTest1()
    {
        Path64 lPoly = new Path64() {
        new Point64(0,0),
        new Point64(0,200000),
        new Point64(200000,200000),
        new Point64(200000,500000),
        new Point64(0,500000),
        new Point64(0,1100000),
        new Point64(1000000,1100000),
        new Point64(1000000,800000),
        new Point64(800000,800000),
        new Point64(800000,600000),
        new Point64(1000000,600000),
        new Point64(1000000,0),
        new Point64(0,0)
        };

        Path64 t = new () {
        new Point64(0,200000),
        new Point64(0,-9800000)
        };

        Clipper c = new();

        c.AddOpenSubject(t);
        c.AddClip(lPoly);

        PolyTree pt = new();
        Paths64 p = new();

        c.Execute(ClipType.Intersection, FillRule.EvenOdd, pt, p);
 
    }
    
    public static void openPath_clipTest2()
    {
        Paths64 rays2 = new()
        {
            new Path64() {new Point64(100000, 200000), new Point64(100000, -9800000)}
        };

        Paths64 rays = new()
        {
            new Path64() {new Point64(100000, 500000), new Point64(10100000, 500000)}
        };

        Paths64 collisionPaths = new()
        {
            new Path64()
            {
                new Point64(0, 0),
                new Point64(0, 500000),
                new Point64(100000, 500000),
                new Point64(100000, 200000),
                new Point64(600000, 200000),
                new Point64(600000, 800000),
                new Point64(1200000, 800000),
                new Point64(1200000, 0),
                new Point64(0, 0)

            }
        };

        Clipper c = new Clipper();
        c.AddOpenSubject(rays);
        c.AddClip(collisionPaths);
        PolyTree pt = new PolyTree();
        Paths64 solution = new();
        c.Execute(ClipType.Intersection, FillRule.EvenOdd, pt, solution);
    }
    
    public static void keyHole_test1()
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
        
        ClipperOffset co = new();
        co.AddPaths(kHSource, JoinType.Miter, EndType.Polygon);
        // ClipperLib2 specifies full width in offset for open path, unlike version 1
        Paths64 out_ = ClipperFunc.Paths(co.Execute(2*keyhole_sizing));
        
        Console.WriteLine("Out count: " + out_.Count);

    }
    
    public static void offsetTest()
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

        ClipperOffset co = new ClipperOffset();
        co.AddPaths(newEdges, JoinType.Miter, EndType.Square);
        PolyTree tp = new();
        // ClipperLib2 specifies full width in offset for open path, unlike version 1
        Paths64 cutters = ClipperFunc.Paths(co.Execute( 2.0));
        
        Clipper c = new Clipper();
        c.AddSubject(lPoly);
        c.AddClip(cutters);
        c.Execute(ClipType.Difference, FillRule.EvenOdd, tp);
        Paths64 solution = ClipperFunc.PolyTreeToPaths(tp);

    }
}