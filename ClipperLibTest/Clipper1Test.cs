﻿
using ClipperLib1;

namespace ClipperLib1Test;

using Path = List<IntPoint>;
using Paths = List<List<IntPoint>>;

public static class Clipper1Test
{
    const double keyhole_sizing = 500;
    public static void test1()
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
        
        ClipperOffset co = new();
        co.AddPaths(kHSource, JoinType.jtMiter, EndType.etClosedPolygon);
        Paths out_ = new ();
        co.Execute(ref out_, keyhole_sizing);

        Console.WriteLine("Out count: " + out_.Count);
        
    }

    public static void test2()
    {
        Paths rays = new()
        {
            new Path() {new IntPoint(100000, 200000), new IntPoint(100000, -9800000)}
        };

        Paths collisionPaths = new()
        {
            new Path()
            {
                new IntPoint(0, 0),
                new IntPoint(0, 500000),
                new IntPoint(100000, 500000),
                new IntPoint(100000, 200000),
                new IntPoint(600000, 200000),
                new IntPoint(600000, 800000),
                new IntPoint(1200000, 800000),
                new IntPoint(1200000, 0),
                new IntPoint(0, 0)

            }
        };

        Clipper c = new Clipper();
        c.AddPaths(rays, PolyType.ptSubject, false);
        c.AddPaths(collisionPaths, PolyType.ptClip, true);
        PolyTree pt = new PolyTree();
        c.Execute(ClipType.ctIntersection, pt);
        Paths solution = Clipper.OpenPathsFromPolyTree(pt);
    }
    public static void offsetTest()
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

        ClipperOffset co = new ClipperOffset();
        co.AddPaths(newEdges, JoinType.jtMiter, EndType.etOpenSquare);
        PolyTree tp = new();
        co.Execute(ref tp, 1.0);

        Paths cutters = Clipper.ClosedPathsFromPolyTree(tp);

        Clipper c = new Clipper();
        c.AddPath(lPoly, PolyType.ptSubject, true);
        c.AddPaths(cutters, PolyType.ptClip, true);
        Paths solution = new();
        c.Execute(ClipType.ctDifference, solution);

    }
    
}