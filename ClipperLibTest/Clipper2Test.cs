﻿
using ClipperLib2;

namespace ClipperLib2Test;

using Path64 = List<Point64>;
using Paths64 = List<List<Point64>>;

public static class Clipper2Test
{
    const double keyhole_sizing = 500;
    
    public static void test2()
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
        c.AddSubject(rays, true);
        c.AddClip(collisionPaths);
        PolyTree pt = new PolyTree();
        Paths64 solution = new();
        c.Execute(ClipType.Intersection, FillRule.EvenOdd, pt, solution);
    }
    
    public static void test1()
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