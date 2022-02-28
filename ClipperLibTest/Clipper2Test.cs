
using ClipperLib2;

namespace ClipperLib2Test;

using Path = List<Point64>;
using Paths = List<List<Point64>>;

public static class Clipper2Test
{
    const double keyhole_sizing = 500;
    public static void test1()
    {
        Console.WriteLine("Clipper2 Test1");
        Path outer = new()
        {
            new Point64(-200000, -200000),
            new Point64(200000, -200000),
            new Point64(200000, 200000),
            new Point64(-200000, 200000),
            new Point64(-200000, -200000)
        };

        Path inner1 = new()
        {
            new Point64(-100000, -100000),
            new Point64(-100000, 100000),
            new Point64(100000, 100000),
            new Point64(100000, -100000),
            new Point64(-100000, -100000)
        };

        Paths kHSource = new()
        {
            outer,
            inner1
        };
        
        ClipperOffset co = new();
        co.AddPaths(kHSource, JoinType.Miter, EndType.Polygon);
        Paths out_ = ClipperFunc.PathsFromPathsD(co.Execute(keyhole_sizing));
        
        Console.WriteLine("Out count: " + out_.Count);

    }
    
    public static void offsetTest()
    {
        Path lPoly = new ()
        {
            new Point64(0, 0),
            new Point64(0, 500000),
            new Point64(100000, 500000),
            new Point64(100000, 200000),
            new Point64(600000, 200000),
            new Point64(600000, 0),
            new Point64(0, 0)
        };

        Path newEdge = new()
        {
            new Point64(100000, 200000),
            new Point64(100000, 0)
        };

        Paths newEdges = new()
        {
            newEdge
        };

        ClipperOffset co = new ClipperOffset();
        co.AddPaths(newEdges, JoinType.Miter, EndType.Square);
        PolyTree tp = new();
        Paths cutters = ClipperFunc.PathsFromPathsD(co.Execute( 2));
        
        Clipper c = new Clipper();
        c.AddSubject(lPoly);
        c.AddClip(cutters);
        c.Execute(ClipType.Difference, FillRule.EvenOdd, tp);
        Paths solution = ClipperFunc.PolyTreeToPaths(tp);

    }
}