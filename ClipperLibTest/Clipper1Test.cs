
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
}