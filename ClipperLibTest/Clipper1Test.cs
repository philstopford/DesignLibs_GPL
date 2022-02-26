
using ClipperLib1;
using PolyTree = ClipperLib2.PolyTree;

namespace ClipperLib1Test;

using Path = List<IntPoint>;
using Paths = List<List<IntPoint>>;

public static class Clipper1Test
{
    const double keyhole_sizing = 500;
    public enum type { outer, cutter }

    private static Paths[] pGetDecomposed(Paths source)
    {
        Paths[] ret = new Paths[2];
        ret[0] = new Paths();
        ret[1] = new Paths();

        foreach (Path t in source)
        {
            int r = (int)type.outer;
            if (!Clipper.Orientation(t))
            {
                r = (int)type.cutter;
            }

            ret[r].Add(new Path(t));
        }

        return ret;
    }

    public static void test1()
    {
        Console.WriteLine("Clipper1 Test1");
        Path outer = new();
        /*
        {
            new IntPoint(-200000, -200000),
            new IntPoint(200000, -200000),
            new IntPoint(200000, 200000),
            new IntPoint(-200000, 200000),
            new IntPoint(-200000, -200000)
        };
        */
        int x = -200000;
        int y = -200000;
        for (x = -200000; x <= 200000; x += 10000)
        {
            outer.Add(new IntPoint(x, y));
        }
        for (; y <= 200000; y += 10000)
        {
            outer.Add(new IntPoint(x, y));
        }
        for (; x >= -200000; x -= 10000)
        {
            outer.Add(new IntPoint(x, y));
        }
        for (; y >= -200000; y -= 10000)
        {
            outer.Add(new IntPoint(x, y));
        }

        Path inner1 = new();
        /*
        {
            new IntPoint(-100000, -100000),
            new IntPoint(-100000, 100000),
            new IntPoint(100000, 100000),
            new IntPoint(100000, -100000),
            new IntPoint(-100000, -100000)
        };
        */        
        x = -100000;
        y = -100000;
        for (; y <= 100000; y += 10000)
        {
            inner1.Add(new IntPoint(x, y));
        }
        for (; x <= 100000; x += 10000)
        {
            inner1.Add(new IntPoint(x, y));
        }
        for (; y >= -100000; y -= 10000)
        {
            inner1.Add(new IntPoint(x, y));
        }
        for (; x >= -100000; x -= 10000)
        {
            inner1.Add(new IntPoint(x, y));
        }
        
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