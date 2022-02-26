
using ClipperLib2;

namespace ClipperLib2Test;

using Path = List<Point64>;
using Paths = List<List<Point64>>;

public static class Clipper2Test
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
            if (!ClipperFunc.Orientation(t))
            {
                r = (int)type.cutter;
            }

            ret[r].Add(new Path(t));
        }

        return ret;
    }
    public static void test1()
    {
        Console.WriteLine("Clipper2 Test1");
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
            outer.Add(new Point64(x, y));
        }
        for (; y <= 200000; y += 10000)
        {
            outer.Add(new Point64(x, y));
        }
        for (; x >= -200000; x -= 10000)
        {
            outer.Add(new Point64(x, y));
        }
        for (; y >= -200000; y -= 10000)
        {
            outer.Add(new Point64(x, y));
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
            inner1.Add(new Point64(x, y));
        }
        for (; x <= 100000; x += 10000)
        {
            inner1.Add(new Point64(x, y));
        }
        for (; y >= -100000; y -= 10000)
        {
            inner1.Add(new Point64(x, y));
        }
        for (; x >= -100000; x -= 10000)
        {
            inner1.Add(new Point64(x, y));
        }
     
        Paths kHSource = new()
        {
            outer,
            inner1
        };
        
        ClipperOffset co = new();
        co.AddPaths(kHSource, JoinType.Miter, EndType.Closed);
        Paths out_ = ClipperFunc.PathsFromPathsD(co.Execute(keyhole_sizing));
        
        Console.WriteLine("Out count: " + out_.Count);

    }
}