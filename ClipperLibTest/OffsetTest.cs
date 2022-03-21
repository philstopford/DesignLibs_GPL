namespace ClipperLibTest;

using Path64 = List<Clipper2Lib.Point64>;
using Paths64 = List<List<Clipper2Lib.Point64>>;
using Path = List<ClipperLib1.IntPoint>;
using Paths = List<List<ClipperLib1.IntPoint>>;

public static class OffsetTest
{
    private static Path64 test = new()
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
    
    public static void compare()
    {
        Path test1 = new();
        for (int pt = 0; pt < test.Count; pt++)
        {
            test1.Add(new ClipperLib1.IntPoint(test[pt].X, test[pt].Y));
        }
        ClipperLib1.ClipperOffset co1 = new();
        co1.AddPath(test1, ClipperLib1.JoinType.jtMiter, ClipperLib1.EndType.etClosedPolygon);
        Paths c1up = new();
        co1.Execute(ref c1up, 2.0);
        co1.Clear();
        co1.AddPaths(c1up, ClipperLib1.JoinType.jtMiter, ClipperLib1.EndType.etClosedPolygon);
        Paths c1down = new();
        co1.Execute(ref c1down, -2.0);
        
        Clipper2Lib.ClipperOffset co2 = new();
        co2.AddPath(test, Clipper2Lib.JoinType.Miter, Clipper2Lib.EndType.Polygon);
        Paths64 c2up = Clipper2Lib.ClipperFunc.Paths64(co2.Execute(2.0));
        co2.Clear();
        co2.AddPaths(c2up, Clipper2Lib.JoinType.Miter, Clipper2Lib.EndType.Polygon);
        Paths64 c2down = Clipper2Lib.ClipperFunc.Paths64(co2.Execute(-2.0));
    }
}