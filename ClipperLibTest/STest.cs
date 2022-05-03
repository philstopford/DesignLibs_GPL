namespace ClipperLibTest;

using Paths = List<List<ClipperLib1.IntPoint>>;

public class STest
{
    public static void compare()
    {
        Clipper2Lib.Path64 BP = new()
        {
            new(1000, 27000),
            new(1000, 2000),
            new(27000, 2000),
            new(27000, 27000),
            new(1000, 27000)
        };

        Clipper2Lib.Paths64 iPoly = new()
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

        Clipper2Lib.ClipperOffset co = new();
        co.AddPaths(iPoly, Clipper2Lib.JoinType.Miter, Clipper2Lib.EndType.Polygon);
        iPoly = Clipper2Lib.ClipperFunc.Paths64(co.Execute(1.0001));

        ClipperLib1.Clipper c1 = new() {PreserveCollinear = false};
        c1.AddPath(BP1, ClipperLib1.PolyType.ptSubject, true);
        c1.AddPaths(iPoly1, ClipperLib1.PolyType.ptClip, true);

        Paths o1 = new();
        c1.Execute(ClipperLib1.ClipType.ctDifference, o1);
        
        Clipper2Lib.Clipper c2 = new() {PreserveCollinear = false};
        c2.AddSubject(BP);
        c2.AddClip(iPoly);

        Clipper2Lib.Paths64 o2 = new();
        c2.Execute(Clipper2Lib.ClipType.Difference, Clipper2Lib.FillRule.EvenOdd, o2);

    }
}