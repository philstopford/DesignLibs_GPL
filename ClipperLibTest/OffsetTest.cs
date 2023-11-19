using Clipper2Lib;
using NUnit.Framework;

namespace ClipperLibTest;

using Path = List<ClipperLib1.IntPoint>;
using Paths = List<List<ClipperLib1.IntPoint>>;

public static class OffsetTest
{
    public static void test3()
    {
        int[] pointData = new[]
        {
            -360,-100,
            120,260,
            300,-140,
            -340,160
        }

;

        Path64 sourcePath = Clipper.MakePath(pointData);
        
        ClipperOffset co = new(miterLimit:2, arcTolerance:0.25);
        co.AddPath(sourcePath, JoinType.Miter, EndType.Square);
        Paths64 resizedPolyData = new();
        co.Execute(50, resizedPolyData);

        double area = resizedPolyData.Sum(t => Clipper.Area(t));
        Assert.AreEqual(1, area);
        Assert.AreEqual(1, resizedPolyData.Count);
        Assert.AreEqual(14, resizedPolyData[0].Count);
    }


    public static void test2()
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

        Path64 sourcePath = Clipper.MakePath(pointData);
        
        ClipperOffset co = new() {PreserveCollinear = true, ReverseSolution = true};
        co.AddPath(sourcePath, JoinType.Miter, EndType.Polygon);
        Paths64 resizedPolyData = new();
        co.Execute(Convert.ToDouble(6 * 10000), resizedPolyData);
        
        double area = resizedPolyData.Sum(t => Clipper.Area(t));
        Assert.AreEqual(1, area);
        Assert.AreEqual(1, resizedPolyData.Count);
        Assert.AreEqual(14, resizedPolyData[0].Count);
    }
    
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
        
        double area = c1up.Sum(t => ClipperLib1.Clipper.Area(t));
        Assert.AreEqual(1, area);
        Assert.AreEqual(1, c1up.Count);
        Assert.AreEqual(14, c1up[0].Count);
        
        co1.Clear();
        co1.AddPaths(c1up, ClipperLib1.JoinType.jtMiter, ClipperLib1.EndType.etClosedPolygon);
        Paths c1down = new();
        co1.Execute(ref c1down, -2.0);

        double area1 = c1down.Sum(t => ClipperLib1.Clipper.Area(t));
        Assert.AreEqual(1, area1);
        Assert.AreEqual(1, c1down.Count);
        Assert.AreEqual(14, c1down[0].Count);
        
        ClipperOffset co2 = new() {PreserveCollinear = true, ReverseSolution = true};
        co2.AddPath(test, Clipper2Lib.JoinType.Miter, Clipper2Lib.EndType.Polygon);
        Paths64 c2up = new();
        co2.Execute(2.0, c2up);
        
        double area2 = c2up.Sum(t => Clipper.Area(t));
        Assert.AreEqual(1, area2);
        Assert.AreEqual(1, c2up.Count);
        Assert.AreEqual(14, c2up[0].Count);
        
        co2.Clear();
        co2.AddPaths(c2up, Clipper2Lib.JoinType.Miter, Clipper2Lib.EndType.Polygon);
        Paths64 c2down = new();
        co2.Execute(-2.0, c2down);
        
        double area3 = c2down.Sum(t => Clipper.Area(t));
        Assert.AreEqual(1, area3);
        Assert.AreEqual(1, c2down.Count);
        Assert.AreEqual(14, c2down[0].Count);
    }
}