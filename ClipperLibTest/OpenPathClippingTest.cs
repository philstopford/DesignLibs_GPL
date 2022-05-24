using Clipper2Lib;

namespace ClipperLibTest;

using Path64 = List<Point64>;
using Paths64 = List<List<Point64>>;

public class OpenPathClippingTest
{
    public static void test()
    {
        Path64 testPath = ClipperFunc.MakePath(new[]
        {
            -50000, -550000,
            -50000, -150000,
            650000, -150000
        });
        
        Path64 b = ClipperFunc.MakePath(new[]
        {
            300000,-800000,
            300000,0,
            500000,0,
            500000,-800000
        });
        
        Clipper c = new() {PreserveCollinear = true};
        c.AddOpenSubject(testPath);
        c.AddClip(b);
        Paths64 unused = new();
        Paths64 topChords = new();
        c.Execute(ClipType.Intersection, FillRule.EvenOdd, unused, topChords);

        Path64 testPath2 = ClipperFunc.MakePath(new[]
        {
            650000,-150000,
            650000,-550000,
            -50000,-550000
        });
        
        c.Clear();
        c.AddOpenSubject(testPath2);
        c.AddClip(b);
        Paths64 bottomChords = new();
        c.Execute(ClipType.Intersection, FillRule.EvenOdd, unused, bottomChords);
        
        Path64 testPath3 = ClipperFunc.MakePath(new[]
        {
            300000,-800000,
            300000,0
        });

        Path64 a = ClipperFunc.MakePath(new[]
        {
            -50000, -550000,
            -50000, -150000,
            650000, -150000,
            650000, -550000
        });

        c.Clear();
        c.AddOpenSubject(testPath3);
        c.AddClip(a);
        Paths64 leftChords = new();
        c.Execute(ClipType.Intersection, FillRule.EvenOdd, unused, leftChords);

        Path64 testPath4 = ClipperFunc.MakePath(new[]
        {
            300000,0,
            500000,0,
            500000,-800000,
            300000,-800000
        });

        c.Clear();
        c.AddOpenSubject(testPath4);
        c.AddClip(a);
        Paths64 rightChords = new();
        c.Execute(ClipType.Intersection, FillRule.EvenOdd, unused, rightChords);

    }
}