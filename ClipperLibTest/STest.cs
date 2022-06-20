using System.Collections.Generic;

namespace ClipperLibTest;

using Paths64 = List<List<Clipper2Lib.Point64>>;
using Paths = List<List<ClipperLib1.IntPoint>>;

public class STest
{
    public static void compare()
    {
        List<Clipper2Lib.Point64> BP = new()
        {
            new(1000, 27000),
            new(1000, 2000),
            new(27000, 2000),
            new(27000, 27000),
            new(1000, 27000)
        };

        Paths64 iPoly = new()
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

        Clipper2Lib.ClipperOffset co = new() {PreserveCollinear = true, ReverseSolution = true};
        co.AddPaths(iPoly, Clipper2Lib.JoinType.Miter, Clipper2Lib.EndType.Polygon);
        iPoly = co.Execute(1.0001);
        
        /* Expected output
        iPoly = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 12
          [0] = {Point64} 999,1999,0 
          [1] = {Point64} 999,7001,0 
          [2] = {Point64} 5999,7001,0 
          [3] = {Point64} 5999,13999,0 
          [4] = {Point64} 999,13999,0 
          [5] = {Point64} 999,27001,0 
          [6] = {Point64} 27001,27001,0 
          [7] = {Point64} 27001,22999,0 
          [8] = {Point64} 16001,22999,0 
          [9] = {Point64} 16001,16001,0 
          [10] = {Point64} 27001,16001,0 
          [11] = {Point64} 27001,1999,0 
           */

        ClipperLib1.Clipper c1 = new() {PreserveCollinear = false};
        c1.AddPath(BP1, ClipperLib1.PolyType.ptSubject, true);
        c1.AddPaths(iPoly1, ClipperLib1.PolyType.ptClip, true);

        Paths o1 = new();
        c1.Execute(ClipperLib1.ClipType.ctDifference, o1);
        
        Clipper2Lib.Clipper c2 = new() {PreserveCollinear = false};
        c2.AddSubject(BP);
        c2.AddClip(iPoly);

        Paths64 o2 = new();
        c2.Execute(Clipper2Lib.ClipType.Difference, Clipper2Lib.FillRule.EvenOdd, o2);
        
        /* Expected output
        o2 = {List<List<Point64>>} Count = 2
         [0] = {List<Point64>} Count = 4
          [0] = {Point64} 27000,16001,0 
          [1] = {Point64} 16001,16001,0 
          [2] = {Point64} 16001,22999,0 
          [3] = {Point64} 27000,22999,0 
         [1] = {List<Point64>} Count = 4
          [0] = {Point64} 1000,7001,0 
          [1] = {Point64} 1000,13999,0 
          [2] = {Point64} 5999,13999,0 
          [3] = {Point64} 5999,7001,0 
           */

    }
}