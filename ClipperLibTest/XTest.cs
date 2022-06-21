using System.Collections.Generic;

namespace ClipperLibTest;

public static class XTest
{
    public static void compare()
    {
        List<List<Clipper2Lib.Point64>> xShape = new()
        {
            new()
            {
                new(0, 13000),
                new(0, 21000),
                new(13000, 21000),
                new(13000, 30000),
                new(24000, 30000),
                new(24000, 21000),
                new(48000, 21000),
                new(48000, 13000),
                new(24000, 13000),
                new(24000, 0),
                new(13000, 0),
                new(13000, 13000)
            }
        };

        List<Clipper2Lib.Point64> bounds = new()
        {
            new (0,30000),
            new (0,0),
            new (48000,0),
            new (48000,30000),
            new (0,30000)
        };

        Clipper2Lib.Clipper64 c2 = new();
        c2.AddSubject( bounds );
        c2.AddClip( xShape);
        List<List<Clipper2Lib.Point64>> c2out = new();
        c2.Execute(Clipper2Lib.ClipType.Difference, Clipper2Lib.FillRule.EvenOdd, c2out);
        
        /* Expected output
        c2out = {List<List<Point64>>} Count = 4
         [0] = {List<Point64>} Count = 4
          [0] = {Point64} 24000,13000,0 
          [1] = {Point64} 48000,13000,0 
          [2] = {Point64} 48000,0,0 
          [3] = {Point64} 24000,0,0 
         [1] = {List<Point64>} Count = 4
          [0] = {Point64} 24000,21000,0 
          [1] = {Point64} 24000,30000,0 
          [2] = {Point64} 48000,30000,0 
          [3] = {Point64} 48000,21000,0 
         [2] = {List<Point64>} Count = 4
          [0] = {Point64} 0,30000,0 
          [1] = {Point64} 13000,30000,0 
          [2] = {Point64} 13000,21000,0 
          [3] = {Point64} 0,21000,0 
         [3] = {List<Point64>} Count = 4
          [0] = {Point64} 0,0,0 
          [1] = {Point64} 0,13000,0 
          [2] = {Point64} 13000,13000,0 
          [3] = {Point64} 13000,0,0 
           */

        List<List<ClipperLib1.IntPoint>> xShape1 = new();
        for (int p = 0; p < xShape.Count; p++)
        {
            List<ClipperLib1.IntPoint> t = new();
            for (int pt = 0; pt < xShape[p].Count; pt++)
            {
                t.Add(new (xShape[p][pt].X, xShape[p][pt].Y));
            }
            xShape1.Add(t);
        }

        List<ClipperLib1.IntPoint> bounds1 = new();
        for (int pt = 0; pt < bounds.Count; pt++)
        {
            bounds1.Add(new (bounds[pt].X, bounds[pt].Y));
        }

        ClipperLib1.Clipper c1 = new() {PreserveCollinear = false};
        c1.AddPath(bounds1, ClipperLib1.PolyType.ptSubject, true);
        c1.AddPaths(xShape1, ClipperLib1.PolyType.ptClip, true);

        List<List<ClipperLib1.IntPoint>> o1 = new();
        c1.Execute(ClipperLib1.ClipType.ctDifference, o1);

        /* Expected output
         */
        
    }
}