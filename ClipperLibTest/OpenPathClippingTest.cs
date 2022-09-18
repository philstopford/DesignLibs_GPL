using System.Collections.Generic;
using Clipper2Lib;

namespace ClipperLibTest;

using Path = Path64;
using Paths = Paths64;

public class OpenPathClippingTest
{
    public static void test()
    {
        Path64 testPath = Clipper.MakePath(new[]
        {
            -50000, -550000,
            -50000, -150000,
            650000, -150000
        });
        
        Path64 b = Clipper.MakePath(new[]
        {
            300000,-800000,
            300000,0,
            500000,0,
            500000,-800000
        });
        
        Clipper64 c = new() {PreserveCollinear = true};
        c.AddOpenSubject(testPath);
        c.AddClip(b);
        Paths64 unused = new();
        Paths64 topChords = new();
        c.Execute(ClipType.Intersection, FillRule.EvenOdd, unused, topChords);

        Path64 testPath2 = Clipper.MakePath(new[]
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
        
        Path64 testPath3 = Clipper.MakePath(new[]
        {
            300000,-800000,
            300000,0
        });

        Path64 a = Clipper.MakePath(new[]
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

        Path64 testPath4 = Clipper.MakePath(new[]
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

        /* Expected output
        bottomChords = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 2
          [0] = {Point64} 500000,-550000,0 
          [1] = {Point64} 300000,-550000,0 
        leftChords = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 2
          [0] = {Point64} 300000,-550000,0 
          [1] = {Point64} 300000,-150000,0 
        rightChords = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 2
          [0] = {Point64} 500000,-150000,0 
          [1] = {Point64} 500000,-550000,0 
        topChords = {List<List<Point64>>} Count = 1
         [0] = {List<Point64>} Count = 2
          [0] = {Point64} 300000,-150000,0 
          [1] = {Point64} 500000,-150000,0 
             */
        
    }
}