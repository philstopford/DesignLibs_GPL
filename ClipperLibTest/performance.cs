using System.Diagnostics;

namespace ClipperLibTest;

using Path64 = Clipper2Lib.Path64;
using Paths64 = Clipper2Lib.Paths64;
using Path = List<ClipperLib1.IntPoint>;
using Paths = List<List<ClipperLib1.IntPoint>>;

public class performance
{
    private static Clipper2Lib.Point64 MakeRandomPt(int maxWidth, int maxHeight, Random rand)
    {
        Int64 x  = rand.Next(maxWidth);
        Int64 y = rand.Next(maxHeight);
        return new Clipper2Lib.Point64(x,y);
    }
        
    public static Path64 MakeRandomPath(int width, int height, int count, Random rand)
    {
        rand.Next();
        Path64 result = new Path64(count);
        for (int i = 0; i < count; ++i)
            result.Add(MakeRandomPt(width, height, rand));
        return result;
    }

    public static void compare()
    {
        const int displayWidth = 800, displayHeight = 600, edgeCnt = 1000;
        Random rand = new ();
        Stopwatch stopwatch = new ();
        Paths64 subj2 = new ();
        Paths64 clip2 = new ();
        Paths64 solution2 = new ();

        subj2.Add(MakeRandomPath(displayWidth, displayHeight, edgeCnt, rand));
        clip2.Add(MakeRandomPath(displayWidth, displayHeight, edgeCnt, rand));

        Clipper2Lib.Clipper64 c2 = new ();
        c2.AddSubject(subj2);
        c2.AddClip(clip2);
        stopwatch.Start();
        c2.Execute(Clipper2Lib.ClipType.Intersection, Clipper2Lib.FillRule.NonZero, solution2);
        stopwatch.Stop();
        
        double c2_elapsed = stopwatch.Elapsed.TotalSeconds;
        
        stopwatch.Reset();

        Paths subj1 = new();
        for (int p = 0; p < subj2.Count; p++)
        {
            Path tmp = new();
            for (int pt = 0; pt < subj2[p].Count; pt++)
            {
                tmp.Add(new (subj2[p][pt].X, subj2[p][pt].Y));
            }
            subj1.Add(tmp);
        }
        
        Paths clip1 = new();
        for (int p = 0; p < clip2.Count; p++)
        {
            Path tmp = new();
            for (int pt = 0; pt < clip2[p].Count; pt++)
            {
                tmp.Add(new (clip2[p][pt].X, clip2[p][pt].Y));
            }
            clip1.Add(tmp);
        }

        Paths solution1 = new();

        ClipperLib1.Clipper c1 = new ();
        c1.AddPaths(subj1, ClipperLib1.PolyType.ptSubject, true);
        c1.AddPaths(clip1, ClipperLib1.PolyType.ptClip, true);
        stopwatch.Start();
        c1.Execute(ClipperLib1.ClipType.ctIntersection, solution1, ClipperLib1.PolyFillType.pftNonZero);
        stopwatch.Stop();

        double c1_elapsed = stopwatch.Elapsed.TotalSeconds;

        Console.WriteLine("Clipper1: " + c1_elapsed);
        Console.WriteLine("Clipper2: " + c2_elapsed);
    }
}