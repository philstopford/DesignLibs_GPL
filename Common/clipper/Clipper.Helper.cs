using System;

namespace Clipper2Lib;

public static class Helper
{
    public static Path64 initedPath64(int count)
    {
        Path64 ret = new(count);
        for (int i = 0; i < count; i++)
        {
            ret.Add(new Point64());
        }

        return ret;
    }
    
    public static Paths64 initedPaths64(int count)
    {
        Paths64 ret = new(count);
        for (int i = 0; i < count; i++)
        {
            ret.Add([]);
        }

        return ret;
    }
    
    public static PathD initedPathD(int count)
    {
        PathD ret = new(count);
        for (int i = 0; i < count; i++)
        {
            ret.Add(new PointD());
        }

        return ret;
    }
    
    public static PathsD initedPathsD(int count)
    {
        PathsD ret = new(count);
        for (int i = 0; i < count; i++)
        {
            ret.Add([]);
        }

        return ret;
    }
    
    public static bool PointsEqual(PointD a, PointD b)
    {
        return Math.Abs(a.x - b.x) < 1e-9 && Math.Abs(a.y - b.y) < 1e-9;
    }

    public static PointD Add(PointD a, PointD b) => new PointD(a.x + b.x, a.y + b.y);
    public static PointD Minus(PointD a, PointD b) => new PointD(a.x - b.x, a.y - b.y);
    public static PointD Mult(PointD a, double s) => new PointD(a.x * s, a.y * s);
    public static double Length(PointD p) => Math.Sqrt(p.x * p.x + p.y * p.y);
    public static PointD Mul(PointD a, double s) => new PointD(a.x * s, a.y * s);

    public static PointD Neg(PointD a) => new PointD(-a.x, -a.y);
    public static double Dot(PointD a, PointD b) => a.x * b.x + a.y * b.y;
    public static PointD Mid(PointD a, PointD b) => new PointD((a.x + b.x) / 2.0, (a.y + b.y) / 2.0);

    public static PointD Normalized(PointD p)
    {
        double len = Length(p);
        return len > 0 ? new PointD(p.x / len, p.y / len) : new PointD(0, 0);
    }    
}