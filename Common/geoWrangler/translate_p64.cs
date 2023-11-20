using System;
using System.Threading.Tasks;
using Clipper2Lib;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static Point64 move(Point64 source, int x, int y)
    {
        return new (source.X + x, source.Y + y, source.Z);
    }

    public static Point64 move(Point64 source, Int64 x, Int64 y)
    {
        return new (source.X + x, source.Y + y, source.Z);
    }

    public static Point64 move(Point64 source, double x, double y)
    {
        return new (source.X + x, source.Y + y, source.Z);
    }

    public static Paths64 move(Paths64 source, decimal x, decimal y)
    {
        return pMove(source, Convert.ToDouble(x), Convert.ToDouble(y));
    }
    public static Paths64 move(Paths64 source, double x, double y)
    {
        return pMove(source, x, y);
    }

    private static Paths64 pMove(Paths64 source, double x, double y)
    {
        int sLength = source.Count;
        Paths64 ret = Helper.initedPaths64(sLength);

#if !GWSINGLETHREADED
        Parallel.For(0, sLength, i =>
#else
            for (int i = 0; i < source.Length; i++)
#endif
            {
                ret[i] = new(pMove(source[i], x, y));
            }
#if !GWSINGLETHREADED
        );
#endif

        return ret;
    }
    
    public static Path64 move(Path64 source, double x, double y)
    {
        return pMove(source, x, y);
    }

    private static Path64 pMove(Path64 source, double x, double y)
    {
        int sLength = source.Count;
        Path64 ret = Helper.initedPath64(sLength);
#if !GWSINGLETHREADED
        Parallel.For(0, sLength, i =>
#else
            for (int i = 0; i < source.Length; i++)
#endif
            {
                ret[i] = new(move(source[i], x, y));
            }
#if !GWSINGLETHREADED
        );
#endif
        return ret;
    }
}