using System;
using System.Collections.Generic;
using geoLib;
using System.Linq;
using System.Threading.Tasks;
using Clipper2Lib;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static Point64 move(Point64 source, decimal x, decimal y)
    {
        return new Point64((Int64)(source.X + x), (Int64)(source.Y + y), source.Z);
    }
    public static Point64 move(Point64 source, double x, double y)
    {
        return new Point64((Int64)(source.X + x), (Int64)(source.Y + y), source.Z);
    }
    
    public static Point64 move(Point64 source, int x, int y)
    {
        return new Point64((source.X + x), (source.Y + y), source.Z);
    }

    public static Point64 move(Point64 source, Int64 x, Int64 y)
    {
        return new Point64((source.X + x), (source.Y + y), source.Z);
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
                ret[i] = new(move(source[i], x, y));
            }
#if !GWSINGLETHREADED
        );
#endif

        return ret;
    }
    
    public static Path64 move(Path64 source, decimal x, decimal y)
    {
        return pMove(source, (Int64)(x), (Int64)(y));
    }

    public static Path64 move(Path64 source, double x, double y)
    {
        return pMove(source, (Int64)(x), (Int64)(y));
    }

    public static Path64 move(Path64 source, Int64 x, Int64 y)
    {
        return pMove(source, x, y);
    }

    private static Path64 pMove(Path64 source, Int64 x, Int64 y)
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

    public static PointD move(PointD source, decimal x, decimal y)
    {
        return move(source, (double)x, (double)y);
    }

    public static PointD move(PointD source, int x, int y)
    {
        return move(source, x, y);
    }

    public static PointD move(PointD source, Int64 x, Int64 y)
    {
        return move(source, x, y);
    }

    public static PointD move(PointD source, double x, double y)
    {
        return new PointD(source.x + x, source.y + y, source.z);
    }

    public static PathsD move(PathsD source, decimal x, decimal y)
    {
        return pMove(source, Convert.ToDouble(x), Convert.ToDouble(y));
    }
    public static PathsD move(PathsD source, double x, double y)
    {
        return pMove(source, x, y);
    }

    private static PathsD pMove(PathsD source, double x, double y)
    {
        int sLength = source.Count;
        PathsD ret = Helper.initedPathsD(sLength);

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
    
    public static PathD move(PathD source, decimal x, decimal y)
    {
        return pMove(source, Convert.ToDouble(x), Convert.ToDouble(y));
    }

    public static PathD move(PathD source, double x, double y)
    {
        return pMove(source, x, y);
    }

    private static PathD pMove(PathD source, double x, double y)
    {
        int sLength = source.Count;
        PathD ret = Helper.initedPathD(sLength);
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