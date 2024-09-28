using System;
using System.Threading.Tasks;
using Clipper2Lib;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static PointD move(PointD source, decimal x, decimal y)
    {
        return new PointD(source.x + (double)x, source.y + (double)y);
    }

    public static PointD move(PointD source, int x, int y)
    {
        return new PointD(source.x + x, source.y + y);
    }

    public static PointD move(PointD source, long x, long y)
    {
        return new PointD(source.x + x, source.y + y);
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
                ret[i] = new PathD(pMove(source[i], x, y));
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
                ret[i] = new PointD(move(source[i], x, y));
            }
#if !GWSINGLETHREADED
        );
#endif
        return ret;
    }
}