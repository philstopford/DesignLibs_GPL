using System;
using Clipper2Lib;
using System.Linq;
using System.Threading.Tasks;

namespace geoWrangler;

public static partial class GeoWrangler
{

    private static Paths64 _pPaths64FromPathsD(PathsD source, double scaling = 1.0)
    {
        Paths64 ret = new();
        try
        {
            ret.AddRange(source.Select(t => _pPath64FromPathD(t, scaling)));
        }
        catch
        {
            // ignored
        }

        return ret;
    }

    private static Path64 _pPath64FromPathD(PathD source, double scaling = 1.0)
    {
        int length = source.Count;
        if (Math.Abs(source[0].x - source[^1].x) > constants.tolerance && Math.Abs(source[0].y - source[^1].y) > constants.tolerance)
        {
            length++; // close the geometry
        }
        Path64 returnPath = new();
        try
        {
            returnPath.AddRange(source.Select(t => new Point64(t.x * scaling, t.y * scaling)));
        }
        catch
        {
            // ignored
        }
        
        // Close the shape
        if (length != source.Count)
        {
            returnPath.Add(new Point64(returnPath[0]));
        }
        return returnPath;
    }

    private static PathsD _pPathsDFromPaths64(Paths64 source, double scaling = 1.0)
    {
        return new (source.Select(t => _pPathDFromPath64(t, scaling)));
    }

    public static PathD PathDFromPath64(Path64 source, double scaling = 1.0)
    {
        return _pPathDFromPath64(source, scaling);
    }
    private static PathD _pPathDFromPath64(Path64 source, double scaling = 1.0)
    {
        int length = source.Count;
        int sCount = length;
        switch (length)
        {
            case > 1:
            {
                if (source[0].X != source[sCount - 1].X && source[0].Y != source[sCount - 1].Y)
                {
                    length++; // close the geometry
                }

                break;
            }
        }
        PathD returnPath = Helper.initedPathD(length);
#if !GWSINGLETHREADED
        Parallel.For(0, sCount, pt =>
#else
            for (int pt = 0; pt < source.Count; pt++)
#endif
            {
                
                returnPath[pt] = new (source[pt]);
            }
#if !GWSINGLETHREADED
        );
#endif

        if (scaling != 1.0)
        {
            returnPath = Clipper.ScalePath(returnPath, scaling);
        }
        
        // Close the shape.
        if (length != sCount)
        {
            returnPath[length - 1] = new (returnPath[0]);
        }
        return returnPath;
    }

    public static Paths64 paths64FromPathsD(PathsD source, double scaling = 1.0)
    {
        return _pPaths64FromPathsD(source, scaling);
    }

    public static Path64 path64FromPathD(PathD source, double scaling = 1.0)
    {
        return _pPath64FromPathD(source, scaling);
    }


    public static PathsD pathsDFromPaths64(Paths64 source, double scaling = 1.0)
    {
        return _pPathsDFromPaths64(source, scaling);
    }

    public static PathD pathDFromPath64(Path64 source, double scaling = 1.0)
    {
        return _pPathDFromPath64(source, scaling);
    }

}