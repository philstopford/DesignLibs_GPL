using Clipper2Lib;
using geoLib;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static Paths64 paths64FromPathsD(PathsD source, double scaling = 1.0)
    {
        return pPaths64FromPathsD(source, scaling);
    }

    private static Paths64 pPaths64FromPathsD(PathsD source, double scaling = 1.0)
    {
        Paths64 ret = new();
        try
        {
            ret.AddRange(source.Select(t => path64FromPathD(t, scaling)));
        }
        catch
        {
            // ignored
        }

        return ret;
    }

    public static Path64 path64FromPathD(PathD source, double scaling = 1.0)
    {
        return pPath64FromPathD(source, scaling);
    }

    private static Path64 pPath64FromPathD(PathD source, double scaling = 1.0)
    {
        int length = source.Count;
        if (source[0].x != source[^1].x && source[0].y != source[^1].y)
        {
            length++; // close the geometry
        }
        Path64 returnPath = new();
        try
        {
            returnPath.AddRange(source.Select(t => new Point64(t)));
        }
        catch
        {
            // ignored
        }

        if (scaling != 1.0)
        {
            returnPath = Clipper.ScalePath(returnPath, scaling);
        }

        // Close the shape
        if (length != source.Count)
        {
            returnPath.Add(new Point64(returnPath[0]));
        }
        return returnPath;
    }

    public static PathsD pathsDFromPaths64(Paths64 source, double scaling = 1.0)
    {
        return pPathsDFromPaths64(source, scaling);
    }

    private static PathsD pPathsDFromPaths64(Paths64 source, double scaling = 1.0)
    {
        return new (source.Select(t => pPathDFromPath64(t, scaling)));
    }

    public static PathD pathDFromPath64(Path64 source, double scaling = 1.0)
    {
        return pPathDFromPath64(source, scaling);
    }

    private static PathD pPathDFromPath64(Path64 source, double scaling = 1.0)
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
}