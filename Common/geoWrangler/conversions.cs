using Clipper2Lib;
using geoLib;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static Paths64 pathsFromPoints(Paths64 source, long scaling)
    {
        return pPathsFromPoints(source, scaling);
    }

    private static Paths64 pPathsFromPoints(Paths64 source, long scaling)
    {
        Paths64 ret = new();
        try
        {
            ret.AddRange(source.Select(t => pathFromPoint(t, scaling)));
        }
        catch
        {
            // ignored
        }

        return ret;
    }

    public static Path64 pathFromPoint(Path64 source, long scaling)
    {
        return pPathFromPoint(source, scaling);
    }

    private static Path64 pPathFromPoint(Path64 source, long scaling)
    {
        int length = source.Count;
        if (source[0].X != source[^1].X && source[0].Y != source[^1].Y)
        {
            length++; // close the geometry
        }
        Path64 returnPath = new();
        try
        {
            returnPath.AddRange(source.Select(t => new Point64(t.X * scaling, t.Y * scaling)));
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

    public static Paths64 pointsFromPaths(Paths64 source, long scaling)
    {
        return pPointsFromPaths(source, scaling);
    }

    private static Paths64 pPointsFromPaths(Paths64 source, long scaling)
    {
        return new (source.Select(t => pPointFromPath(t, scaling)));
    }

    public static Path64 pointFromPath(Path64 source, long scaling)
    {
        return pPointFromPath(source, scaling);
    }

    private static Path64 pPointFromPath(Path64 source, long scaling)
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
        Path64 returnPoint = new (length);
#if !GWSINGLETHREADED
        Parallel.For(0, sCount, pt =>
#else
            for (int pt = 0; pt < source.Count; pt++)
#endif
            {
                
                returnPoint[pt] = new ((long)Math.Round(Convert.ToDecimal(source[pt].X) / scaling), (long)Math.Round(Convert.ToDecimal(source[pt].Y) / scaling));
            }
#if !GWSINGLETHREADED
        );
#endif
        // Close the shape.
        if (length != sCount)
        {
            returnPoint[length - 1] = new (returnPoint[0]);
        }
        return returnPoint;
    }

    public static Paths64 pathsFromPointFs(PathsD source, long scaling)
    {
        return pPathsFromPointFs(source, scaling);
    }

    private static Paths64 pPathsFromPointFs(PathsD source, long scaling)
    {
        Paths64 ret = new();
        try
        {
            ret.AddRange(source.Select(t => pathFromPointF(t, scaling)));
        }
        catch
        {
            // ignored
        }

        return ret;
    }

    public static Path64 pathFromPointF(PathD source, long scaling)
    {
        return pPathFromPointF(source, scaling);
    }

    private static Path64 pPathFromPointF(PathD source, long scaling)
    {
        int length = source.Count;
        switch (Math.Abs(source[0].x - source[^1].x))
        {
            case > double.Epsilon when Math.Abs(source[0].y - source[^1].y) > double.Epsilon:
                length++; // close the geometry
                break;
        }
        Path64 returnPath = new();
        try
        {
            returnPath.AddRange(source.Select(t => new Point64(Convert.ToInt64(t.x * scaling), Convert.ToInt64(t.y * scaling))));
        }
        catch
        {
            // ignored
        }

        // Close the shape
        if (length != source.Count)
        {
            returnPath.Add(new (returnPath[0]));
        }
        return returnPath;
    }

    public static PathsD pointFsFromPaths(Paths64 source, long scaling)
    {
        return pPointFsFromPaths(source, scaling);
    }

    private static PathsD pPointFsFromPaths(Paths64 source, long scaling)
    {
        return new(source.Select(t => pPointFFromPath(t, scaling)));
    }

    public static PathD pointFFromPath(Path64 source, long scaling)
    {
        return pPointFFromPath(source, scaling);
    }

    private static PathD pPointFFromPath(Path64 source, long scaling)
    {
        int length = source.Count;
        int sourceCount = length;
        if (source[0].X != source[length - 1].X && source[0].Y != source[length - 1].Y)
        {
            length++; // close the geometry
        }
        PathD returnPointF = new (length);
#if !GWSINGLETHREADED
        Parallel.For(0, sourceCount, pt =>
#else
            for (int pt = 0; pt < source.Count(); pt++)
#endif
            {
                returnPointF[pt] = new ((double)source[pt].X / scaling,
                    (double)source[pt].Y / scaling);
            }
#if !GWSINGLETHREADED
        );
#endif
        // Close the shape.
        if (length != sourceCount)
        {
            returnPointF[length - 1] = new (returnPointF[0]);
        }
        return returnPointF;
    }

    public static Paths64 pointsFromPointFs(PathsD source, long scaling)
    {
        return pPointsFromPointFs(source, scaling);
    }

    private static Paths64 pPointsFromPointFs(PathsD source, long scaling)
    {
        Paths64 ret = new();
        try
        {
            ret.AddRange(source.Select(t => pPointsFromPointF(t, scaling)));
        }
        catch
        {
            // ignored
        }

        return ret;
    }

    public static Path64 pointsFromPointF(PathD source, long scaling)
    {
        return pPointsFromPointF(source, scaling);
    }

    private static Path64 pPointsFromPointF(PathD source, long scaling)
    {
        Path64 ret = new (source.Count);
        for (int pt = 0; pt < source.Count; pt++)
        {
            try
            {
                ret[pt] = new (Convert.ToInt64(source[pt].x * scaling),
                    Convert.ToInt64(source[pt].y * scaling));
            }
            catch
            {
                // ignored
            }
        }

        return ret;
    }

    public static PathsD pointFsFromPoints(Paths64 source, long scaling)
    {
        return pPointFsFromPoints(source, scaling);
    }

    private static PathsD pPointFsFromPoints(Paths64 source, long scaling)
    {
        PathsD ret = new();
        try
        {
            ret.AddRange(source.Select(t => pPointFsFromPoint(t, scaling)));
        }
        catch
        {
            // ignored
        }

        return ret;
    }

    public static PathD pointFsFromPoint(Path64 source, long scaling)
    {
        return pPointFsFromPoint(source, scaling);
    }

    private static PathD pPointFsFromPoint(Path64 source, long scaling)
    {
        PathD ret = new (source.Count);
        for (int pt = 0; pt < source.Count; pt++)
        {
            try
            {
                ret[pt] = new (Convert.ToDouble(source[pt].X) / scaling,
                    Convert.ToDouble(source[pt].Y) / scaling);
            }
            catch
            {
                // ignored
            }
        }

        return ret;
    }

}