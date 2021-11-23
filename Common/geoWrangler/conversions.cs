using ClipperLib;
using geoLib;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace geoWrangler;

using Path = List<IntPoint>;
using Paths = List<List<IntPoint>>;

public static partial class GeoWrangler
{
    public static Paths pathsFromPoints(List<GeoLibPoint[]> source, long scaling)
    {
        return pPathsFromPoints(source, scaling);
    }

    private static Paths pPathsFromPoints(List<GeoLibPoint[]> source, long scaling)
    {
        Paths ret = new();
        try
        {
            foreach (GeoLibPoint[] t in source)
            {
                ret.Add(pathFromPoint(t, scaling));
            }
        }
        catch (Exception)
        {
        }
        return ret;
    }

    public static Path pathFromPoint(GeoLibPoint[] source, long scaling)
    {
        return pPathFromPoint(source, scaling);
    }

    private static Path pPathFromPoint(GeoLibPoint[] source, long scaling)
    {
        int length = source.Length;
        if (source[0].X != source[^1].X && source[0].Y != source[^1].Y)
        {
            length++; // close the geometry
        }
        Path returnPath = new();
        try
        {
            foreach (GeoLibPoint t in source)
            {
                returnPath.Add(new IntPoint(t.X * scaling, t.Y * scaling));
            }
        }
        catch (Exception)
        {
        }

        // Close the shape
        if (length != source.Length)
        {
            returnPath.Add(new IntPoint(returnPath[0]));
        }
        return returnPath;
    }

    public static List<GeoLibPoint[]> pointsFromPaths(Paths source, long scaling)
    {
        return pPointsFromPaths(source, scaling);
    }

    private static List<GeoLibPoint[]> pPointsFromPaths(Paths source, long scaling)
    {
        List<GeoLibPoint[]> ret = new();
        foreach (Path t in source)
        {
            ret.Add(pPointFromPath(t, scaling));
        }

        return ret;
    }

    public static GeoLibPoint[] pointFromPath(Path source, long scaling)
    {
        return pPointFromPath(source, scaling);
    }

    private static GeoLibPoint[] pPointFromPath(Path source, long scaling)
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
        GeoLibPoint[] returnPoint = new GeoLibPoint[length];
#if !GWSINGLETHREADED
        Parallel.For(0, sCount, pt =>
#else
            for (int pt = 0; pt < source.Count; pt++)
#endif
            {
                
                returnPoint[pt] = new GeoLibPoint((long)Math.Round(Convert.ToDecimal(source[pt].X) / scaling), (long)Math.Round(Convert.ToDecimal(source[pt].Y) / scaling));
            }
#if !GWSINGLETHREADED
        );
#endif
        // Close the shape.
        if (length != sCount)
        {
            returnPoint[length - 1] = new GeoLibPoint(returnPoint[0]);
        }
        return returnPoint;
    }

    public static Paths pathsFromPointFs(List<GeoLibPointF[]> source, long scaling)
    {
        return pPathsFromPointFs(source, scaling);
    }

    private static Paths pPathsFromPointFs(List<GeoLibPointF[]> source, long scaling)
    {
        Paths ret = new();
        try
        {
            foreach (GeoLibPointF[] t in source)
            {
                ret.Add(pathFromPointF(t, scaling));
            }
        }
        catch (Exception)
        {
        }
        return ret;
    }

    public static Path pathFromPointF(GeoLibPointF[] source, long scaling)
    {
        return pPathFromPointF(source, scaling);
    }

    private static Path pPathFromPointF(GeoLibPointF[] source, long scaling)
    {
        int length = source.Length;
        switch (Math.Abs(source[0].X - source[^1].X))
        {
            case > double.Epsilon when Math.Abs(source[0].Y - source[^1].Y) > double.Epsilon:
                length++; // close the geometry
                break;
        }
        Path returnPath = new();
        try
        {
            foreach (GeoLibPointF t in source)
            {
                returnPath.Add(new IntPoint(Convert.ToInt64(t.X * scaling),
                    Convert.ToInt64(t.Y * scaling)));
            }
        }
        catch (Exception)
        {
        }

        // Close the shape
        if (length != source.Length)
        {
            returnPath.Add(new IntPoint(returnPath[0]));
        }
        return returnPath;
    }

    public static List<GeoLibPointF[]> pointFsFromPaths(Paths source, long scaling)
    {
        return pPointFsFromPaths(source, scaling);
    }

    private static List<GeoLibPointF[]> pPointFsFromPaths(Paths source, long scaling)
    {
        List<GeoLibPointF[]> ret = new();
        foreach (Path t in source)
        {
            ret.Add(pPointFFromPath(t, scaling));
        }

        return ret;
    }

    public static GeoLibPointF[] pointFFromPath(Path source, long scaling)
    {
        return pPointFFromPath(source, scaling);
    }

    private static GeoLibPointF[] pPointFFromPath(Path source, long scaling)
    {
        int length = source.Count;
        int sourceCount = length;
        if (source[0].X != source[length - 1].X && source[0].Y != source[length - 1].Y)
        {
            length++; // close the geometry
        }
        GeoLibPointF[] returnPointF = new GeoLibPointF[length];
#if !GWSINGLETHREADED
        Parallel.For(0, sourceCount, pt =>
#else
            for (int pt = 0; pt < source.Count(); pt++)
#endif
            {
                returnPointF[pt] = new GeoLibPointF((double)source[pt].X / scaling,
                    (double)source[pt].Y / scaling);
            }
#if !GWSINGLETHREADED
        );
#endif
        // Close the shape.
        if (length != sourceCount)
        {
            returnPointF[length - 1] = new GeoLibPointF(returnPointF[0]);
        }
        return returnPointF;
    }

    public static List<GeoLibPoint[]> pointsFromPointFs(List<GeoLibPointF[]> source, long scaling)
    {
        return pPointsFromPointFs(source, scaling);
    }

    private static List<GeoLibPoint[]> pPointsFromPointFs(List<GeoLibPointF[]> source, long scaling)
    {
        List<GeoLibPoint[]> ret = new();
        try
        {
            foreach (GeoLibPointF[] t in source)
            {
                ret.Add(pPointsFromPointF(t, scaling));
            }
        }
        catch (Exception)
        {
        }
        return ret;
    }

    public static GeoLibPoint[] pointsFromPointF(GeoLibPointF[] source, long scaling)
    {
        return pPointsFromPointF(source, scaling);
    }

    private static GeoLibPoint[] pPointsFromPointF(GeoLibPointF[] source, long scaling)
    {
        GeoLibPoint[] ret = new GeoLibPoint[source.Length];
        for (int pt = 0; pt < source.Length; pt++)
        {
            try
            {
                ret[pt] = new GeoLibPoint(Convert.ToInt64(source[pt].X * scaling),
                    Convert.ToInt64(source[pt].Y * scaling));
            }
            catch (Exception)
            {
            }
        }

        return ret;
    }

    public static List<GeoLibPointF[]> pointFsFromPoints(List<GeoLibPoint[]> source, long scaling)
    {
        return pPointFsFromPoints(source, scaling);
    }

    private static List<GeoLibPointF[]> pPointFsFromPoints(List<GeoLibPoint[]> source, long scaling)
    {
        List<GeoLibPointF[]> ret = new();
        try
        {
            foreach (GeoLibPoint[] t in source)
            {
                ret.Add(pPointFsFromPoint(t, scaling));
            }
        }
        catch (Exception)
        {
        }
        return ret;
    }

    public static GeoLibPointF[] pointFsFromPoint(GeoLibPoint[] source, long scaling)
    {
        return pPointFsFromPoint(source, scaling);
    }

    private static GeoLibPointF[] pPointFsFromPoint(GeoLibPoint[] source, long scaling)
    {
        GeoLibPointF[] ret = new GeoLibPointF[source.Length];
        for (int pt = 0; pt < source.Length; pt++)
        {
            try
            {
                ret[pt] = new GeoLibPointF(Convert.ToDouble(source[pt].X) / scaling,
                    Convert.ToDouble(source[pt].Y) / scaling);
            }
            catch (Exception)
            {
            }
        }

        return ret;
    }

}