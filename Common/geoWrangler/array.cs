using System;
using Clipper2Lib;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static PathsD makeArray2(PathD source, int xCount, double xPitch, int yCount, double yPitch)
    {
        return pMakeArray(source, xCount, xPitch, yCount, yPitch);
    }
    
    public static PathsD makeArray(PathD source, int xCount, decimal xPitch, int yCount, decimal yPitch)
    {
        return pMakeArray(source, xCount, Convert.ToDouble(xPitch), yCount, Convert.ToDouble(yPitch));
    }

    public static PathsD makeArray(PathD source, int xCount, double xPitch, int yCount, double yPitch)
    {
        return pMakeArray(source, xCount, xPitch, yCount, yPitch);
    }

    private static PathsD pMakeArray(PathD source, int xCount, double xPitch, int yCount, double yPitch)
    {
        PathsD ret = new();
        for (int x = 0; x < xCount; x++)
        {
            for (int y = 0; y < yCount; y++)
            {
                ret.Add(pMove(source, x * xPitch, y * yPitch));
            }
        }

        return ret;
    }

    public static PathsD makeArray(PathsD source, int xCount, decimal xPitch, int yCount, decimal yPitch)
    {
        return pMakeArray(source, xCount, Convert.ToDouble(xPitch), yCount, Convert.ToDouble(yPitch));
    }

    public static PathsD makeArray(PathsD source, int xCount, double xPitch, int yCount, double yPitch)
    {
        return pMakeArray(source, xCount, xPitch, yCount, yPitch);
    }

    private static PathsD pMakeArray(PathsD source, int xCount, double xPitch, int yCount, double yPitch)
    {
        PathsD ret = new();
        for (int x = 0; x < xCount; x++)
        {
            for (int y = 0; y < yCount; y++)
            {
                ret.AddRange(pMove(source, x * xPitch, y * yPitch));
            }
        }

        return ret;
    }
}