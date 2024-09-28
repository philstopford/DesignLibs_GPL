using System;
using System.Linq;
using Clipper2Lib;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static PathsD close(PathsD source)
    {
        return pClose(source);
    }
    
    public static PathD close(PathD source)
    {
        return pClose(source);
    }

    private static PathD pClose(PathD source)
    {
        switch (source.Count)
        {
            case < 1:
                return new PathD(source);
        }

        PathD ret = new(source);
        if (Math.Abs(source[0].x - source[^1].x) > Constants.tolerance || Math.Abs(source[0].y - source[^1].y) > Constants.tolerance)
        {
            ret.Add(new PointD(source[0]));
        }
        return ret;
    }
    
    private static PathsD pClose(PathsD source)
    {
        PathsD ret = [];
        ret.AddRange(source.Select(pClose));

        return ret;
    }
}