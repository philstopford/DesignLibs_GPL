using System;
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
                return source;
        }
        if (Math.Abs(source[0].x - source[^1].x) > double.Epsilon || Math.Abs(source[0].y - source[^1].y) > double.Epsilon)
        {
            source.Add(new (source[0]));
        }
        return source;
    }
    
    private static PathsD pClose(PathsD source)
    {
        PathsD ret = new();
        foreach (PathD p in source)
        {
            ret.Add(pClose(p));
        }

        return ret;
    }
}