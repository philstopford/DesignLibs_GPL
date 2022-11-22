using System;
using Clipper2Lib;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static Paths64 close(Paths64 source)
    {
        return pClose(source);
    }
    
    public static Path64 close(Path64 source)
    {
        return pClose(source);
    }

    private static Path64 pClose(Path64 source)
    {
        switch (source.Count)
        {
            case < 1:
                return source;
        }
        if (Math.Abs(source[0].X - source[^1].X) > double.Epsilon || Math.Abs(source[0].Y - source[^1].Y) > double.Epsilon)
        {
            source.Add(new (source[0]));
        }
        return source;
    }
    
    private static Paths64 pClose(Paths64 source)
    {
        Paths64 ret = new();
        foreach (Path64 p in source)
        {
            ret.Add(pClose(p));
        }

        return ret;
    }
}