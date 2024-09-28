using System;
using System.Linq;
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

        Path64 ret = new(source);
        if (Math.Abs(source[0].X - source[^1].X) > Constants.tolerance || Math.Abs(source[0].Y - source[^1].Y) > Constants.tolerance)
        {
            ret.Add(new Point64(source[0]));
        }
        return ret;
    }
    
    private static Paths64 pClose(Paths64 source)
    {
        Paths64 ret = [];
        ret.AddRange(source.Select(pClose));

        return ret;
    }
}