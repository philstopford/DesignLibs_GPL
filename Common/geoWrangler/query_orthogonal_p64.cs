using System;
using System.Linq;
using Clipper2Lib;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static bool pOrthogonal(Paths64 sourcePoly, double angularTolerance)
    {
        bool ret = true;
        foreach (Path64 p in sourcePoly)
        {
            ret = ret && pOrthogonal(p, angularTolerance);
        }

        return ret;
    }

    public static bool orthogonal(Path64 sourcePoly, double angularTolerance)
    {
        return pOrthogonal(sourcePoly, angularTolerance);
    }

    private static bool pOrthogonal(Path64 sourcePoly, double angularTolerance)
    {
        double[] _angles = angles(sourcePoly, allowNegative: false);

        return _angles.All(t => !(Math.Abs(Math.Abs(t) - 90.0) > angularTolerance));
    }
}