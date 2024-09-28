using System.Threading.Tasks;
using Clipper2Lib;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static PathsD simplify(PathsD source)
    {
        return pSimplify(source);
    }

    private static PathsD pSimplify(PathsD source)
    {
        int sLength = source.Count;
        PathsD ret = Helper.initedPathsD(sLength);
#if !GWSINGLETHREADED
        Parallel.For(0, sLength, pt =>
#else
            for (int pt = 0; pt < sLength; pt++)
#endif
            {
                ret[pt] = pSimplify(source[pt]);
            }
#if !GWSINGLETHREADED
        );
#endif
        return ret;
    }
   
    public static PathD simplify(PathD iPoints)
    {
        return pSimplify(iPoints);
    }

    private static PathD pSimplify(PathD iPoints)
    {
        ClipperD c = new(Constants.roundingDecimalPrecision)
        {
            PreserveCollinear = false
        };
        c.AddSubject(iPoints);
        PathsD oPoly = [];
        c.Execute(ClipType.Union, FillRule.EvenOdd, oPoly);

        oPoly = pReorderXY(oPoly);

        PathD working = new(oPoly[0]);

        return working;
    }
}