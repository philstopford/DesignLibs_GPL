using System.Threading.Tasks;
using Clipper2Lib;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static Paths64 simplify(Paths64 source)
    {
        return pSimplify(source);
    }

    private static Paths64 pSimplify(Paths64 source)
    {
        int sLength = source.Count;
        Paths64 ret = Helper.initedPaths64(sLength);
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
   
    public static Path64 simplify(Path64 iPoints)
    {
        return pSimplify(iPoints);
    }

    private static Path64 pSimplify(Path64 iPoints)
    {
        Clipper64 c = new();
        c.PreserveCollinear = false;
        c.AddSubject(iPoints);
        Paths64 oPoly = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, oPoly);

        oPoly = pReorderXY(oPoly);

        Path64 working = new(oPoly[0]);

        return working;
    }
}