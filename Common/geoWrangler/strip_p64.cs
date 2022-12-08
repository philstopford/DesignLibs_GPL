using System;
using System.Threading.Tasks;
using Clipper2Lib;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static Path64 removeDuplicates(Path64 source, double threshold = constants.tolerance)
    {
        return pRemoveDuplicates(source, threshold);
    }
    
    private static Paths64 pRemoveDuplicates(Paths64 source, double threshold = constants.tolerance)
    {
        Paths64 ret = new ();
        foreach (var t in source)
        {
            ret.Add(removeDuplicates(t, threshold));
        }

        return ret;
    }
    
    private static Path64 pRemoveDuplicates(Path64 source, double threshold = constants.tolerance)
    {
        Path64 ret = new();
        switch (source.Count)
        {
            case > 0:
            {
                ret.Add(new (source[0]));
                int retIndex = 1;
                for (int i = 1; i < source.Count - 1; i++)
                {
                    if (!(Math.Abs(source[i].X - ret[retIndex - 1].X) > threshold) &&
                        !(Math.Abs(source[i].Y - ret[retIndex - 1].Y) > threshold))
                    {
                        continue;
                    }

                    ret.Add(new (source[i]));
                    retIndex++;
                }

                break;
            }
        }

        return ret;
    }
    
    public static Paths64 stripCollinear(Paths64 source)
    {
        return pStripCollinear(source);
    }

    private static Paths64 pStripCollinear(Paths64 source)
    {
        int sLength = source.Count;
        Paths64 ret = Helper.initedPaths64(sLength);
#if !GWSINGLETHREADED
        Parallel.For(0, sLength, pt =>
#else
            for (int pt = 0; pt < sLength; pt++)
#endif
            {
                ret[pt] = pStripCollinear(source[pt]);
            }
#if !GWSINGLETHREADED
        );
#endif

        return ret;
    }

    public static Path64 stripCollinear(Path64 source)
    {
        return pStripCollinear(source);
    }

    private static Path64 pStripCollinear(Path64 source)
    {
        switch (source.Count)
        {
            case < 3:
                return new(source);
        }

        Path64 ret = Clipper.TrimCollinear(source, distanceBetweenPoints(source[0], source[^1]) == 0);

        return ret;
    }
    
    public static Path64 stripTerminators(Path64 source, bool keepLast)
    {
        return pStripTerminators(source, keepLast);
    }
    
    private static Path64 pStripTerminators(Path64 source_, bool keepLast)
    {
        Path64 ret = new(source_);
        bool firstLast_same = false;
        int pt_Check = ret.Count - 1;
        if (distanceBetweenPoints(ret[pt_Check], ret[0]) < 0.01)
        {
            firstLast_same = true; // remove duplicated points. The shape will be closed later.
        }
        while (firstLast_same)
        {
            ret.RemoveAt(pt_Check); // remove duplicated points. The shape will be closed later
            pt_Check--;
            if (distanceBetweenPoints(ret[pt_Check], ret[0]) > 0.01)
            {
                firstLast_same = false; // stop at the first unmatched point.
            }
        }

        ret = keepLast switch
        {
            true => pClose(ret),
            _ => ret
        };

        return ret;
    }
}