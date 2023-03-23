using System.Threading.Tasks;
using Clipper2Lib;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static PathD removeDuplicates(PathD source, double threshold = constants.tolerance)
    {
        return pRemoveDuplicates(source, threshold);
    }
    
    private static PathsD pRemoveDuplicates(PathsD source, double threshold = constants.tolerance)
    {
        PathsD ret = new ();
        foreach (var t in source)
        {
            ret.Add(removeDuplicates(t, threshold));
        }

        return ret;
    }
    
    private static PathD pRemoveDuplicates(PathD source, double threshold = constants.tolerance)
    {
        // Experiment to see if we can use Clipper instead.
        PathD ret = Clipper.StripNearDuplicates(source, threshold, distanceBetweenPoints(source[0], source[^1]) <= constants.tolerance);
        /*
        switch (source.Count)
        {
            case > 0:
            {
                ret.Add(new (source[0]));
                int retIndex = 1;
                for (int i = 1; i < source.Count - 1; i++)
                {
                    if (!(Math.Abs(source[i].x - ret[retIndex - 1].x) > threshold) &&
                        !(Math.Abs(source[i].y - ret[retIndex - 1].y) > threshold))
                    {
                        continue;
                    }

                    ret.Add(new (source[i]));
                    retIndex++;
                }

                break;
            }
        }
        */
        return pClose(ret);
    }
    
    public static PathsD stripCollinear(PathsD source)
    {
        return pStripCollinear(source);
    }

    private static PathsD pStripCollinear(PathsD source)
    {
        int sLength = source.Count;
        PathsD ret = Helper.initedPathsD(sLength);
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

    public static PathD stripCollinear(PathD source, int precision = 2)
    {
        return pStripCollinear(source, precision);
    }

    private static PathD pStripCollinear(PathD source, int precision = 2)
    {
        switch (source.Count)
        {
            case < 3:
                return new(source);
        }

        PathD ret = Clipper.TrimCollinear(source, precision, distanceBetweenPoints(source[0], source[^1]) < constants.tolerance);

        return ret;
    }

    public static PathD stripTerminators(PathD source, bool keepLast)
    {
        return pStripTerminators(source, keepLast);
    }
    
    private static PathD pStripTerminators(PathD source_, bool keepLast)
    {
        PathD ret = new(source_);
        bool firstLast_same = false;
        int pt_Check = ret.Count - 1;
        if (distanceBetweenPoints(ret[pt_Check], ret[0]) < constants.tolerance)
        {
            firstLast_same = true; // remove duplicated points. The shape will be closed later.
        }
        while (firstLast_same)
        {
            ret.RemoveAt(pt_Check); // remove duplicated points. The shape will be closed later
            pt_Check--;
            if (distanceBetweenPoints(ret[pt_Check], ret[0]) > constants.tolerance)
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