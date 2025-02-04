using System;
using System.Linq;
using System.Threading.Tasks;
using Clipper2Lib;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static PathD removeDuplicates(PathD source, double threshold = Constants.tolerance)
    {
        return pRemoveDuplicates(source, threshold);
    }
    
    private static PathsD pRemoveDuplicates(PathsD source, double threshold = Constants.tolerance)
    {
        PathsD ret = [];
        ret.AddRange(source.Select(t => removeDuplicates(t, threshold)));

        return ret;
    }
    
    private static PathD pRemoveDuplicates(PathD source, double threshold = Constants.tolerance)
    {
        // Experiment to see if we can use Clipper instead.
        PathD ret = Clipper.StripNearDuplicates(source, threshold, distanceBetweenPoints(source[0], source[^1]) <= Constants.tolerance);
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
                return new PathD(source);
        }

        PathD ret = Clipper.TrimCollinear(source, precision, distanceBetweenPoints(source[0], source[^1]) < Constants.tolerance);

        return ret;
    }

    public static PathD stripTerminators(PathD source, bool keepLast)
    {
        return pStripTerminators(source, keepLast);
    }
    
    private static PathD pStripTerminators(PathD source_, bool keepLast)
    {
        PathD ret = new(source_);
        try
        {
            bool firstLast_same = false;
        int pt_Check = ret.Count - 1;
        double dist = distanceBetweenPoints(ret[pt_Check], ret[0]);
        if (dist < Constants.tolerance)
        {
            firstLast_same = true; // remove duplicated points. The shape will be closed later.
        }
        while (firstLast_same)
        {
            ret.RemoveAt(pt_Check); // remove duplicated points. The shape will be closed later
            pt_Check--;

            // We ought not to get here - this is something for debug checks.
            /*
            if (pt_Check >= ret.Count || pt_Check < 0)
            {
                int x = 2;
            }
            */
            dist = distanceBetweenPoints(ret[pt_Check], ret[0]);
            if (dist > Constants.tolerance)
            {
                firstLast_same = false; // stop at the first unmatched point.
            }
        }
        }
        catch (Exception e)
        {
            Console.WriteLine(e);
            throw;
        }

        ret = keepLast switch
        {
            true => pClose(ret),
            _ => ret
        };

        return ret;
    }
}