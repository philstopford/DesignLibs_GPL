using System;
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
        PathD ret = new();
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

        return ret;
    }
    
    public static PathsD stripColinear(PathsD source, double angularTolerance = 0.0f)
    {
        return pStripColinear(source, angularTolerance);
    }

    private static PathsD pStripColinear(PathsD source, double angularTolerance = 0.0f)
    {
        int sLength = source.Count;
        PathsD ret = Helper.initedPathsD(sLength);
#if !GWSINGLETHREADED
        Parallel.For(0, sLength, pt =>
#else
            for (int pt = 0; pt < sLength; pt++)
#endif
            {
                ret[pt] = pStripColinear(source[pt], angularTolerance);
            }
#if !GWSINGLETHREADED
        );
#endif

        return ret;
    }

    public static PathD stripColinear(PathD source, double angularTolerance = 0.0f)
    {
        return pStripColinear(source, angularTolerance);
    }

    private static PathD pStripColinear(PathD source, double angularTolerance = 0.0f)
    {
        switch (source.Count)
        {
            case < 3:
                return new(source);
        }

        PathD ret = new(source);
        
        for (int pt = 0; pt < ret.Count; pt++)
        {
            PointD interSection_A, interSection_B, interSection_C;
            switch (pt)
            {
                // Assess angle.
                case 0:
                    interSection_B = new(ret[^1]); // map to last point
                    interSection_C = new(ret[pt]);
                    interSection_A = new(ret[pt + 1]);
                    break;
                default:
                {
                    if (pt == ret.Count - 1) // last point in the list
                    {
                        interSection_B = new(ret[pt - 1]);
                        interSection_C = new(ret[pt]);
                        interSection_A = new(ret[0]); // map to the first point
                    }
                    else
                    {
                        interSection_B = new(ret[pt - 1]);
                        interSection_C = new(ret[pt]);
                        interSection_A = new(ret[pt + 1]);
                    }

                    break;
                }
            }

            double theta = pAngleBetweenPoints(interSection_A, interSection_B, interSection_C, false);

            bool addPoint = true;
            if (pt != 0 && pt != ret.Count - 1)
            {
                if (Math.Abs(theta - 180) <= angularTolerance)
                {
                    addPoint = false;
                }
            }

            switch (addPoint)
            {
                case true:
                    ret.Add(new (ret[pt]));
                    break;
            }
        }
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