using System;
using System.Threading.Tasks;
using Clipper2Lib;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static Path64 removeDuplicates(Path64 source, double threshold = Double.Epsilon)
    {
        return pRemoveDuplicates(source, threshold);
    }
    
    private static Paths64 pRemoveDuplicates(Paths64 source, double threshold = Double.Epsilon)
    {
        Paths64 ret = new ();
        foreach (var t in source)
        {
            ret.Add(removeDuplicates(t, threshold));
        }

        return ret;
    }
    
    private static Path64 pRemoveDuplicates(Path64 source, double threshold = Double.Epsilon)
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
    
    public static Paths64 stripColinear(Paths64 source, double angularTolerance = 0.0f)
    {
        return pStripColinear(source, angularTolerance);
    }

    private static Paths64 pStripColinear(Paths64 source, double angularTolerance = 0.0f)
    {
        int sLength = source.Count;
        Paths64 ret = Helper.initedPaths64(sLength);
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

    public static Path64 stripColinear(Path64 source, double angularTolerance = 0.0f)
    {
        return pStripColinear(source, angularTolerance);
    }

    private static Path64 pStripColinear(Path64 source, double angularTolerance = 0.0f)
    {
        switch (source.Count)
        {
            case < 3:
                return new(source);
        }

        Path64 ret = new(source);

        for (int pt = 0; pt < source.Count; pt++)
        {
            Point64 interSection_A, interSection_B, interSection_C;
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
                    if (pt == source.Count - 1) // last point in the list
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
            if (pt != 0 && pt != source.Count - 1)
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