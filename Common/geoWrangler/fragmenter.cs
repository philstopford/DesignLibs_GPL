using Clipper2Lib;
using System;
using System.Linq;
using System.Threading.Tasks;

namespace geoWrangler;
public class Fragmenter
{
    private double resolution;

    private const double def_res = 1.0;

    public Fragmenter(double res = def_res)
    {
        pFragmenter(res);
    }

    private void pFragmenter(double res = def_res)
    {
        resolution = res;
    }

    public PathsD fragmentPaths(PathsD source)
    {
        return pFragmentPaths(source);
    }

    private PathsD pFragmentPaths(PathsD source)
    {
        switch (source.Count)
        {
            case 0:
                return new (source);
        }

        PathsD ret = new(source.Select(pFragmentPath));

        return ret;
    }

    public PathD fragmentPath(PathD pointList)
    {
        return pFragmentPath(pointList);
    }

    public PathsD refragmentPaths(PathsD pointList, double res)
    {
        return pRefragmentPaths(pointList, res);
    }

    private PathsD pRefragmentPaths(PathsD pointList, double res)
    {
        int count = pointList.Count;
        PathsD ret = Helper.initedPathsD(count);
#if !GWSINGLETHREADED
        Parallel.For(0, count, i =>
#else
                for (int i = 0; i < count; i++)
#endif
            {
                ret[i] = pRefragmentPath(pointList[i], res);
            }
#if !GWSINGLETHREADED
        );
#endif
        return ret;
    }

    public PathD refragmentPath(PathD pointList, double res)
    {
        return pRefragmentPath(pointList, res);
    }

    private PathD pRefragmentPath(PathD pointList, double res)
    {
        bool closed = GeoWrangler.distanceBetweenPoints(pointList[0], pointList[^1]) <= constants.tolerance;
        PathD ret = Clipper.TrimCollinear(pointList, constants.roundingDecimalPrecision, closed);
        ret = fragmentPath(ret, res);
        return ret;
    }

    public PathD fragmentPath(PathD pointList, double res)
    {
        resolution = res;
        return pFragmentPath(pointList);
    }

    private PathD pFragmentPath(PathD pointList)
    {
        PathD returnList = new();
        for (int pt = 0; pt < pointList.Count; pt++)
        {
            returnList.Add(new(pointList[pt]));
            if (pt == pointList.Count - 1)
            {
                continue;
            }

            // Fragment path doesn't return start and end points - just points between.
            PathD newPtsList = fragmentPath(pointList[pt], pointList[pt + 1]);
            returnList.AddRange(newPtsList);
        }
        return returnList;
    }

    // Internals that don't return start and end points.

    private PathD fragmentPath(PointD pt1, PointD pt2)
    {
        PathD returnList = new();
        double fragmentCount = Math.Floor(GeoWrangler.distanceBetweenPoints(pt1, pt2) / resolution);

        if (fragmentCount > 0)
        {
            double x_Distance = pt2.x - pt1.x;
            double y_Distance = pt2.y - pt1.y;
            double x_Step = x_Distance / fragmentCount;
            double y_Step = y_Distance / fragmentCount;

            // Avoid sending back the first and last vertices of the path.
            for (int i = 1; i < fragmentCount; i++)
            {
                returnList.Add(new(pt1.x + i * x_Step, pt1.y + i * y_Step, pt1.z));
            }
        }

        return returnList;
    }
}