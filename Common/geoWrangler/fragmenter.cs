using Clipper2Lib;
using System;
using System.Linq;

namespace geoWrangler;
public class Fragmenter
{
    private double resolution;
    private double scaleFactor;

    private const double def_res = 1.0;
    private const double def_scale = 1.0;

    public Fragmenter(double res = def_res, double scaling = def_scale)
    {
        pFragmenter(res, scaling);
    }

    private void pFragmenter(double res = def_res, double scaling = def_scale)
    {
        resolution = res;
        scaleFactor = scaling;
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
    
    public PathD fragmentPath(PathD pointList, double res, double scaling)
    {
        resolution = res;
        scaleFactor = scaling;
        return pFragmentPath(pointList);
    }

    private PathD pFragmentPath(PathD pointList)
    {
        PathD returnList = new();
        for (int pt = 0; pt < pointList.Count; pt++)
        {
            returnList.Add(pointList[pt]);
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
        double x_Distance = pt2.x - pt1.x;
        double y_Distance = pt2.y - pt1.y;
        double fragmentCount = Math.Floor(GeoWrangler.distanceBetweenPoints(pt1, pt2) / resolution);

        double x_Step = x_Distance / fragmentCount;
        double y_Step = y_Distance / fragmentCount;

        // Avoid sending back the first and last vertices of the path.
        for (int i = 1; i < fragmentCount; i++)
        {
            returnList.Add(new (pt1.x + i * x_Step, pt1.y + i * y_Step, pt1.z));
        }
        return returnList;
    }
}