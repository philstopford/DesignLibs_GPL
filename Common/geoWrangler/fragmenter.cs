using Clipper2Lib;
using geoLib;
using System;
using System.Collections.Generic;
using System.Linq;

namespace geoWrangler;
public class Fragmenter
{
    private double resolution;
    private double scaleFactor;

    private const double def_res = 1.0;
    private const double def_scale = 1E4;

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
    
    public Paths64 fragmentPaths(Paths64 source)
    {
        return pFragmentPaths(source);
    }

    private Paths64 pFragmentPaths(Paths64 source)
    {
        return new (source.Select(pFragmentPath));
    }

    public Path64 fragmentPath(Path64 source)
    {
        return pFragmentPath(source);
    }

    private Path64 pFragmentPath(Path64 source)
    {
        Path64 ret = new();
        bool closed = source[0].X == source[^1].X && source[0].Y == source[^1].Y;
        // Decouple the geometry.
        Path64 t = new();
        for (int p = 0; p < source.Count; p++)
        {
            t.Add(new Point64(source[p].X, source[p].Y, source[p].Z));
        }

        t = closed switch
        {
            true => GeoWrangler.stripTerminators(t, false),
            _ => t
        };

        for (int p = 0; p < t.Count - 1; p++)
        {
            ret.Add(new Point64(t[p]));
            ret.AddRange(pFragmentPath(t[p], t[p + 1]));
        }

        switch (closed)
        {
            case true:
                ret.Add(new Point64(t[^1]));
                ret.AddRange(pFragmentPath(t[^1], t[0]));
                ret = GeoWrangler.close(ret);
                break;
        }

        return ret;
    }

    public Path64 fragmentPath(Point64 pt1, Point64 pt2)
    {
        return pFragmentPath(pt1, pt2);
    }

    private Path64 pFragmentPath(Point64 pt1, Point64 pt2, bool startAndEndPoints = false)
    {
        Path64 returnPath = new();
        long x_Distance = pt2.X - pt1.X;
        long y_Distance = pt2.Y - pt1.Y;
        // We do some manipulation here because the points in the call are scaled by CentralProperties.scaleFactorForOperation.
        // That would lead to too many points being injected. We therefore scale resolution by 1E4
        int fragmentCount = Convert.ToInt32(Math.Floor(GeoWrangler.distanceBetweenPoints(pt1, pt2) / (resolution * (scaleFactor / 1E4))));

        switch (startAndEndPoints)
        {
            case true:
                returnPath.Add(new Point64(pt1));
                break;
        }
        switch (fragmentCount)
        {
            case > 0:
            {
                double x_Step = Convert.ToDouble(x_Distance) / fragmentCount;
                double y_Step = Convert.ToDouble(y_Distance) / fragmentCount;

                for (int i = 1; i < fragmentCount; i++)
                {
                    returnPath.Add(new Point64(pt1.X + i * x_Step, pt1.Y + i * y_Step, pt1.Z));
                }

                break;
            }
        }
        switch (startAndEndPoints)
        {
            case true:
                returnPath.Add(new Point64(pt2));
                break;
        }
        return returnPath;
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