using System.Linq;
using Clipper2Lib;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static PathsD skeleton(PathsD source)
    {
        PathsD ret = [];
        ret.AddRange(source.Select(skeleton));

        return ret;
    }
    
    public static PathD skeleton(PathD source)
    {
        PathD working = new(source);
        working.Reverse();
        
        RayCast rc = new(working, source, 1000, false);
        PathsD rays = rc.getRays();

        ClipperD c = new();
        c.AddOpenSubject(rays);
        c.AddClip(source);

        PathsD unused = [];
        PathsD rays_c1 = [];

        c.Execute(ClipType.Intersection, FillRule.EvenOdd, unused, rays_c1);

        PathD median = [];
        foreach (PathD ray in rays_c1)
        {
            PointD v0 = new(ray[0]);
            PointD dist = distanceBetweenPoints_point(ray[1], ray[0]);
            PointD v1 = new(ray[0].x + (0.5 * dist.x), ray[0].y + (0.5 * dist.y));
            median.Add(v1);
        }

        median = removeDuplicates(median);

        // Seem to be getting a trailing entry (duplicate of the first). Not sure why.
        median = stripTerminators(median, false);
        return median;
    }
}