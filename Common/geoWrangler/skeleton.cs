using System;
using System.Collections.Generic;
using Clipper2Lib;

namespace geoWrangler;

public static partial class GeoWrangler
{
    private static void ZFillCallback(PointD bot1, PointD top1, PointD bot2, PointD top2, ref PointD pt)
    {
        if (pt.z > -100)
        {
            pt.z = -1; // Tag our intersection points.
        }
    }

    public static PathsD skeleton(PathsD source)
    {
        PathsD ret = new();
        foreach (PathD s in source)
        {
            ret.Add(skeleton(s));
        }

        return ret;
    }
    
    public static PathD skeleton(PathD source_)
    {
        PathD source = new(source_);
        PathD working = new(source);
        working.Reverse();

        const int multirays = 1;
        const int totalrays_per_point = 1 + (2 * multirays);
        
        RayCast rc = new(working, source, 1000, false, multisampleRayCount: multirays);
        PathsD rays = rc.getRays();

        ClipperD c = new();
        c.ZCallback = ZFillCallback;
        c.AddOpenSubject(rays);
        c.AddClip(source);

        PathsD unused = new();
        PathsD rays_clipped = new();

        c.Execute(ClipType.Intersection, FillRule.EvenOdd, unused, rays_clipped);

        rays_clipped.RemoveAt(rays_clipped.Count - 1);
        PathsD rays_c1 = new();
        
        foreach (PathD ray in rays_clipped)
        {
            // Get source index from ray
            int index = (int)ray[0].z;
            if (index >= -1)
            {
                index = (int)ray[1].z;
            }

            // Invalid ray - probably due to fully overlapping with an edge. Reject this ray and move on.
            if (index > -100)
            {
                continue;
            }

            index = index + 100;
            // Working is reversed from source, so let's compensate.
            index = (source.Count-1) + index;

            PathD ray_corrected = new();
            ray_corrected.Add(new (ray[0].x, ray[0].y, index));
            ray_corrected.Add(new (ray[1].x, ray[1].y, index));
                
            // Do we need to reverse the ray?
            if ((ray_corrected[1].x == source[index].x) && (ray_corrected[1].y == source[index].y))
            {
                ray_corrected.Reverse();
            }
            
            rays_c1.Add(ray_corrected);
        }

        Dictionary<string, PointD> initial = new();
        Dictionary<string, PathD> displacements = new();
        Dictionary<string, PathD> displacements_by_initial_xy = new();
        
        foreach (PathD ray in rays_c1)
        {
            string key = "P" + ray[0].z;
            string key_xy = "X" + ray[0].x + "Y" + ray[0].y;

            if (!initial.ContainsKey(key))
            {
                initial.Add(key, new(ray[0]));
            }
            PointD dist = distanceBetweenPoints_point(ray[1], ray[0]);
            if (!displacements.ContainsKey(key))
            {
                displacements.Add(key, new PathD());
            }
            displacements[key].Add(dist);
            if (!displacements_by_initial_xy.ContainsKey(key_xy))
            {
                displacements_by_initial_xy.Add(key_xy, new PathD());
            }
            displacements_by_initial_xy[key_xy].Add(dist);
        }

        PathD median = new();

        for (int i = 0; i < initial.Count; i++)
        {
            string key = "P" + i;
            PointD start = new(initial[key]);
            string key_xy = "X" + start.x + "Y" + start.y;
            try
            {

                // Figure out our average displacement.
                PointD displacement = new(0, 0);
                //int displacement_count = displacements[key].Count;
                int displacement_count = displacements_by_initial_xy[key_xy].Count;
                PathD displacements_path = displacements_by_initial_xy[key_xy];
                for (int d = 0; d < displacement_count; d++)
                {
                    displacement.x += displacements_path[d].x;
                    displacement.y += displacements_path[d].y;
                }

                displacement.x /= totalrays_per_point;
                displacement.y /= totalrays_per_point;

                median.Add(new(start.x + (0.5 * displacement.x), start.y + (0.5 * displacement.y)));

            }
            catch (Exception e)
            {
            }
        }

        median = removeDuplicates(median);

        // Seem to be getting a trailing entry (duplicate of the first). Not sure why.
        median = stripTerminators(median, false);
        return median;
    }
}