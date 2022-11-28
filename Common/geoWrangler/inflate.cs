using System;
using Clipper2Lib;
using System.Threading.Tasks;
using utility;

namespace geoWrangler;

public static partial class GeoWrangler
{
    public static PathsD extendEdges(PathsD edges, double sizing)
    {
        return pExtendEdges(edges, sizing);
    }

    private static PathsD pExtendEdges(PathsD edges_, double sizing)
    {
        PathsD edges = new(edges_);
        int sLength = edges.Count;
#if !GWSINGLETHREADED
        Parallel.For(0, sLength, i =>
#else
            for (int i = 0; i < sLength; i++)
#endif
            {
                edges[i] = pExtendEdge(edges[i], sizing);
            }
#if !GWSINGLETHREADED
        );
#endif
        return edges;
    }

    public static PathD extendEdge(PathD edge, double sizing)
    {
        return pExtendEdge(edge, sizing);
    }
    private static PathD pExtendEdge(PathD edge, double sizing)
    {
        // Get sorted out for dx, dy and normalization.
        double dx = edge[0].x - edge[1].x;
        double dy = edge[0].y - edge[1].y;

        double length = Math.Sqrt(Utils.myPow(dx, 2) + Utils.myPow(dy, 2));

        dx /= length;
        dy /= length;
        
        // Extend the line slightly.
        double edge0_newX = edge[0].x;
        double edge0_newY = edge[0].y;
        double edge1_newX = edge[1].x;
        double edge1_newY = edge[1].y;

        // Move the edge according to the keyhole sizing, to extend it. Then add 1 to ensure an overlap.
        double X_shift = Math.Abs(dx * sizing) + 1;
        double Y_shift = Math.Abs(dy * sizing) + 1;
        if (edge[1].x - edge[0].x > constants.tolerance )
        {
            edge0_newX -= X_shift;
            edge1_newX += X_shift;
        }
        else
        {
            edge0_newX += X_shift;
            edge1_newX -= X_shift;
        }
        
        if (edge[1].y - edge[0].y > constants.tolerance )
        {
            edge0_newY -= Y_shift;
            edge1_newY += Y_shift;
        }
        else
        {
            edge0_newY += Y_shift;
            edge1_newY -= Y_shift;
        }
        
        edge[0] = new (edge0_newX, edge0_newY, edge[0].z);
        edge[1] = new (edge1_newX, edge1_newY, edge[1].z);

        return edge;
    }
    public static PathD inflatePath(PathD source, int width)
    {
        return pInflatePath(source, width);
    }

    private static PathD pInflatePath(PathD source, int width)
    {
        switch (width)
        {
            case 0:
                return new(source);
        }

        
        // Need to workaround missing PathD support in ClipperOffset...
        double scalar = 1000;
        double scalar_inv = 1.0 / scalar;
        Path64 rescaledSource = _pPath64FromPathD(source, scalar);

        ClipperOffset co = new() {PreserveCollinear = true};
        // Path from Point auto-closes the input for historical reasons. We may not want this....
        if (pDistanceBetweenPoints(rescaledSource[0], rescaledSource[^1]) > constants.tolerance)
        {
            rescaledSource = pStripTerminators(rescaledSource, false);
        }
        co.AddPath(rescaledSource, JoinType.Miter, EndType.Square);

        Paths64 output = co.Execute(width);
        PathD ret = _pPathDFromPath64(output[0], scalar_inv);
        ret = pReorderXY(ret);
        return pClose(ret);
    }

    public static PathsD resize(PathsD source, double factor)
    {
        return pResize(source, factor);
    }

    private static PathsD pResize(PathsD source, double factor)
    {
        PathsD ret = new();
        foreach (PathD p in source)
        {
            ret.Add(pResize(p, factor));
        }

        return ret;
    }

    public static PathD resize(PathD source, double factor)
    {
        return pResize(source, factor);
    }

    private static PathD pResize(PathD source, double factor)
    {
        int sLength = source.Count;
        PathD ret = Helper.initedPathD(sLength);
#if !GWSINGLETHREADED
        Parallel.For(0, sLength, pt =>
#else
            for (int pt = 0; pt < sLength; pt++)
#endif
            {
                ret[pt] = new (source[pt].x * factor, source[pt].y * factor, source[pt].z);
            }
#if !GWSINGLETHREADED
        );
#endif
        return ret;
    }

    public static Path64 resize_to_int(PathD source, double factor)
    {
        return pResize_to_int(source, factor);
    }

    private static Path64 pResize_to_int(PathD source, double factor)
    {
        int sLength = source.Count;
        Path64 ret = Helper.initedPath64(sLength);

#if !GWSINGLETHREADED
        Parallel.For(0, sLength, pt =>
#else
            for (int pt = 0; pt < sLength; pt++)
#endif
            {
                ret[pt] = new ((long)(source[pt].x * factor), (long)(source[pt].y * factor), (long)source[pt].z);
            }
#if !GWSINGLETHREADED
        );
#endif
        return ret;
    }

    public static Path64 resize(Path64 source, double factor)
    {
        return pResize(source, factor);
    }

    private static Path64 pResize(Path64 source, double factor)
    {
        int sLength = source.Count;
        Path64 ret = Helper.initedPath64(sLength);
#if !GWSINGLETHREADED
        Parallel.For(0, sLength, pt =>
#else
            for (int pt = 0; pt < sLength; pt++)
#endif
            {
                ret[pt] = new (source[pt].X * factor, source[pt].Y * factor, source[pt].Z);
            }
#if !GWSINGLETHREADED
        );
#endif
        return ret;
    }

    public static Path64 resize(Point64 pivot, Path64 source, double factor)
    {
        return pResize(pivot, source, factor);
    }

    private static Path64 pResize(Point64 pivot, Path64 source, double factor)
    {
        Path64 pointarray = Helper.initedPath64(source.Count);
#if !GWSINGLETHREADED
        Parallel.For(0, pointarray.Count, i => 
#else
            for (int i = 0; i < pointarray.Length; i++)
#endif
            {
                pointarray[i] = new (pivot.X + (source[i].X - pivot.X) * factor, pivot.Y + (source[i].Y - pivot.Y) * factor, source[i].Z);
            }
#if !GWSINGLETHREADED
        );
#endif
        return pointarray;
    }
}