using System;
using Clipper2Lib;
using geoLib;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using utility;

namespace geoWrangler;

using Path = List<Point64>;
using Paths = List<List<Point64>>;

public static partial class GeoWrangler
{
    public static Paths extendEdges(Paths edges, double sizing)
    {
        return pExtendEdges(edges, sizing);
    }

    private static Paths pExtendEdges(Paths edges, double sizing)
    {
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

    public static Path extendEdge(Path edge, double sizing)
    {
        return pExtendEdge(edge, sizing);
    }
    private static Path pExtendEdge(Path edge, double sizing)
    {
        // Get sorted out for dx, dy and normalization.
        double dx = edge[0].X - edge[1].X;
        double dy = edge[0].Y - edge[1].Y;

        double length = Math.Sqrt(Utils.myPow(dx, 2) + Utils.myPow(dy, 2));

        dx /= length;
        dy /= length;
        
        // Extend the line slightly.
        double edge0_newX = edge[0].X;
        double edge0_newY = edge[0].Y;
        double edge1_newX = edge[1].X;
        double edge1_newY = edge[1].Y;

        // Move the edge according to the keyhole sizing, to extend it. Then add 1 to ensure an overlap.
        double X_shift = (long) Math.Abs(dx * sizing) + 1;
        double Y_shift = (long) Math.Abs(dy * sizing) + 1;
        if (edge[0].X <= edge[1].X)
        {
            edge0_newX -= X_shift;
            edge1_newX += X_shift;
        }
        else
        {
            edge0_newX += X_shift;
            edge1_newX -= X_shift;
        }
        
        if (edge[0].Y <= edge[1].Y)
        {
            edge0_newY -= Y_shift;
            edge1_newY += Y_shift;
        }
        else
        {
            edge0_newY += Y_shift;
            edge1_newY -= Y_shift;
        }
        
        edge[0] = new Point64(edge0_newX, edge0_newY);
        edge[1] = new Point64(edge1_newX, edge1_newY);

        return edge;
    }
    public static GeoLibPoint[] inflatePath(GeoLibPoint[] source, int width)
    {
        return pInflatePath(source, width);
    }

    private static GeoLibPoint[] pInflatePath(GeoLibPoint[] source, int width)
    {
        switch (width)
        {
            case 0:
                return source;
        }
        
        ClipperOffset co = new();
        Path a = GeoWrangler.pathFromPoint(source, 1);
        // Path from Point auto-closes the input for historical reasons. We may not want this....
        if (pDistanceBetweenPoints(source[0], source[^1]) > Double.Epsilon)
        {
            pStripTerminators(a, false);
        }
        co.AddPath(a, JoinType.Miter, EndType.Square);
        Paths output = ClipperFunc.Paths64(co.Execute(width));

        output = pReorder(output);

        return pPointFromPath(pClose(output[0]), 1);

    }

    public static GeoLibPointF[] resize(GeoLibPointF[] source, double factor)
    {
        return pResize(source, factor);
    }

    private static GeoLibPointF[] pResize(GeoLibPointF[] source, double factor)
    {
        int sLength = source.Length;
        GeoLibPointF[] ret = new GeoLibPointF[sLength];
#if !GWSINGLETHREADED
        Parallel.For(0, sLength, pt =>
#else
            for (int pt = 0; pt < sLength; pt++)
#endif
            {
                ret[pt] = new GeoLibPointF(source[pt].X * factor, source[pt].Y * factor);
            }
#if !GWSINGLETHREADED
        );
#endif
        return ret;
    }

    public static GeoLibPoint[] resize_to_int(GeoLibPointF[] source, double factor)
    {
        return pResize_to_int(source, factor);
    }

    private static GeoLibPoint[] pResize_to_int(GeoLibPointF[] source, double factor)
    {
        int sLength = source.Length;
        GeoLibPoint[] ret = new GeoLibPoint[sLength];

#if !GWSINGLETHREADED
        Parallel.For(0, sLength, pt =>
#else
            for (int pt = 0; pt < sLength; pt++)
#endif
            {
                ret[pt] = new GeoLibPoint((int)(source[pt].X * factor), (int)(source[pt].Y * factor));
            }
#if !GWSINGLETHREADED
        );
#endif
        return ret;
    }

    public static GeoLibPoint[] resize(GeoLibPoint[] source, double factor)
    {
        return pResize(source, factor);
    }

    private static GeoLibPoint[] pResize(GeoLibPoint[] source, double factor)
    {
        int sLength = source.Length;
        GeoLibPoint[] ret = new GeoLibPoint[sLength];
#if !GWSINGLETHREADED
        Parallel.For(0, sLength, pt =>
#else
            for (int pt = 0; pt < sLength; pt++)
#endif
            {
                ret[pt] = new GeoLibPoint(source[pt].X * factor, source[pt].Y * factor);
            }
#if !GWSINGLETHREADED
        );
#endif
        return ret;
    }

    public static GeoLibPoint[] resize(GeoLibPoint pivot, GeoLibPoint[] source, double factor)
    {
        return pResize(pivot, source, factor);
    }

    private static GeoLibPoint[] pResize(GeoLibPoint pivot, GeoLibPoint[] source, double factor)
    {
        GeoLibPoint[] pointarray = new GeoLibPoint[source.Length];
#if !GWSINGLETHREADED
        Parallel.For(0, pointarray.Length, i => 
#else
            for (int i = 0; i < pointarray.Length; i++)
#endif
            {
                pointarray[i] = new GeoLibPoint(pivot.X + (source[i].X - pivot.X) * factor, pivot.Y + (source[i].Y - pivot.Y) * factor);
            }
#if !GWSINGLETHREADED
        );
#endif
        return pointarray;
    }
}