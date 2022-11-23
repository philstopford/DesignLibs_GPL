using System;
using Clipper2Lib;
using System.Threading.Tasks;
using utility;

namespace geoWrangler;
public static partial class GeoWrangler
{
    public static Paths64 extendEdges(Paths64 edges, double sizing)
    {
        return pExtendEdges(edges, sizing);
    }

    private static Paths64 pExtendEdges(Paths64 edges_, double sizing)
    {
        Paths64 edges = new(edges_);
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

    public static Path64 extendEdge(Path64 edge, double sizing)
    {
        return pExtendEdge(edge, sizing);
    }
    private static Path64 pExtendEdge(Path64 edge, double sizing)
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
    
    public static Path64 inflatePath(Path64 source, int width)
    {
        return pInflatePath(source, width);
    }

    private static Path64 pInflatePath(Path64 source, int width)
    {
        switch (width)
        {
            case 0:
                return source;
        }
        
        ClipperOffset co = new() {PreserveCollinear = true};
        Path64 a = new(source);
        // Path from Point auto-closes the input for historical reasons. We may not want this....
        if (pDistanceBetweenPoints(source[0], source[^1]) > Double.Epsilon)
        {
            pStripTerminators(a, false);
        }
        co.AddPath(a, JoinType.Miter, EndType.Square);
        Paths64 output = co.Execute(width);

        output = pReorderXY(output);

        return pClose(output[0]);
    }
    
}