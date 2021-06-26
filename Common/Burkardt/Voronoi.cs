using System;
using Burkardt.Types;

namespace Burkardt
{
    public static class Voronoi
    {
        public static double voronoi_polygon_area(int node, int neighbor_num,
            int[] neighbor_index, int node_num, double[] node_xy )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VORONOI_POLYGON_AREA computes the area of a Voronoi polygon.
        //
        //  Formula:
        //
        //    It is assumed that the Voronoi polygon is finite!  Every Voronoi
        //    diagram includes some regions which are infinite, and for those,
        //    this formula is not appropriate.
        //
        //    The routine is given the indices of the nodes that are
        //    Voronoi "neighbors" of a given node.  These are also the nodes
        //    that are paired to form edges of Delaunay triangles.
        //
        //    The assumption that the polygon is a Voronoi polygon is
        //    used to determine the location of the boundaries of the polygon,
        //    which are the perpendicular bisectors of the lines connecting
        //    the center point to each of its neighbors.
        //
        //    The finiteness assumption is employed in part in the
        //    assumption that the polygon is bounded by the finite
        //    line segments from point 1 to 2, 2 to 3, ...,
        //    M-1 to M, and M to 1, where M is the number of neighbors.
        //
        //    It is assumed that this routine is being called by a
        //    process which has computed the Voronoi diagram of a large
        //    set of nodes, so the arrays X and Y are dimensioned by
        //    NODE_NUM, which may be much greater than the number of neighbor
        //    nodes.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 February 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Atsuyuki Okabe, Barry Boots, Kokichi Sugihara, Sung Nok Chiu,
        //    Spatial Tessellations: Concepts and Applications of Voronoi Diagrams,
        //    Second Edition,
        //    Wiley, 2000, page 485.
        //
        //  Parameters:
        //
        //    Input, int NODE, the index of the node whose Voronoi
        //    polygon is to be measured. 0 <= NODE < NODE_NUM.
        //
        //    Input, int NEIGHBOR_NUM, the number of neighbor nodes of
        //    the given node.
        //
        //    Input, int NEIGHBOR_INDEX[NEIGHBOR_NUM], the indices
        //    of the neighbor nodes (used to access X and Y).  The neighbor
        //    nodes should be listed in the (counter-clockwise) order in
        //    which they occur as one circles the center node.
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
        //
        //    Output, double VORONOI_POLYGON_AREA, the area of the Voronoi polygon.
        //
        {
            double a;
            double area;
            double b;
            double c;
            int i;
            int ip1;
            double ui;
            double uip1;
            double vi;
            double vip1;
            double xc;
            double xi;
            double xip1;
            double yc;
            double yi;
            double yip1;

            area = 0.0;

            if (node < 0 || node_num <= node)
            {
                Console.WriteLine("");
                Console.WriteLine("  VORONOI_POLYGON_AREA - Fatal error!");
                Console.WriteLine("  Illegal value of input parameter NODE.");
                return 1;
            }

            xc = node_xy[0 + node * 2];
            yc = node_xy[1 + node * 2];

            i = neighbor_num - 1;
            i = neighbor_index[i];

            xi = node_xy[0 + i * 2];
            yi = node_xy[1 + i * 2];

            ip1 = 0;
            ip1 = neighbor_index[ip1];

            xip1 = node_xy[0 + ip1 * 2];
            yip1 = node_xy[1 + ip1 * 2];
            a = (xi * xi + yi * yi - xc * xc - yc * yc);
            b = (xip1 * xip1 + yip1 * yip1 - xc * xc - yc * yc);
            c = 2.0 * ((xi - xc) * (yip1 - yc) - (xip1 - xc) * (yi - yc));
            uip1 = (a * (yip1 - yc) - b * (yi - yc)) / c;
            vip1 = (a * (xip1 - xc) - b * (xi - xc)) / c;

            for (i = 0; i < neighbor_num; i++)
            {
                xi = xip1;
                yi = yip1;
                ui = uip1;
                vi = vip1;

                ip1 = typeMethods.i4_wrap(i + 1, 0, neighbor_num - 1);
                ip1 = neighbor_index[ip1];

                xip1 = node_xy[0 + ip1 * 2];
                yip1 = node_xy[1 + ip1 * 2];
                a = (xi * xi + yi * yi - xc * xc - yc * yc);
                b = (xip1 * xip1 + yip1 * yip1 - xc * xc - yc * yc);
                c = 2.0 * ((xi - xc) * (yip1 - yc) - (xip1 - xc) * (yi - yc));
                uip1 = (a * (yip1 - yc) - b * (yi - yc)) / c;
                vip1 = (a * (xip1 - xc) - b * (xi - xc)) / c;

                area = area + uip1 * vi - ui * vip1;
            }

            area = 0.5 * area;

            return area;
        }

        public static double[] voronoi_polygon_centroid(int node, int neighbor_num,
            int[] neighbor_index, int node_num, double[] node_xy )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VORONOI_POLYGON_CENTROID_2D computes the centroid of a Voronoi polygon.
        //
        //  Formula:
        //
        //    It is assumed that the Voronoi polygon is finite!  Every Voronoi
        //    diagram includes some regions which are infinite, and for those,
        //    this formula is not appropriate.
        //
        //    The routine is given the indices of the nodes that are
        //    Voronoi "neighbors" of a given node.  These are also the nodes
        //    that are paired to form edges of Delaunay triangles.
        //
        //    The assumption that the polygon is a Voronoi polygon is
        //    used to determine the location of the boundaries of the polygon,
        //    which are the perpendicular bisectors of the lines connecting
        //    the center point to each of its neighbors.
        //
        //    The finiteness assumption is employed in part in the
        //    assumption that the polygon is bounded by the finite
        //    line segments from point 1 to 2, 2 to 3, ...,
        //    M-1 to M, and M to 1, where M is the number of neighbors.
        //
        //    It is assumed that this routine is being called by a
        //    process which has computed the Voronoi diagram of a large
        //    set of nodes, so the arrays X and Y are dimensioned by
        //    NODE_NUM, which may be much greater than the number of neighbor
        //    nodes.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 February 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Atsuyuki Okabe, Barry Boots, Kokichi Sugihara, Sung Nok Chiu,
        //    Spatial Tessellations: Concepts and Applications of Voronoi Diagrams,
        //    Second Edition,
        //    Wiley, 2000, page 490.
        //
        //  Parameters:
        //
        //    Input, int NODE, the index of the node whose Voronoi
        //    polygon is to be analyzed.  1 <= NODE <= NODE_NUM.
        //
        //    Input, int NEIGHBOR_NUM, the number of neighbor nodes of
        //    the given node.
        //
        //    Input, int NEIGHBOR_INDEX[NEIGHBOR_NUM], the indices
        //    of the neighbor nodes.  These indices are used to access the
        //    X and Y arrays.  The neighbor nodes should be listed in the
        //    (counter-clockwise) order in which they occur as one circles
        //    the center node.
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
        //
        //    Output, double *VORONOI_POLYGON_CENTROID_2D, a pointer to a 2D array
        //    containing the coordinates of the centroid of the Voronoi polygon
        //    of node NODE.
        //
        {
            double a;
            double area;
            double b;
            double c;
            double[] centroid;
            int i;
            int ip1;
            double ui;
            double uip1;
            double vi;
            double vip1;
            double xc;
            double xi;
            double xip1;
            double yc;
            double yi;
            double yip1;

            centroid = new double[2];

            centroid[0] = 0.0;
            centroid[1] = 0.0;

            if (node < 0 || node_num <= node)
            {
                Console.WriteLine("");
                Console.WriteLine("VORONOI_POLYGON_CENTROID - Fatal error!");
                Console.WriteLine("  Illegal value of input parameter NODE.");
                return null;
            }

            xc = node_xy[0 + node * 2];
            yc = node_xy[1 + node * 2];

            i = neighbor_num - 1;
            i = neighbor_index[i];

            xi = node_xy[0 + i * 2];
            yi = node_xy[1 + i * 2];

            ip1 = 0;
            ip1 = neighbor_index[ip1];

            xip1 = node_xy[0 + ip1 * 2];
            yip1 = node_xy[1 + ip1 * 2];
            a = (xi * xi + yi * yi - xc * xc - yc * yc);
            b = (xip1 * xip1 + yip1 * yip1 - xc * xc - yc * yc);
            c = 2.0 * ((xi - xc) * (yip1 - yc) - (xip1 - xc) * (yi - yc));
            uip1 = (a * (yip1 - yc) - b * (yi - yc)) / c;
            vip1 = (a * (xip1 - xc) - b * (xi - xc)) / c;

            for (i = 0; i < neighbor_num; i++)
            {
                xi = xip1;
                yi = yip1;
                ui = uip1;
                vi = vip1;

                ip1 = typeMethods.i4_wrap(i + 1, 0, neighbor_num - 1);
                ip1 = neighbor_index[ip1];

                xip1 = node_xy[0 + ip1 * 2];
                yip1 = node_xy[1 + ip1 * 2];
                a = (xi * xi + yi * yi - xc * xc - yc * yc);
                b = (xip1 * xip1 + yip1 * yip1 - xc * xc - yc * yc);
                c = 2.0 * ((xi - xc) * (yip1 - yc) - (xip1 - xc) * (yi - yc));
                uip1 = (a * (yip1 - yc) - b * (yi - yc)) / c;
                vip1 = (a * (xip1 - xc) - b * (xi - xc)) / c;

                centroid[0] = centroid[0] + (vi - vip1)
                    * ((uip1 + ui) * (uip1 + ui) - uip1 * ui);
                centroid[1] = centroid[1] + (ui - uip1)
                    * ((vip1 + vi) * (vip1 + vi) - vip1 * vi);
            }

            area = voronoi_polygon_area(node, neighbor_num, neighbor_index,
                node_num, node_xy);

            centroid[0] = centroid[0] / (6.0 * area);
            centroid[1] = centroid[1] / (6.0 * area);

            return centroid;
        }

        public static void voronoi_polygon_vertices(int node, int neighbor_num,
            int[] neighbor_index, int node_num, double[] node_xy, ref double[] v )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VORONOI_POLYGON_VERTICES_2D computes the vertices of a Voronoi polygon.
        //
        //  Formula:
        //
        //    This routine is only appropriate for Voronoi polygons that are finite.
        //
        //    The routine is given the indices of the nodes that are neighbors of a
        //    given "center" node.  A node is a neighbor of the center node if the
        //    Voronoi polygons of the two nodes share an edge.  The triangles of the
        //    Delaunay triangulation are formed from successive pairs of these neighbor
        //    nodes along with the center node.
        //
        //    Given only the neighbor node information, it is possible to determine
        //    the location of the vertices of the polygonal Voronoi region by computing
        //    the circumcenters of the Delaunay triangles.
        //
        //    It is assumed that this routine is being called by a process which has
        //    computed the Voronoi diagram of a large set of nodes, so the arrays X and
        //    Y are dimensioned by NODE_NUM, which may be much greater than the number
        //    of neighbor nodes.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 February 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Atsuyuki Okabe, Barry Boots, Kokichi Sugihara, Sung Nok Chiu,
        //    Spatial Tessellations: Concepts and Applications of Voronoi Diagrams,
        //    Second Edition,
        //    Wiley, 2000.
        //
        //  Parameters:
        //
        //    Input, int NODE, the index of the node whose Voronoi
        //    polygon is to be analyzed.  1 <= NODE <= NODE_NUM.
        //
        //    Input, int NEIGHBOR_NUM, the number of neighbor nodes of
        //    the given node.
        //
        //    Input, int NEIGHBOR_INDEX(NEIGHBOR_NUM), the indices
        //    of the neighbor nodes.  These indices are used to access the
        //    X and Y arrays.  The neighbor nodes should be listed in the
        //    (counter-clockwise) order in which they occur as one circles
        //    the center node.
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
        //
        //    Output, double V[2*NEIGHBOR_NUM], the vertices of the Voronoi polygon
        //    around node NODE.
        //
        {
            int DIM_NUM = 2;

            double[] center;
            int i;
            int ip1;
            double[] t = new double[DIM_NUM * 3];

            if (node < 0 || node_num <= node)
            {
                Console.WriteLine("");
                Console.WriteLine("VORONOI_POLYGON_VERTICES - Fatal error!");
                Console.WriteLine("  Illegal value of input parameter NODE.");
                return;
            }

            t[0 + 0 * 2] = node_xy[0 + node * 2];
            t[1 + 0 * 2] = node_xy[1 + node * 2];

            ip1 = neighbor_index[0];
            t[0 + 2 * 2] = node_xy[0 + ip1 * 2];
            t[1 + 2 * 2] = node_xy[1 + ip1 * 2];

            for (i = 0; i < neighbor_num; i++)
            {
                t[0 + 1 * 2] = t[0 + 2 * 2];
                t[1 + 1 * 2] = t[1 + 2 * 2];

                ip1 = typeMethods.i4_wrap(i + 1, 0, neighbor_num - 1);
                ip1 = neighbor_index[ip1];

                t[0 + 2 * 2] = node_xy[0 + ip1 * 2];
                t[1 + 2 * 2] = node_xy[1 + ip1 * 2];

                center = typeMethods.triangle_circumcenter_2d(t);

                v[0 + i * 2] = center[0];
                v[1 + i * 2] = center[1];

            }
        }

    }
}