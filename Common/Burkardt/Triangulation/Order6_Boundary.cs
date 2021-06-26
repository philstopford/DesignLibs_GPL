using System;
using Burkardt.Types;

namespace Burkardt.TriangulationNS
{
    public static partial class Boundary
    {
        public static int triangulation_order6_boundary_edge_count(int triangle_num,
                int[] triangle_node)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGULATION_ORDER6_BOUNDARY_EDGE_COUNT counts the boundary edges.
            //
            //  Discussion:
            //
            //    This routine is given a triangulation, a set of 6-node triangles.
            //    It is assumed that, in each list of 6 nodes, the vertices are listed
            //    first, in counterclockwise order, followed by the three midside nodes,
            //    in counterclockwise order, starting with the node between vertices
            //    1 and 2.
            //
            //    It is assumed that each edge of the triangulation is either
            //    * an INTERIOR edge, which is listed twice, once with positive
            //      orientation and once with negative orientation, or;
            //    * a BOUNDARY edge, which will occur only once.
            //
            //    This routine should work even if the region has holes - as long
            //    as the boundary of the hole comprises more than 3 edges!
            //
            //    Except for the dimension of TRIANGLE, this routine is identical
            //    to the routine for the order 3 case.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 June 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int TRIANGLE_NUM, the number of triangles.
            //
            //    Input, int TRIANGLE_NODE[6*TRIANGLE_NUM], the nodes that make up the
            //    triangles.  These should be listed in counterclockwise order.
            //
            //    Output, integer TRIANGULATION_ORDER6_BOUNDARY_EDGE_COUNT, the number
            //    of boundary edges.
            //
        {
            int boundary_edge_num;
            int e1;
            int e2;
            int[] edge;
            int interior_edge_num;
            int j;
            int m;
            int n;
            int unique_num;

            m = 2;
            n = 3 * triangle_num;
            //
            //  Set up the edge array.
            //
            edge = new int[m * n];

            for (j = 0; j < triangle_num; j++)
            {
                edge[0 + (j) * m] = triangle_node[0 + j * 6];
                edge[1 + (j) * m] = triangle_node[1 + j * 6];
                edge[0 + (j + triangle_num) * m] = triangle_node[1 + j * 6];
                edge[1 + (j + triangle_num) * m] = triangle_node[2 + j * 6];
                edge[0 + (j + 2 * triangle_num) * m] = triangle_node[2 + j * 6];
                edge[1 + (j + 2 * triangle_num) * m] = triangle_node[0 + j * 6];
            }

            //
            //  In each column, force the smaller entry to appear first.
            //
            for (j = 0; j < n; j++)
            {
                e1 = Math.Min(edge[0 + j * m], edge[1 + j * m]);
                e2 = Math.Max(edge[0 + j * m], edge[1 + j * m]);
                edge[0 + j * m] = e1;
                edge[1 + j * m] = e2;
            }

            //
            //  Ascending sort the column array.
            //
            typeMethods.i4col_sort_a(m, n, ref edge);
            //
            //  Get the number of unique columns in EDGE.
            //
            unique_num = typeMethods.i4col_sorted_unique_count(m, n, edge);

            interior_edge_num = 3 * triangle_num - unique_num;

            boundary_edge_num = 3 * triangle_num - 2 * interior_edge_num;
            
            return boundary_edge_num;
        }

        public static int triangulation_order6_boundary_edge_count_euler(int node_num,
                int triangle_num, int hole_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGULATION_ORDER6_BOUNDARY_EDGE_COUNT_EULER counts boundary edges.
            //
            //  Discussion:
            //
            //    We assume we are given information about an order 6 triangulation
            //    of a set of nodes in the plane.
            //
            //    By ignoring the midside nodes, we can determine the corresponding
            //    information for an order 3 triangulation, and apply
            //    Euler's formula to determine the number of edges that lie on the
            //    boundary of the set of nodes.
            //
            //    Thus, if we have TRIANGLE_NUM triangles, and NODE_NUM nodes, we
            //    imagine that each triangle is replaced by 4 triangles, created
            //    by adding the edges created by joining the midside nodes.
            //
            //    Thus, for 4 * TRIANGLE_NUM triangles, we can apply Euler's formula
            //    to compute the number of boundary edges.
            //
            //    Now, to adjust the data to our order 6 triangles, we divide the
            //    number of boundary edges by 2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 June 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Marc deBerg, Marc Krevald, Mark Overmars, Otfried Schwarzkopf,
            //    Computational Geometry,
            //    Springer, 2000,
            //    ISBN: 3-540-65620-0.
            //
            //  Parameters:
            //
            //    Input, integer NODE_NUM, the number of nodes.
            //
            //    Input, integer TRIANGLE_NUM, the number of triangles.
            //
            //    Input, integer HOLE_NUM, the number of internal nodes.
            //
            //    Output, int TRIANGULATION_ORDER6_BOUNDARY_EDGE_COUNT, the number of
            //    edges that lie on the boundary of the triangulation.
            //
        {
            int boundary_num;

            boundary_num = (2 * node_num + 2 * hole_num - 4 * triangle_num - 2) / 2;

            return boundary_num;
        }

        public static bool[] triangulation_order6_boundary_node(int node_num, int triangle_num,
                int[] triangle_node)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGULATION_ORDER6_BOUNDARY_NODE indicates nodes on the boundary.
            //
            //  Discussion:
            //
            //    This routine is given an order 6 triangulation, an abstract list of
            //    sets of six nodes.  The vertices are listed clockwise, then the
            //    midside nodes.
            //
            //    It is assumed that each edge of the triangulation is either
            //    * an INTERIOR edge, which is listed twice, once with positive
            //      orientation and once with negative orientation, or;
            //    * a BOUNDARY edge, which will occur only once.
            //
            //    This routine should work even if the region has holes - as long
            //    as the boundary of the hole comprises more than 3 edges!
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    25 January 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NODE_NUM, the number of nodes.
            //
            //    Input, int TRIANGLE_NUM, the number of triangles.
            //
            //    Input, int TRIANGLE_NODE[6*TRIANGLE_NUM], the nodes that make up the
            //    triangles.
            //
            //    Output, bool TRIANGULATION_ORDER6_BOUNDARY_NODE[NODE_NUM],
            //    is TRUE if the node is on a boundary edge.
            //
        {
            int e1;
            int e2;
            int[] edge;
            bool equal;
            int i;
            int j;
            int m;
            int n;
            bool[] node_boundary;

            m = 3;
            n = 3 * triangle_num;
            //
            //  Set up the edge array.
            //
            edge = new int[m * n];

            for (j = 0; j < triangle_num; j++)
            {
                edge[0 + (j) * m] = triangle_node[0 + j * 6];
                edge[1 + (j) * m] = triangle_node[3 + j * 6];
                edge[2 + (j) * m] = triangle_node[1 + j * 6];

                edge[0 + (j + triangle_num) * m] = triangle_node[1 + j * 6];
                edge[1 + (j + triangle_num) * m] = triangle_node[4 + j * 6];
                edge[2 + (j + triangle_num) * m] = triangle_node[2 + j * 6];

                edge[0 + (j + 2 * triangle_num) * m] = triangle_node[2 + j * 6];
                edge[1 + (j + 2 * triangle_num) * m] = triangle_node[5 + j * 6];
                edge[2 + (j + 2 * triangle_num) * m] = triangle_node[0 + j * 6];
            }

            //
            //  In each column, force the smaller entry to appear first.
            //
            for (j = 0; j < n; j++)
            {
                e1 = Math.Min(edge[0 + j * m], edge[2 + j * m]);
                e2 = Math.Max(edge[0 + j * m], edge[2 + j * m]);
                edge[0 + j * m] = e1;
                edge[2 + j * m] = e2;
            }

            //
            //  Ascending sort the column array.
            //
            typeMethods.i4col_sort_a(m, n, ref edge);
            //
            //  Records which appear twice are internal edges and can be ignored.
            //
            node_boundary = new bool[node_num];

            for (i = 0; i < node_num; i++)
            {
                node_boundary[i] = false;
            }

            j = 0;

            while (j < 3 * triangle_num)
            {
                j = j + 1;

                if (j == 3 * triangle_num)
                {
                    for (i = 0; i < m; i++)
                    {
                        node_boundary[edge[i + (j - 1) * m] - 1] = true;
                    }

                    break;
                }

                equal = true;

                for (i = 0; i < m; i++)
                {
                    if (edge[i + (j - 1) * m] != edge[i + j * m])
                    {
                        equal = false;
                    }
                }

                if (equal)
                {
                    j = j + 1;
                }
                else
                {
                    for (i = 0; i < m; i++)
                    {
                        node_boundary[edge[i + (j - 1) * m] - 1] = true;
                    }
                }

            }
            
            return node_boundary;
        }

    }
}