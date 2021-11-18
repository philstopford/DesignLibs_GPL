using System;
using Burkardt.Types;

namespace Burkardt.TriangulationNS;

public static partial class Boundary
{
    public static int triangulation_order3_boundary_edge_count(int triangle_num,
            int[] triangle_node)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULATION_ORDER3_BOUNDARY_EDGE_COUNT counts the boundary edges.
        //
        //  Discussion:
        //
        //    This routine is given a triangulation, an abstract list of triples
        //    of nodes.  It is assumed that the nodes in each triangle are listed
        //    in a counterclockwise order, although the routine should work
        //    if the nodes are consistently listed in a clockwise order as well.
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
        //    12 June 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int TRIANGLE_NUM, the number of triangles.
        //
        //    Input, int TRIANGLE_NODE[3*TRIANGLE_NUM], the nodes that make up the
        //    triangles.  These should be listed in counterclockwise order.
        //
        //    Output, integer TRIANGULATION_ORDER3_BOUNDARY_EDGE_COUNT, the number
        //    of boundary edges.
        //
    {
        int j;

        const int m = 2;
        int n = 3 * triangle_num;
        //
        //  Set up the edge array.
        //
        int[] edge = new int[m * n];

        for (j = 0; j < triangle_num; j++)
        {
            edge[0 + j * m] = triangle_node[0 + j * 3];
            edge[1 + j * m] = triangle_node[1 + j * 3];
            edge[0 + (j +     triangle_num) * m] = triangle_node[1 + j * 3];
            edge[1 + (j +     triangle_num) * m] = triangle_node[2 + j * 3];
            edge[0 + (j + 2 * triangle_num) * m] = triangle_node[2 + j * 3];
            edge[1 + (j + 2 * triangle_num) * m] = triangle_node[0 + j * 3];
        }

        //
        //  In each column, force the smaller entry to appear first.
        //
        for (j = 0; j < n; j++)
        {
            int e1 = Math.Min(edge[0 + j * m], edge[1 + j * m]);
            int e2 = Math.Max(edge[0 + j * m], edge[1 + j * m]);
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
        int unique_num = typeMethods.i4col_sorted_unique_count(m, n, edge);

        int interior_edge_num = 3 * triangle_num - unique_num;

        int boundary_edge_num = 3 * triangle_num - 2 * interior_edge_num;

        return boundary_edge_num;
    }

    public static int triangulation_order3_boundary_edge_count_euler(int node_num,
            int triangle_num, int hole_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULATION_ORDER3_BOUNDARY_EDGE_COUNT_EULER counts boundary edges.
        //
        //  Discussion:
        //
        //    We assume we are given information about a triangulation
        //    of a set of nodes in the plane.
        //
        //    Given the number of nodes and triangles, we are going to apply
        //    Euler's formula to determine the number of edges that lie on the
        //    boundary of the set of nodes.
        //
        //    The number of faces, including the infinite face and internal holes,
        //    is TRIANGLE_NUM + HOLE_NUM + 1.
        //
        //    Let BOUNDARY_NUM denote the number of edges on the boundary.
        //    Each of the TRIANGLE_NUM triangles uses three edges.  Every edge
        //    occurs in two different faces, so the number of edges must be
        //    ( 3 * TRIANGLE_NUM + BOUNDARY_NUM ) / 2.
        //
        //    The number of nodes used in the triangulation is NODE_NUM.
        //
        //    Euler's formula asserts that, for a simple connected figure in the
        //    plane with no edge crossings, NODE_NUM nodes, EDGE_NUM edges and
        //    FACE_NUM faces:
        //
        //      NODE_NUM - EDGE_NUM + FACE_NUM = 2
        //
        //    In our context, this becomes
        //
        //      NODE_NUM - ( 3 * TRIANGLE_NUM + BOUNDARY_NUM ) / 2
        //      + TRIANGLE_NUM + HOLE_NUM + 1 = 2
        //
        //    or
        //
        //      BOUNDARY_NUM = 2 * NODE_NUM + 2 * HOLE_NUM - TRIANGLE_NUM - 2
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
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, int TRIANGLE_NUM, the number of triangles.
        //
        //    Input, int HOLE_NUM, the number of holes.
        //
        //    Output, int TRIANGULATION_BOUNDARY_COUNT, the number of edges that
        //    lie on the convex hull of the triangulation.
        //
    {
        return 2 * node_num + 2 * hole_num - triangle_num - 2;
    }

    public static bool[] triangulation_order3_boundary_node(int node_num, int triangle_num,
            int[] triangle_node)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULATION_ORDER3_BOUNDARY_NODE indicates nodes on the boundary.
        //
        //  Discussion:
        //
        //    This routine is given a triangulation, an abstract list of triples
        //    of nodes.  It is assumed that the nodes in each triangle are listed
        //    in a counterclockwise order, although the routine should work
        //    if the nodes are consistently listed in a clockwise order as well.
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
        //    Input, int TRIANGLE_NODE[3*TRIANGLE_NUM], the nodes that make up the
        //    triangles.  These should be listed in counterclockwise order.
        //
        //    Output, bool TRIANGULATION_ORDER3_BOUNDARY_NODE[NODE_NUM],
        //    is TRUE if the node is on a boundary edge.
        //
    {
        int i;
        int j;

        const int m = 2;
        int n = 3 * triangle_num;
        //
        //  Set up the edge array.
        //
        int[] edge = new int[m * n];

        for (j = 0; j < triangle_num; j++)
        {
            edge[0 + j * m] = triangle_node[0 + j * 3];
            edge[1 + j * m] = triangle_node[1 + j * 3];
            edge[0 + (j + triangle_num) * m] = triangle_node[1 + j * 3];
            edge[1 + (j + triangle_num) * m] = triangle_node[2 + j * 3];
            edge[0 + (j + 2 * triangle_num) * m] = triangle_node[2 + j * 3];
            edge[1 + (j + 2 * triangle_num) * m] = triangle_node[0 + j * 3];
        }

        //
        //  In each column, force the smaller entry to appear first.
        //
        for (j = 0; j < n; j++)
        {
            int e1 = Math.Min(edge[0 + j * m], edge[1 + j * m]);
            int e2 = Math.Max(edge[0 + j * m], edge[1 + j * m]);
            edge[0 + j * m] = e1;
            edge[1 + j * m] = e2;
        }

        //
        //  Ascending sort the column array.
        //
        typeMethods.i4col_sort_a(m, n, ref edge);
        //
        //  Records which appear twice are internal edges and can be ignored.
        //
        bool[] node_boundary = new bool[node_num];

        for (i = 0; i < node_num; i++)
        {
            node_boundary[i] = false;
        }

        j = 0;

        while (j < 3 * triangle_num)
        {
            j += 1;

            if (j == 3 * triangle_num)
            {
                for (i = 0; i < m; i++)
                {
                    node_boundary[edge[i + (j - 1) * m] - 1] = true;
                }

                break;
            }

            bool equal = true;

            for (i = 0; i < m; i++)
            {
                if (edge[i + (j - 1) * m] != edge[i + j * m])
                {
                    equal = false;
                }
            }

            switch (equal)
            {
                case true:
                    j += 1;
                    break;
                default:
                {
                    for (i = 0; i < m; i++)
                    {
                        node_boundary[(node_boundary.Length + (edge[(edge.Length + i + (j - 1) * m) % edge.Length] - 1)) % node_boundary.Length] = true;
                    }

                    break;
                }
            }

        }
            
        return node_boundary;
    }
}