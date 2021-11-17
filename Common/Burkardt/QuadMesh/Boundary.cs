using System;
using Burkardt.Types;

namespace Burkardt.QuadMesh;

public static class Boundary
{
    public static int boundary_edge_count_q4_mesh(int element_num, int[] element_node)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BOUNDARY_EDGE_COUNT_Q4_MESH counts the boundary edges.
        //
        //  Discussion:
        //
        //    This routine is given a Q4 mesh, an abstract list of sets of 4 nodes.
        //    It is assumed that the nodes in each Q4 are listed
        //    in a counterclockwise order, although the routine should work 
        //    if the nodes are consistently listed in a clockwise order as well.
        //
        //    It is assumed that each edge of the mesh is either 
        //    * an INTERIOR edge, which is listed twice, once with positive
        //      orientation and once with negative orientation, or;
        //    * a BOUNDARY edge, which will occur only once.
        //
        //    This routine should work even if the region has holes.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 February 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ELEMENT_NUM, the number of elements.
        //
        //    Input, int ELEMENT_NODE[4*ELEMENT_NUM], the nodes 
        //    that make up the elements.  These should be listed in counterclockwise 
        //    order.
        //
        //    Output, int BOUNDARY_EDGE_COUNT_Q4_MESH, the number of boundary 
        //    edges.
        //
    {
        int boundary_edge_num;
        int e1;
        int e2;
        int[] edge;
        int element;
        int interior_edge_num;
        int j;
        int m;
        int n;
        int unique_num;

        m = 2;
        n = 4 * element_num;
        //
        //  Set up the edge array.
        //
        edge = new int[2 * 4 * element_num];

        for (element = 0; element < element_num; element++)
        {
            edge[0 + element * 2 + element_num * 2 * 0] = element_node[0 + element * 4];
            edge[1 + element * 2 + element_num * 2 * 0] = element_node[1 + element * 4];

            edge[0 + element * 2 + element_num * 2 * 1] = element_node[1 + element * 4];
            edge[1 + element * 2 + element_num * 2 * 1] = element_node[2 + element * 4];

            edge[0 + element * 2 + element_num * 2 * 2] = element_node[2 + element * 4];
            edge[1 + element * 2 + element_num * 2 * 2] = element_node[3 + element * 4];

            edge[0 + element * 2 + element_num * 2 * 3] = element_node[3 + element * 4];
            edge[1 + element * 2 + element_num * 2 * 3] = element_node[0 + element * 4];
        }

        //
        //  In each column, force the smaller entry to appear first.
        //
        for (j = 0; j < n; j++)
        {
            e1 = Math.Min(edge[0 + 2 * j], edge[1 + 2 * j]);
            e2 = Math.Max(edge[0 + 2 * j], edge[1 + 2 * j]);
            edge[0 + 2 * j] = e1;
            edge[1 + 2 * j] = e2;
        }

        //
        //  Ascending sort the column array.
        //
        typeMethods.i4col_sort_a(m, n, ref edge);
        //
        //  Get the number of unique columns in EDGE.
        //
        unique_num = typeMethods.i4col_sorted_unique_count(m, n, edge);

        interior_edge_num = 4 * element_num - unique_num;

        boundary_edge_num = 4 * element_num - 2 * interior_edge_num;


        return boundary_edge_num;
    }

    public static int boundary_edge_count_euler_q4_mesh(int node_num, int element_num,
            int hole_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BOUNDARY_EDGE_COUNT_EULER_Q4_MESH counts boundary edges.
        //
        //  Discussion:
        //
        //    We assume we are given information about a quadrilateral mesh
        //    of a set of nodes in the plane.
        //
        //    Given the number of nodes, elements and holes, we are going to apply
        //    Euler's formula to determine the number of edges that lie on the
        //    boundary of the set of nodes.
        //
        //    The number of faces, including the infinite face and internal holes, 
        //    is ELEMENT_NUM + HOLE_NUM + 1.
        //
        //    Let BOUNDARY_NUM denote the number of edges on the boundary.
        //    Each of the ELEMENT_NUM quadrilaterals uses four edges.  Every edge
        //    occurs in two different elements, so the number of edges must be
        //    ( 4 * ELEMENT_NUM + BOUNDARY_NUM ) / 2.
        //
        //    The number of nodes used in the mesh is NODE_NUM.
        //
        //    Euler's formula asserts that, for a simple connected figure in the
        //    plane with no edge crossings, NODE_NUM nodes, EDGE_NUM edges and
        //    FACE_NUM faces:
        //
        //      NODE_NUM - EDGE_NUM + FACE_NUM = 2
        //
        //    In our context, this becomes
        //
        //      NODE_NUM - ( 4 * ELEMENT_NUM + BOUNDARY_NUM ) / 2 
        //      + ELEMENT_NUM + HOLE_NUM + 1 = 2
        //
        //    or
        //
        //      BOUNDARY_NUM = 2 * NODE_NUM + 2 * HOLE_NUM - 2 * ELEMENT_NUM - 2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 February 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Marc de Berg, Marc Krevald, Mark Overmars, Otfried Schwarzkopf,
        //    Computational Geometry, Section 9.1,
        //    Springer, 2000.
        //
        //  Parameters:
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, int ELEMENT_NUM, the number of elements.
        //
        //    Input, int HOLE_NUM, the number of internal holes.
        //
        //    Output, int BOUNDARY_EDGE_COUNT_EULER_Q4_MESH, the number of edges that 
        //    lie on the boundary of the mesh.
        //
    {
        int boundary_num;

        boundary_num = 2 * node_num + 2 * hole_num - 2 * element_num - 2;

        return boundary_num;
    }
}