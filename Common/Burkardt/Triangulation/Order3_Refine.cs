using System;
using Burkardt.Types;

namespace Burkardt.TriangulationNS;

public static partial class Refine
{
    public static void triangulation_order3_refine_compute(int node_num1, int triangle_num1,
            double[] node_xy1, int[] triangle_node1, int node_num2, int triangle_num2,
            int[] edge_data, ref double[] node_xy2, ref int[] triangle_node2 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULATION_ORDER3_REFINE_COMPUTE computes a refined order 3 triangulation.
        //
        //  Discussion:
        //
        //    Given a triangle defined by nodes 1, 2, 3, we need to generate
        //    nodes 12, 23, and 13, and create 4 new subtriangles, T1, T2, T3
        //    and T4.
        //
        //    The task is more complicated by the fact that we are working with
        //    a mesh of triangles, so that we want to create a node only once,
        //    even though it may be shared by other triangles.
        //
        //          3
        //         . .
        //        .T3 .
        //      13----23
        //      . .T4 . .
        //     .T1 . .T2 .
        //    1----12-----2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NODE_NUM1, the number of nodes.
        //
        //    Input, int TRIANGLE_NUM1, the number of triangles.
        //
        //    Input, double NODE_XY1[2*NODE_NUM1], the nodes.
        //
        //    Input, int TRIANGLE_NODE1[3*TRIANGLE_NUM1], the nodes that make up the
        //    triangles.  These should be listed in counterclockwise order.
        //
        //    Input, int NODE_NUM2, the number of nodes in the refined mesh.
        //
        //    Input, int TRIANGLE_NUM2, the number of triangles in the refined mesh.
        //
        //    Input, int EDGE_DATA[5*(3*TRIANGLE_NUM1)], edge information computed
        //    by TRIANGULATION_ORDER3_REFINE_SIZE.
        //
        //    Output, double NODE_XY2[2*NODE_NUM2], the refined nodes.
        //
        //    Output, int TRIANGLE_NODE2[3*TRIANGLE_NUM2], the nodes that make up the
        //    triangles in the refined mesh.
        //
    {
        int edge;
        int i;
        int j;
        int triangle1;
        //
        //  Copy the old nodes.
        //
        for (j = 0; j < node_num1; j++)
        {
            for (i = 0; i < 2; i++)
            {
                node_xy2[i + j * 2] = node_xy1[i + j * 2];
            }
        }

        for (j = 0; j < triangle_num2; j++)
        {
            for (i = 0; i < 3; i++)
            {
                triangle_node2[i + j * 3] = -1;
            }
        }

        //
        //  We can assign the existing nodes to the new triangles.
        //
        for (triangle1 = 0; triangle1 < triangle_num1; triangle1++)
        {
            triangle_node2[0 + (triangle1 * 4 + 0) * 3] = triangle_node1[0 + triangle1 * 3];
            triangle_node2[1 + (triangle1 * 4 + 1) * 3] = triangle_node1[1 + triangle1 * 3];
            triangle_node2[2 + (triangle1 * 4 + 2) * 3] = triangle_node1[2 + triangle1 * 3];
        }

        int node = node_num1;

        int n1_old = -1;
        int n2_old = -1;

        for (edge = 0; edge < 3 * triangle_num1; edge++)
        {
            int n1 = edge_data[0 + edge * 5] - 1;
            int n2 = edge_data[1 + edge * 5] - 1;
            //
            //  If this edge is new, create the coordinates and index for this node.
            //
            if (n1 != n1_old || n2 != n2_old)
            {

                if (node_num2 < node)
                {
                    Console.WriteLine("");
                    Console.WriteLine("TRIANGLE_MESH_ORDER3_REFINE - Fatal error!");
                    Console.WriteLine("  Node index exceeds NODE_NUM2.");
                    return;
                }

                for (i = 0; i < 2; i++)
                {
                    node_xy2[(i + node * 2 + node_xy2.Length) % node_xy2.Length] = (node_xy2[(i + n1 * 2 + node_xy2.Length) % node_xy2.Length] + node_xy2[(i + n2 * 2 + node_xy2.Length) % node_xy2.Length]) / 2.0;
                }

                node += 1;

                n1_old = n1;
                n2_old = n2;
            }

            //
            //  Assign the node to triangles.
            //
            int v1 = edge_data[2 + edge * 5];
            int v2 = edge_data[3 + edge * 5];
            triangle1 = edge_data[4 + edge * 5];

            switch (v1)
            {
                case 1 when v2 == 2:
                    triangle_node2[0 + (triangle1 * 4 + 1) * 3] = node;
                    triangle_node2[1 + (triangle1 * 4 + 0) * 3] = node;
                    triangle_node2[2 + (triangle1 * 4 + 3) * 3] = node;
                    break;
                case 1 when v2 == 3:
                    triangle_node2[0 + (triangle1 * 4 + 2) * 3] = node;
                    triangle_node2[1 + (triangle1 * 4 + 3) * 3] = node;
                    triangle_node2[2 + (triangle1 * 4 + 0) * 3] = node;
                    break;
                case 2 when v2 == 3:
                    triangle_node2[0 + (triangle1 * 4 + 3) * 3] = node;
                    triangle_node2[1 + (triangle1 * 4 + 2) * 3] = node;
                    triangle_node2[2 + (triangle1 * 4 + 1) * 3] = node;
                    break;
            }
        }
    }

    public static void triangulation_order3_refine_size(int node_num1, int triangle_num1,
            int[] triangle_node1, ref int node_num2, ref int triangle_num2, ref int[] edge_data )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULATION_ORDER3_REFINE_SIZE sizes a refined order 3 triangulation.
        //
        //  Discussion:
        //
        //    Given a triangle defined by nodes 1, 2, 3, we need to generate
        //    nodes 12, 23, and 13, and create 4 new subtriangles, T1, T2, T3
        //    and T4.
        //
        //    The task is more complicated by the fact that we are working with
        //    a mesh of triangles, so that we want to create a node only once,
        //    even though it may be shared by other triangles.
        //
        //          3
        //         . .
        //        .T3 .
        //      13----23
        //      . .T4 . .
        //     .T1 . .T2 .
        //    1----12-----2
        //
        //    This routine simply determines the sizes of the resulting node
        //    and triangle arrays.
        //
        //    The primary amount of work occurs in sorting a list of 3 * TRIANGLE_NUM
        //    data items, one item for every edge of every triangle.  Each
        //    data item records, for a given edge, the global indices
        //    of the two endpoints, the local indices of the two endpoints,
        //    and the index of the triangle.
        //
        //    Through careful sorting, it is possible to arrange this data in
        //    a way that allows the proper generation of the interpolated nodes.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NODE_NUM1, the number of nodes in the original mesh.
        //
        //    Input, int  TRIANGLE_NUM1, the number of triangles in the
        //    original mesh.
        //
        //    Input, int TRIANGLE_NODE1[3*TRIANGLE_NUM1], the indices of the nodes
        //    that form the triangles in the input mesh.
        //
        //    Output, int *NODE_NUM2, the number of nodes in the refined mesh.
        //
        //    Output, int *TRIANGLE_NUM2, the number of triangles in the
        //    refined mesh.
        //
        //    Output, int EDGE_DATA[5*(3*TRIANGLE_NUM1)], edge data that will
        //    be needed by TRIANGULATION_ORDER3_REFINE_COMPUTE.
        //
    {
        int edge;
        int triangle;
        //
        //  Step 1.
        //  From the list of nodes for triangle T, of the form: (I,J,K)
        //  construct the edge relations:
        //
        //    (I,J,1,2,T)
        //    (I,K,1,3,T)
        //    (J,K,2,3,T)
        //
        //  In order to make matching easier, we reorder each pair of nodes
        //  into ascending order.
        //
        for (triangle = 0; triangle < triangle_num1; triangle++)
        {
            int i = triangle_node1[0 + triangle * 3];
            int j = triangle_node1[1 + triangle * 3];
            int k = triangle_node1[2 + triangle * 3];

            int a = Math.Min(i, j);
            int b = Math.Max(i, j);

            edge_data[0 + 5 * (3 * triangle + 0)] = a;
            edge_data[1 + 5 * (3 * triangle + 0)] = b;
            edge_data[2 + 5 * (3 * triangle + 0)] = 1;
            edge_data[3 + 5 * (3 * triangle + 0)] = 2;
            edge_data[4 + 5 * (3 * triangle + 0)] = triangle;

            a = Math.Min(i, k);
            b = Math.Max(i, k);

            edge_data[0 + 5 * (3 * triangle + 1)] = a;
            edge_data[1 + 5 * (3 * triangle + 1)] = b;
            edge_data[2 + 5 * (3 * triangle + 1)] = 1;
            edge_data[3 + 5 * (3 * triangle + 1)] = 3;
            edge_data[4 + 5 * (3 * triangle + 1)] = triangle;

            a = Math.Min(j, k);
            b = Math.Max(j, k);

            edge_data[0 + 5 * (3 * triangle + 2)] = a;
            edge_data[1 + 5 * (3 * triangle + 2)] = b;
            edge_data[2 + 5 * (3 * triangle + 2)] = 2;
            edge_data[3 + 5 * (3 * triangle + 2)] = 3;
            edge_data[4 + 5 * (3 * triangle + 2)] = triangle;
        }

        //
        //  Step 2. Perform an ascending dictionary sort on the neighbor relations.
        //  We only intend to sort on rows 1:2; the routine we call here
        //  sorts on the full column but that won't hurt us.
        //
        //  What we need is to find all cases where triangles share an edge.
        //  By sorting the columns of the EDGE_DATA array, we will put shared edges
        //  next to each other.
        //
        typeMethods.i4col_sort_a(5, 3 * triangle_num1, ref edge_data);
        //
        //  Step 3. All the triangles which share an edge show up as consecutive
        //  columns with identical first two entries.  Figure out how many new
        //  nodes there are, and allocate space for their coordinates.
        //
        node_num2 = node_num1;

        int n1_old = -1;
        int n2_old = -1;

        for (edge = 0; edge < 3 * triangle_num1; edge++)
        {
            int n1 = edge_data[0 + edge * 5];
            int n2 = edge_data[1 + edge * 5];
            if (n1 != n1_old || n2 != n2_old)
            {
                node_num2 += 1;
                n1_old = n1;
                n2_old = n2;
            }
        }

        triangle_num2 = 4 * triangle_num1;
    }
}