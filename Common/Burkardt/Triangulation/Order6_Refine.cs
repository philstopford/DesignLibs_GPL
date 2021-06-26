using System;
using Burkardt.Types;

namespace Burkardt.TriangulationNS
{
    public static partial class Refine
    {
        public static void triangulation_order6_refine_compute(int node_num1, int triangle_num1,
            double[] node_xy1, int[] triangle_node1, int node_num2, int triangle_num2,
        int[] edge_data, ref double[] node_xy2, ref int[] triangle_node2 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULATION_ORDER6_REFINE_COMPUTE computes a refined order 6 triangulation.
        //
        //  Discussion:
        //
        //    Given a quadratic triangle defined by nodes 1, 2, 3, 4, 5, 6, we
        //    need to generate nodes 14, 16, 24, 25, 35, 36, 45, 46, 56, and 4 new
        //    quadratic subtriangles T1, T2, T3 and T4.
        //
        //    The task is more complicated by the fact that we are working with
        //    a mesh of triangles, so that we want to create a node only once,
        //    even though it may be shared by other triangles.  (In fact, only
        //    the new nodes on the edges can be shared, and then only by at most
        //    one other triangle.)
        //
        //            3
        //           . .
        //          36 35
        //         . T3  .
        //        6--56---5
        //       . . T4  . .
        //      16 46  45  25
        //     . T1  . . T2  .
        //    1--14---4--24---2
        //
        //    This routine is given sorted information defining the edges, and uses
        //    it to build the new node and triangle arrays.
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
        //    Input, int TRIANGLE_NODE1[6*TRIANGLE_NUM1], the nodes that make up the
        //    triangles.
        //
        //    Input, int NODE_NUM2, the number of nodes in the refined mesh.
        //
        //    Input, int TRIANGLE_NUM2, the number of triangles in the refined mesh.
        //
        //    Input, int EDGE_DATA[5*(3*TRIANGLE_NUM1)], edge information computed
        //    by TRIANGULATION_ORDER6_REFINE_SIZE.
        //
        //    Output, double NODE_XY2[2*NODE_NUM2], the refined nodes.
        //
        //    Output, int TRIANGLE_NODE2[6*TRIANGLE_NUM2], the nodes that make up the
        //    triangles in the refined mesh.
        //
        {
            int edge;
            int i;
            int j;
            int l1 = 0;
            int l2 = 0;
            int l3 = 0;
            int n1;
            int n1_old;
            int n2;
            int n2_old;
            int node;
            int t1;
            int t2;
            int t3;
            int t4;
            int triangle1;
            int v1 = 0;
            int v2 = 0;
            int v3 = 0;
            int v4 = 0;
            int v5 = 0;
            int v6 = 0;
            //
            //  Step 1:
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
                for (i = 0; i < 6; i++)
                {
                    triangle_node2[i + j * 6] = -1;
                }
            }

            //
            //  We can assign the existing nodes to the new triangles.
            //
            for (triangle1 = 0; triangle1 < triangle_num1; triangle1++)
            {
                t1 = triangle1 * 4 + 0;
                t2 = triangle1 * 4 + 1;
                t3 = triangle1 * 4 + 2;
                t4 = triangle1 * 4 + 3;

                triangle_node2[0 + t1 * 6] = triangle_node1[0 + triangle1 * 6];
                triangle_node2[1 + t1 * 6] = triangle_node1[3 + triangle1 * 6];
                triangle_node2[2 + t1 * 6] = triangle_node1[5 + triangle1 * 6];

                triangle_node2[0 + t2 * 6] = triangle_node1[3 + triangle1 * 6];
                triangle_node2[1 + t2 * 6] = triangle_node1[1 + triangle1 * 6];
                triangle_node2[2 + t2 * 6] = triangle_node1[4 + triangle1 * 6];

                triangle_node2[0 + t3 * 6] = triangle_node1[5 + triangle1 * 6];
                triangle_node2[1 + t3 * 6] = triangle_node1[4 + triangle1 * 6];
                triangle_node2[2 + t3 * 6] = triangle_node1[2 + triangle1 * 6];

                triangle_node2[0 + t4 * 6] = triangle_node1[4 + triangle1 * 6];
                triangle_node2[1 + t4 * 6] = triangle_node1[5 + triangle1 * 6];
                triangle_node2[2 + t4 * 6] = triangle_node1[3 + triangle1 * 6];
            }

            //
            //  Step 2.
            //  Examine sorted edge information.  The first time an edge is encountered,
            //  generate two new nodes, then assign them (usually) to the four subtriangles
            //  of the two triangles that share that edge.
            //
            node = node_num1;

            n1_old = -1;
            n2_old = -1;

            for (edge = 0; edge < 3 * triangle_num1; edge++)
            {
                n1 = edge_data[0 + edge * 5] - 1;
                n2 = edge_data[1 + edge * 5] - 1;

                l1 = edge_data[2 + edge * 5];
                l3 = edge_data[3 + edge * 5];

                if (l1 == 1 && l3 == 2)
                {
                    l2 = 4;
                }
                else if (l1 == 1 && l3 == 3)
                {
                    l2 = 6;
                }
                else if (l1 == 2 && l3 == 3)
                {
                    l2 = 5;
                }

                triangle1 = edge_data[4 + edge * 5];
                //
                //  If this is the first time we've encountered this edge,
                //  create the new nodes.
                //
                if (n1 != n1_old || n2 != n2_old)
                {
                    n1_old = n1;
                    n2_old = n2;

                    v1 = triangle_node1[l1 - 1 + triangle1 * 6];
                    v2 = triangle_node1[l2 - 1 + triangle1 * 6];
                    v3 = triangle_node1[l3 - 1 + triangle1 * 6];

                    for (i = 0; i < 2; i++)
                    {
                        node_xy2[i + node * 2] = (node_xy2[i + (v1 - 1) * 2]
                                                  + node_xy2[i + (v2 - 1) * 2]) / 2.0;
                    }

                    node = node + 1;
                    v4 = node;

                    for (i = 0; i < 2; i++)
                    {
                        node_xy2[i + node * 2] = (node_xy2[i + (v2 - 1) * 2]
                                                  + node_xy2[i + (v3 - 1) * 2]) / 2.0;
                    }

                    node = node + 1;
                    v5 = node;
                }

                t1 = triangle1 * 4 + 0;
                t2 = triangle1 * 4 + 1;
                t3 = triangle1 * 4 + 2;

                if (l1 == 1 && l3 == 2)
                {
                    if (triangle_node1[0 + triangle1 * 6] == v1 + 1)
                    {
                        triangle_node2[3 + t1 * 6] = v4;
                        triangle_node2[3 + t2 * 6] = v5;
                    }
                    else
                    {
                        triangle_node2[3 + t1 * 6] = v5;
                        triangle_node2[3 + t2 * 6] = v4;
                    }
                }
                else if (l1 == 1 && l3 == 3)
                {
                    if (triangle_node1[0 + triangle1 * 6] == v1 + 1)
                    {
                        triangle_node2[5 + t1 * 6] = v4;
                        triangle_node2[5 + t3 * 6] = v5;
                    }
                    else
                    {
                        triangle_node2[5 + t1 * 6] = v5;
                        triangle_node2[5 + t3 * 6] = v4;
                    }
                }
                else if (l1 == 2 && l3 == 3)
                {
                    if (triangle_node1[1 + triangle1 * 6] == v1 + 1)
                    {
                        triangle_node2[4 + t3 * 6] = v4;
                        triangle_node2[4 + t2 * 6] = v5;
                    }
                    else
                    {
                        triangle_node2[4 + t3 * 6] = v5;
                        triangle_node2[4 + t2 * 6] = v4;
                    }
                }
            }

            //
            //  Step 3.
            //  Each old triangle has a single central subtriangle, for which we now
            //  need to generate three new "interior" nodes.
            //
            for (triangle1 = 0; triangle1 < triangle_num1; triangle1++)
            {
                v4 = triangle_node1[3 + triangle1 * 6];
                v5 = triangle_node1[4 + triangle1 * 6];
                v6 = triangle_node1[5 + triangle1 * 6];

                t1 = triangle1 * 4 + 0;
                t2 = triangle1 * 4 + 1;
                t3 = triangle1 * 4 + 2;
                t4 = triangle1 * 4 + 3;

                node_xy2[0 + node * 2] = 0.5 * (node_xy1[0 + (v5 - 1) * 2] + node_xy1[0 + (v6 - 1) * 2]);
                node_xy2[1 + node * 2] = 0.5 * (node_xy1[1 + (v5 - 1) * 2] + node_xy1[1 + (v6 - 1) * 2]);
                node = node + 1;
                triangle_node2[3 + t4 * 6] = node;
                triangle_node2[3 + t3 * 6] = node;

                node_xy2[0 + node * 2] = 0.5 * (node_xy1[0 + (v6 - 1) * 2] + node_xy1[0 + (v4 - 1) * 2]);
                node_xy2[1 + node * 2] = 0.5 * (node_xy1[1 + (v6 - 1) * 2] + node_xy1[1 + (v4 - 1) * 2]);
                node = node + 1;
                triangle_node2[4 + t4 * 6] = node;
                triangle_node2[4 + t1 * 6] = node;

                node_xy2[0 + node * 2] = 0.5 * (node_xy1[0 + (v4 - 1) * 2] + node_xy1[0 + (v5 - 1) * 2]);
                node_xy2[1 + node * 2] = 0.5 * (node_xy1[1 + (v4 - 1) * 2] + node_xy1[1 + (v5 - 1) * 2]);
                node = node + 1;
                triangle_node2[5 + t4 * 6] = node;
                triangle_node2[5 + t2 * 6] = node;
            }

            return;
        }

        public static void triangulation_order6_refine_size(int node_num1, int triangle_num1,
            int[] triangle_node1, ref int node_num2, ref int triangle_num2, ref int[] edge_data )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULATION_ORDER6_REFINE_SIZE sizes a refined order 6 triangulation.
        //
        //  Discussion:
        //
        //    Given a quadratic triangle defined by nodes 1, 2, 3, 4, 5, 6, we
        //    need to generate nodes 14, 16, 24, 25, 35, 36, 45, 46, 56, and 4 new
        //    quadratic subtriangles T1, T2, T3 and T4.
        //
        //    The task is more complicated by the fact that we are working with
        //    a mesh of triangles, so that we want to create a node only once,
        //    even though it may be shared by other triangles.  (In fact, only
        //    the new nodes on the edges can be shared, and then only by at most
        //    one other triangle.)
        //
        //            3
        //           . .
        //          36 35
        //         . T3  .
        //        6--56---5
        //       . . T4  . .
        //      16 46  45  25
        //     . T1  . . T2  .
        //    1--14---4--24---2
        //
        //    This routine determines the sizes of the resulting node and
        //    triangles, and constructs an edge array that can be used to
        //    properly number the new nodes.
        //
        //    The primary work occurs in sorting a list related to the edges.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 February 2007
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
        //    Input, int TRIANGLE_NODE1[6*TRIANGLE_NUM1], the indices of the nodes
        //    that form the triangles in the input mesh.
        //
        //    Output, int *NODE_NUM2, the number of nodes in the refined mesh.
        //
        //    Output, int *TRIANGLE_NUM2, the number of triangles in the
        //    refined mesh.
        //
        //    Output, int EDGE_DATA[5*(3*TRIANGLE_NUM1)], edge data that will
        //    be needed by TRIANGULATION_ORDER6_REFINE_COMPUTE.
        //
        {
            int a;
            int b;
            int edge;
            int i;
            int j;
            int k;
            int n1;
            int n1_old;
            int n2;
            int n2_old;
            int triangle1;
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
            for (triangle1 = 0; triangle1 < triangle_num1; triangle1++)
            {
                i = triangle_node1[0 + triangle1 * 6];
                j = triangle_node1[1 + triangle1 * 6];
                k = triangle_node1[2 + triangle1 * 6];

                a = Math.Min(i, j);
                b = Math.Max(i, j);

                edge_data[0 + 5 * (3 * triangle1 + 0)] = a;
                edge_data[1 + 5 * (3 * triangle1 + 0)] = b;
                edge_data[2 + 5 * (3 * triangle1 + 0)] = 1;
                edge_data[3 + 5 * (3 * triangle1 + 0)] = 2;
                edge_data[4 + 5 * (3 * triangle1 + 0)] = triangle1;

                a = Math.Min(i, k);
                b = Math.Max(i, k);

                edge_data[0 + 5 * (3 * triangle1 + 1)] = a;
                edge_data[1 + 5 * (3 * triangle1 + 1)] = b;
                edge_data[2 + 5 * (3 * triangle1 + 1)] = 1;
                edge_data[3 + 5 * (3 * triangle1 + 1)] = 3;
                edge_data[4 + 5 * (3 * triangle1 + 1)] = triangle1;

                a = Math.Min(j, k);
                b = Math.Max(j, k);

                edge_data[0 + 5 * (3 * triangle1 + 2)] = a;
                edge_data[1 + 5 * (3 * triangle1 + 2)] = b;
                edge_data[2 + 5 * (3 * triangle1 + 2)] = 2;
                edge_data[3 + 5 * (3 * triangle1 + 2)] = 3;
                edge_data[4 + 5 * (3 * triangle1 + 2)] = triangle1;
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

            n1_old = -1;
            n2_old = -1;

            for (edge = 0; edge < 3 * triangle_num1; edge++)
            {
                n1 = edge_data[0 + edge * 5];
                n2 = edge_data[1 + edge * 5];
                if (n1 != n1_old || n2 != n2_old)
                {
                    node_num2 = node_num2 + 2;
                    n1_old = n1;
                    n2_old = n2;
                }
            }

            node_num2 = node_num2 + 3 * triangle_num1;

            triangle_num2 = 4 * triangle_num1;
        }
    }
}