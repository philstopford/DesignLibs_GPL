using System;
using Burkardt.Types;

namespace Burkardt.TriangulationNS
{
    public static partial class Neighbor
    {
        public static void triangulation_order3_neighbor(int triangle_num, ref int[] triangle_node,
                int t1, int s1, ref int t2, ref int s2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGULATION_ORDER3_NEIGHBOR determines a neighbor of a given triangle.
            //
            //  Discussion:
            //
            //    A set of nodes is given.  A triangulation of the nodes has been
            //    defined and recorded in TRIANGLE_NODE.  The TRIANGLE_NODE data structure
            //    records triangles as sets of three nodes, N1, N2, N3, that implicitly
            //    define three sides, being the line segments N1-N2, N2-N3, and N3-N1.
            //
            //    The nodes of the triangle are listed in counterclockwise order.
            //    This means that if two triangles share a side, then the nodes
            //    defining that side occur in the order (N1,N2) for one triangle,
            //    and (N2,N1) for the other.
            //
            //    The routine is given a triangle and a side, and asked to find
            //    another triangle (if any) that shares that side.  The routine
            //    simply searches the TRIANGLE_NODE structure for an occurrence of the
            //    nodes in the opposite order.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 October 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int TRIANGLE_NUM, the number of triangles.
            //
            //    Input/output, int TRIANGLE_NODE[3*TRIANGLE_NUM], the nodes that define
            //    each triangle.
            //
            //    Input, int T1, the index of the triangle.
            //
            //    Input, int S1, the index of the triangle side.
            //
            //    Output, int *T2, the index of the triangle which is the neighbor
            //    to T1 on side S1, or -1 if there is no such neighbor.
            //
            //    Output, int *S2, the index of the side of triangle T2 which
            //    is shared with triangle T1, or -1 if there is no such neighbor.
            //
        {
            int n1;
            int n2;
            int s;
            int ss;
            int t;

            n1 = triangle_node[s1 - 1 + (t1 - 1) * 3];
            ss = typeMethods.i4_wrap(s1 + 1, 1, 3);
            n2 = triangle_node[ss - 1 + (t1 - 1) * 3];

            for (t = 0; t < triangle_num; t++)
            {
                for (s = 0; s < 3; s++)
                {
                    if (triangle_node[s + t * 3] == n1)
                    {
                        ss = typeMethods.i4_wrap(s - 1, 0, 2);
                        if (triangle_node[ss + t * 3] == n2)
                        {
                            t2 = t + 1;
                            s2 = ss + 1;
                            return;
                        }
                    }
                }
            }

            t2 = -1;
            s2 = -1;

            return;
        }

        public static void triangulation_order3_neighbor_nodes(int node_num, int triangle_num,
                int[] triangle_node, ref int[] nabes_first, ref int[] nabes_num, int nabes_max,
                ref int nabes_dim, ref int[] nabes)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGULATION_ORDER3_NEIGHBOR_NODES determines node neighbors.
            //
            //  Example:
            //
            //    On input, the triangle data structure is:
            //
            //    Triangle  Nodes
            //    --------  ----------
            //     1        3,   4,   1
            //     2        3,   1,   2
            //     3        3,   2,   6
            //     4        2,   1,   5
            //     5        6,   2,   5
            //
            //  On output, the auxilliary neighbor arrays are:
            //
            //    Node  Num  First
            //    ----  ---  -----
            //     1     4     1
            //     2     4     5
            //     3     4     9
            //     4     2    13
            //     5     3    15
            //     6     3    18
            //
            //  and the neighbor array is:
            //
            //    Position  Node
            //    --------  ----
            //
            //     1        2
            //     2        3
            //     3        4
            //     4        5
            //    -----------
            //     5        1
            //     6        3
            //     7        5
            //     8        6
            //    -----------
            //     9        1
            //    10        2
            //    11        4
            //    12        6
            //    -----------
            //    13        1
            //    14        3
            //    -----------
            //    15        1
            //    16        2
            //    17        6
            //    -----------
            //    18        2
            //    19        3
            //    20        5
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    10 November 2003
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
            //    Input, int TRIANGLE_NODE[3*TRIANGLE_NUM], the nodes that make up
            //    each triangle.
            //
            //    Output, int NABES_FIRST[NODE_NUM], the index in NABES of the first
            //    neighbor in the list for each node.
            //
            //    Output, int NABES_NUM[NODE_NUM], the number of neighbors of each node.
            //
            //    Input, int NABES_MAX, the maximum dimension of NABES.
            //
            //    Output, int *NABES_DIM, the dimension of NABES.
            //
            //    Output, int NABES[*NABES_DIM], a list of the neighbors of all the nodes.
            //    Neighbors of node 1 are listed first, and so on.
            //
        {
            int i;
            int i_current;
            int j;
            int k;
            int n;
            int nabe;
            int[] nabes1;
            int tri;

            nabes = new int[nabes_max];
            nabes1 = new int[nabes_max];
            //
            //  Step 1.  From the triangle list (I,J,K)
            //  construct the neighbor relations: (I,J), (J,K), (K,I), (J,I), (K,J), (I,K).
            //
            n = 0;

            for (tri = 0; tri < triangle_num; tri++)
            {
                i = triangle_node[0 + tri * 3];
                j = triangle_node[1 + tri * 3];
                k = triangle_node[2 + tri * 3];
                nabes1[n] = i;
                nabes1[n + 1] = i;
                nabes1[n + 2] = j;
                nabes1[n + 3] = j;
                nabes1[n + 4] = k;
                nabes1[n + 5] = k;
                nabes[n] = j;
                nabes[n + 1] = k;
                nabes[n + 2] = i;
                nabes[n + 3] = k;
                nabes[n + 4] = i;
                nabes[n + 5] = j;

                n = n + 6;
            }

            //
            //  Step 2. Dictionary sort the neighbor relations.
            //
            typeMethods.i4vec2_sort_a(n, ref nabes1, ref nabes);
            //
            //  Step 3. Remove duplicate entries.
            //
            int nu = 0;
            typeMethods.i4vec2_sorted_unique(n, ref nabes1, ref nabes, ref nu);
            n = nu;
            //
            //  Step 4. Construct the NABES_NUM and NABES_FIRST data.
            //
            for (i = 0; i < node_num; i++)
            {
                nabes_num[i] = 0;
            }

            for (i = 0; i < node_num; i++)
            {
                nabes_first[i] = 0;
            }

            i_current = 0;

            for (nabe = 1; nabe <= n; nabe++)
            {
                i = nabes1[nabe - 1];
                if (i == i_current)
                {
                    nabes_num[i - 1] = nabes_num[i - 1] + 1;
                }
                else
                {
                    i_current = i;
                    nabes_first[i - 1] = nabe;
                    nabes_num[i - 1] = 1;
                }
            }

            nabes_dim = n;
        }

        public static void triangulation_order3_neighbor_nodes_print(int node_num,
                int[] nabes_first, int[] nabes_num, int nabes_dim, int[] nabes)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGULATION_ORDER3_NEIGHBOR_NODES_PRINT prints a node neighbor array.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 July 2001
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NODE_NUM, the number of nodes.
            //
            //    Input, int NABES_FIRST[NODE_NUM], the index in NABES of the first
            //    neighbor in the list for each node.
            //
            //    Input, int NABES_NUM[NODE_NUM], the number of neighbors of each node.
            //
            //    Input, int NABES_DIM, the dimension of NABES.
            //
            //    Input, int NABES[NABES_DIM], a list of the neighbors of all the nodes.
            //    Neighbors of node 1 are listed first, and so on.
            //
        {
            int i;
            int j;
            int k;

            Console.WriteLine("");
            Console.WriteLine("  Node Nabes Index  List");
            Console.WriteLine("");
            string cout = "";

            for (i = 0; i < node_num; i++)
            {
                cout = i.ToString().PadLeft(4) + "  "
                                  + nabes_num[i].ToString().PadLeft(4) + "  "
                                  + nabes_first[i].ToString().PadLeft(4) + "  ";

                k = 0;
                for (j = nabes_first[i] - 1; j < nabes_first[i] + nabes_num[i]; j++)
                {
                    if (k == 10)
                    {
                        Console.WriteLine(cout);
                        cout = "                  ";
                        k = 0;
                    }

                    cout += nabes[j].ToString().PadLeft(4) + "  ";
                    k = k + 1;
                }
            }

            Console.WriteLine(cout);
        }
    }
}