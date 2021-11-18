using Burkardt.Types;

namespace Burkardt.TriangulationNS;

public static partial class Adjacency
{
    public static int triangulation_order6_adj_count(int node_num, int triangle_num,
            int[] triangle_node, ref int[] triangle_neighbor, ref int[] adj_col )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULATION_ORDER6_ADJ_COUNT counts adjacencies in a triangulation.
        //
        //  Discussion:
        //
        //    This routine is called to count the adjacencies, so that the
        //    appropriate amount of memory can be set aside for storage when
        //    the adjacency structure is created.
        //
        //    The triangulation is assumed to involve 6-node triangles.
        //
        //    Two nodes are "adjacent" if they are both nodes in some triangle.
        //    Also, a node is considered to be adjacent to itself.
        //
        //  Diagram:
        //
        //       3
        //    s  |.
        //    i  | .
        //    d  |  .
        //    e  6   5  side 2
        //       |    .
        //    3  |     .
        //       |      .
        //       1---4---2
        //
        //         side 1
        //
        //    The local node numbering
        //
        //
        //   21-22-23-24-25
        //    |.    |.    |
        //    | .   | .   |
        //   16 17 18 19 20
        //    |   . |   . |
        //    |    .|    .|
        //   11-12-13-14-15
        //    |.    |.    |
        //    | .   | .   |
        //    6  7  8  9 10
        //    |   . |   . |
        //    |    .|    .|
        //    1--2--3--4--5
        //
        //    A sample grid.
        //
        //
        //    Below, we have a chart that lists the nodes adjacent to each node, with
        //    an asterisk to indicate the adjacency of the node to itself
        //    (in some cases, you want to count this self adjacency and in some
        //    you don't).
        //
        //    N   Adjacencies
        //
        //    1:  *  2  3  6  7 11
        //    2:  1  *  3  6  7 11
        //    3:  1  2  *  4  5  6  7  8  9 11 12 13
        //    4:  3  *  5  8  9 13
        //    5:  3  4  *  8  9 10 13 14 15
        //    6:  1  2  3  *  7 11
        //    7:  1  2  3  6  *  8 11 12 13
        //    8:  3  4  5  7  *  9 11 12 13
        //    9:  3  4  5  8  * 10 13 14 15
        //   10:  5  9  * 13 14 15
        //   11:  1  2  3  6  7  8  * 12 13 16 17 21
        //   12:  3  7  8 11  * 13 16 17 21
        //   13:  3  4  5  7  8  9 10 11 12  * 14 15 16 17 18 19 21 22 23
        //   14:  5  9 10 13  * 15 18 19 23
        //   15:  5  9 10 13 14  * 18 19 20 23 24 25
        //   16: 11 12 13  * 17 21
        //   17: 11 12 13 16  * 18 21 22 23
        //   18: 13 14 15 17  * 19 21 22 23
        //   19: 13 14 15 18  * 20 23 24 25
        //   20: 15 19  * 23 24 25
        //   21: 11 12 13 16 17 18  * 22 23
        //   22: 13 17 18 21  * 23
        //   23: 13 14 15 17 18 19 20 21 22  * 24 25
        //   24: 15 19 20 23  * 25
        //   25: 15 19 20 23 24  *
        //
        //    Below, we list the number of adjancencies to lower-indexed nodes, to
        //    the node itself, to higher-indexed nodes, the total number of
        //    adjacencies for this node, and the location of the first and last
        //    entries required to list this set of adjacencies in a single list
        //    of all the adjacencies.
        //
        //    N   Below  Self   Above   Total First  Last
        //
        //   --      --    --      --      --   ---     0
        //    1:      0     1       5       6     1     6
        //    2:      1     1       4       6     7    12
        //    3:      2     1       9      12    13    24
        //    4:      1     1       4       6    25    30
        //    5:      2     1       6       9    31    39
        //    6:      3     1       2       6    40    45
        //    7:      4     1       4       9    46    54
        //    8:      4     1       4       9    55    63
        //    9:      4     1       4       9    62    72
        //   10:      2     1       3       6    73    78
        //   11:      6     1       5      12    79    90
        //   12:      4     1       4       9    91    99
        //   13:      9     1       9      19   100   118
        //   14:      4     1       4       9   119   127
        //   15:      5     1       6      12   128   139
        //   16:      3     1       2       6   140   145
        //   17:      4     1       4       9   146   154
        //   18:      4     1       4       9   155   163
        //   19:      4     1       4       9   164   172
        //   20:      2     1       3       6   173   178
        //   21:      6     1       2       9   179   187
        //   22:      4     1       1       6   188   193
        //   23:      9     1       2      12   194   205
        //   24:      4     1       1       6   206   211
        //   25:      5     1       0       6   212   217
        //   --      --    --      --      --   218   ---
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 August 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, int TRIANGLE_NUM, the number of triangles.
        //
        //    Input, int TRIANGLE_NODE[6*TRIANGLE_NUM], lists the nodes that
        //    make up each triangle.  The first three nodes are the vertices,
        //    in counterclockwise order.  The fourth value is the midside
        //    node between nodes 1 and 2; the fifth and sixth values are
        //    the other midside nodes in the logical order.
        //
        //    Input, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], for each side of
        //    a triangle, lists the neighboring triangle, or -1 if there is
        //    no neighbor.
        //
        //    Output, int TRIANGULATION_ORDER6_ADJ_COUNT, the number of adjacencies.
        //
        //    Output, int ADJ_COL[NODE_NUM+1].  Information about column J is stored
        //    in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
        //
    {
        int i;
        int node;
        int triangle;
        const int triangle_order = 6;

        //
        //  Set every node to be adjacent to itself.
        //
        for (node = 0; node < node_num; node++)
        {
            adj_col[node] = 1;
        }

        //
        //  Examine each triangle.
        //
        for (triangle = 0; triangle < triangle_num; triangle++)
        {
            int n1 = triangle_node[0 + triangle * triangle_order];
            int n2 = triangle_node[1 + triangle * triangle_order];
            int n3 = triangle_node[2 + triangle * triangle_order];
            int n4 = triangle_node[3 + triangle * triangle_order];
            int n5 = triangle_node[4 + triangle * triangle_order];
            int n6 = triangle_node[5 + triangle * triangle_order];
            //
            //  For sure, we add the adjacencies:
            //    43 / (34)
            //    51 / (15)
            //    54 / (45)
            //    62 / (26)
            //    64 / (46)
            //    65 / (56)
            //
            adj_col[n3 - 1] += 1;
            adj_col[n4 - 1] += 1;
            adj_col[n1 - 1] += 1;
            adj_col[n5 - 1] += 1;
            adj_col[n4 - 1] += 1;
            adj_col[n5 - 1] += 1;
            adj_col[n2 - 1] += 1;
            adj_col[n6 - 1] += 1;
            adj_col[n4 - 1] += 1;
            adj_col[n6 - 1] += 1;
            adj_col[n5 - 1] += 1;
            adj_col[n6 - 1] += 1;
            //
            //  Add edges (1,2), (1,4), (2,4) if this is the first occurrence,
            //  that is, if the edge (1,4,2) is on a boundary (TRIANGLE2 <= 0)
            //  or if this triangle is the first of the pair in which the edge
            //  occurs (TRIANGLE < TRIANGLE2).
            //
            //  Maybe add
            //    21 / 12
            //    41 / 14
            //    42 / 24
            //
            int triangle2 = triangle_neighbor[0 + triangle * 3];

            if (triangle2 < 0 || triangle < triangle2)
            {
                adj_col[n1 - 1] += 1;
                adj_col[n2 - 1] += 1;
                adj_col[n1 - 1] += 1;
                adj_col[n4 - 1] += 1;
                adj_col[n2 - 1] += 1;
                adj_col[n4 - 1] += 1;
            }

            //
            //  Maybe add
            //    32 / 23
            //    52 / 25
            //    53 / 35
            //
            triangle2 = triangle_neighbor[1 + triangle * 3];

            if (triangle2 < 0 || triangle < triangle2)
            {
                adj_col[n2 - 1] += 1;
                adj_col[n3 - 1] += 1;
                adj_col[n2 - 1] += 1;
                adj_col[n5 - 1] += 1;
                adj_col[n3 - 1] += 1;
                adj_col[n5 - 1] += 1;
            }

            //
            //  Maybe add
            //    31 / 13
            //    61 / 16
            //    63 / 36
            //
            triangle2 = triangle_neighbor[2 + triangle * 3];

            if (triangle2 < 0 || triangle < triangle2)
            {
                adj_col[n1 - 1] += 1;
                adj_col[n3 - 1] += 1;
                adj_col[n1 - 1] += 1;
                adj_col[n6 - 1] += 1;
                adj_col[n3 - 1] += 1;
                adj_col[n6 - 1] += 1;
            }
        }

        //
        //  We used ADJ_COL to count the number of entries in each column.
        //  Convert it to pointers into the ADJ array.
        //
        for (node = node_num; 1 <= node; node--)
        {
            adj_col[node] = adj_col[node - 1];
        }

        adj_col[0] = 1;
        for (i = 1; i <= node_num; i++)
        {
            adj_col[i] = adj_col[i - 1] + adj_col[i];
        }

        int adj_num = adj_col[node_num] - 1;

        return adj_num;
    }

    public static int[] triangulation_order6_adj_set(int node_num, int triangle_num,
            int[] triangle_node, int[] triangle_neighbor, int adj_num, int[] adj_col )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULATION_ORDER6_ADJ_SET sets adjacencies in a triangulation.
        //
        //  Discussion:
        //
        //    This routine is called to count the adjacencies, so that the
        //    appropriate amount of memory can be set aside for storage when
        //    the adjacency structure is created.
        //
        //    The triangulation is assumed to involve 6-node triangles.
        //
        //    Two nodes are "adjacent" if they are both nodes in some triangle.
        //    Also, a node is considered to be adjacent to itself.
        //
        //    This routine can be used to create the compressed column storage
        //    for a quadratic triangle finite element discretization of
        //    Poisson's equation in two dimensions.
        //
        //  Diagram:
        //
        //       3
        //    s  |.
        //    i  | .
        //    d  |  .
        //    e  6   5  side 2
        //       |    .
        //    3  |     .
        //       |      .
        //       1---4---2
        //
        //         side 1
        //
        //    The local node numbering
        //
        //
        //   21-22-23-24-25
        //    |.    |.    |
        //    | .   | .   |
        //   16 17 18 19 20
        //    |   . |   . |
        //    |    .|    .|
        //   11-12-13-14-15
        //    |.    |.    |
        //    | .   | .   |
        //    6  7  8  9 10
        //    |   . |   . |
        //    |    .|    .|
        //    1--2--3--4--5
        //
        //    A sample grid.
        //
        //
        //    Below, we have a chart that lists the nodes adjacent to each node, with
        //    an asterisk to indicate the adjacency of the node to itself
        //    (in some cases, you want to count this self adjacency and in some
        //    you don't).
        //
        //    N   Adjacencies
        //
        //    1:  *  2  3  6  7 11
        //    2:  1  *  3  6  7 11
        //    3:  1  2  *  4  5  6  7  8  9 11 12 13
        //    4:  3  *  5  8  9 13
        //    5:  3  4  *  8  9 10 13 14 15
        //    6:  1  2  3  *  7 11
        //    7:  1  2  3  6  *  8 11 12 13
        //    8:  3  4  5  7  *  9 11 12 13
        //    9:  3  4  5  8  * 10 13 14 15
        //   10:  5  9  * 13 14 15
        //   11:  1  2  3  6  7  8  * 12 13 16 17 21
        //   12:  3  7  8 11  * 13 16 17 21
        //   13:  3  4  5  7  8  9 10 11 12  * 14 15 16 17 18 19 21 22 23
        //   14:  5  9 10 13  * 15 18 19 23
        //   15:  5  9 10 13 14  * 18 19 20 23 24 25
        //   16: 11 12 13  * 17 21
        //   17: 11 12 13 16  * 18 21 22 23
        //   18: 13 14 15 17  * 19 21 22 23
        //   19: 13 14 15 18  * 20 23 24 25
        //   20: 15 19  * 23 24 25
        //   21: 11 12 13 16 17 18  * 22 23
        //   22: 13 17 18 21  * 23
        //   23: 13 14 15 17 18 19 20 21 22  * 24 25
        //   24: 15 19 20 23  * 25
        //   25: 15 19 20 23 24  *
        //
        //    Below, we list the number of adjancencies to lower-indexed nodes, to
        //    the node itself, to higher-indexed nodes, the total number of
        //    adjacencies for this node, and the location of the first and last
        //    entries required to list this set of adjacencies in a single list
        //    of all the adjacencies.
        //
        //    N   Below  Self   Above   Total First  Last
        //
        //   --      --    --      --      --   ---     0
        //    1:      0     1       5       6     1     6
        //    2:      1     1       4       6     7    12
        //    3:      2     1       9      12    13    24
        //    4:      1     1       4       6    25    30
        //    5:      2     1       6       9    31    39
        //    6:      3     1       2       6    40    45
        //    7:      4     1       4       9    46    54
        //    8:      4     1       4       9    55    63
        //    9:      4     1       4       9    62    72
        //   10:      2     1       3       6    73    78
        //   11:      6     1       5      12    79    90
        //   12:      4     1       4       9    91    99
        //   13:      9     1       9      19   100   118
        //   14:      4     1       4       9   119   127
        //   15:      5     1       6      12   128   139
        //   16:      3     1       2       6   140   145
        //   17:      4     1       4       9   146   154
        //   18:      4     1       4       9   155   163
        //   19:      4     1       4       9   164   172
        //   20:      2     1       3       6   173   178
        //   21:      6     1       2       9   179   187
        //   22:      4     1       1       6   188   193
        //   23:      9     1       2      12   194   205
        //   24:      4     1       1       6   206   211
        //   25:      5     1       0       6   212   217
        //   --      --    --      --      --   218   ---
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 August 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, int TRIANGLE_NUM, the number of triangles.
        //
        //    Input, int TRIANGLE_NODE[6*TRIANGLE_NUM], lists the nodes that
        //    make up each triangle.  The first three nodes are the vertices,
        //    in counterclockwise order.  The fourth value is the midside
        //    node between nodes 1 and 2; the fifth and sixth values are
        //    the other midside nodes in the logical order.
        //
        //    Input, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], for each side of
        //    a triangle, lists the neighboring triangle, or -1 if there is
        //    no neighbor.
        //
        //    Input, int ADJ_NUM, the number of adjacencies.
        //
        //    Input, int ADJ_COL[NODE_NUM+1].  Information about column J is stored
        //    in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
        //
        //    Output, int TRIANGULATION_ORDER6_ADJ_SET[ADJ_NUM], the adjacency
        //    information.
        //
    {
        int k;
        int node;
        int triangle;
        const int triangle_order = 6;

        int[] adj = new int[adj_num];
        for (k = 0; k < adj_num; k++)
        {
            adj[k] = -1;
        }

        int[] adj_copy = new int[node_num];
        for (node = 0; node < node_num; node++)
        {
            adj_copy[node] = adj_col[node];
        }

        //
        //  Set every node to be adjacent to itself.
        //
        for (node = 1; node <= node_num; node++)
        {
            adj[adj_copy[node - 1] - 1] = node;
            adj_copy[node - 1] += 1;
        }

        //
        //  Examine each triangle.
        //
        for (triangle = 0; triangle < triangle_num; triangle++)
        {
            int n1 = triangle_node[0 + triangle * triangle_order];
            int n2 = triangle_node[1 + triangle * triangle_order];
            int n3 = triangle_node[2 + triangle * triangle_order];
            int n4 = triangle_node[3 + triangle * triangle_order];
            int n5 = triangle_node[4 + triangle * triangle_order];
            int n6 = triangle_node[5 + triangle * triangle_order];
            //
            //  For sure, we add the adjacencies:
            //    43 / (34)
            //    51 / (15)
            //    54 / (45)
            //    62 / (26)
            //    64 / (46)
            //    65 / (56)
            //
            adj[adj_copy[n3 - 1] - 1] = n4;
            adj_copy[n3 - 1] += 1;
            adj[adj_copy[n4 - 1] - 1] = n3;
            adj_copy[n4 - 1] += 1;

            adj[adj_copy[n1 - 1] - 1] = n5;
            adj_copy[n1 - 1] += 1;
            adj[adj_copy[n5 - 1] - 1] = n1;
            adj_copy[n5 - 1] += 1;

            adj[adj_copy[n4 - 1] - 1] = n5;
            adj_copy[n4 - 1] += 1;
            adj[adj_copy[n5 - 1] - 1] = n4;
            adj_copy[n5 - 1] += 1;

            adj[adj_copy[n2 - 1] - 1] = n6;
            adj_copy[n2 - 1] += 1;
            adj[adj_copy[n6 - 1] - 1] = n2;
            adj_copy[n6 - 1] += 1;

            adj[adj_copy[n4 - 1] - 1] = n6;
            adj_copy[n4 - 1] += 1;
            adj[adj_copy[n6 - 1] - 1] = n4;
            adj_copy[n6 - 1] += 1;

            adj[adj_copy[n5 - 1] - 1] = n6;
            adj_copy[n5 - 1] += 1;
            adj[adj_copy[n6 - 1] - 1] = n5;
            adj_copy[n6 - 1] += 1;
            //
            //  Add edges (1,2), (1,4), (2,4) if this is the first occurrence,
            //  that is, if the edge (1,4,2) is on a boundary (TRIANGLE2 <= 0)
            //  or if this triangle is the first of the pair in which the edge
            //  occurs (TRIANGLE < TRIANGLE2).
            //
            //  Maybe add
            //    21 / 12
            //    41 / 14
            //    42 / 24
            //
            int triangle2 = triangle_neighbor[0 + triangle * 3];

            if (triangle2 < 0 || triangle < triangle2)
            {
                adj[adj_copy[n1 - 1] - 1] = n2;
                adj_copy[n1 - 1] += 1;
                adj[adj_copy[n2 - 1] - 1] = n1;
                adj_copy[n2 - 1] += 1;
                adj[adj_copy[n1 - 1] - 1] = n4;
                adj_copy[n1 - 1] += 1;
                adj[adj_copy[n4 - 1] - 1] = n1;
                adj_copy[n4 - 1] += 1;
                adj[adj_copy[n2 - 1] - 1] = n4;
                adj_copy[n2 - 1] += 1;
                adj[adj_copy[n4 - 1] - 1] = n2;
                adj_copy[n4 - 1] += 1;
            }

            //
            //  Maybe add
            //    32 / 23
            //    52 / 25
            //    53 / 35
            //
            triangle2 = triangle_neighbor[1 + triangle * 3];

            if (triangle2 < 0 || triangle < triangle2)
            {
                adj[adj_copy[n2 - 1] - 1] = n3;
                adj_copy[n2 - 1] += 1;
                adj[adj_copy[n3 - 1] - 1] = n2;
                adj_copy[n3 - 1] += 1;
                adj[adj_copy[n2 - 1] - 1] = n5;
                adj_copy[n2 - 1] += 1;
                adj[adj_copy[n5 - 1] - 1] = n2;
                adj_copy[n5 - 1] += 1;
                adj[adj_copy[n3 - 1] - 1] = n5;
                adj_copy[n3 - 1] += 1;
                adj[adj_copy[n5 - 1] - 1] = n3;
                adj_copy[n5 - 1] += 1;
            }

            //
            //  Maybe add
            //    31 / 13
            //    61 / 16
            //    63 / 36
            //
            triangle2 = triangle_neighbor[2 + triangle * 3];

            if (triangle2 >= 0 && triangle >= triangle2)
            {
                continue;
            }

            adj[adj_copy[n1 - 1] - 1] = n3;
            adj_copy[n1 - 1] += 1;
            adj[adj_copy[n3 - 1] - 1] = n1;
            adj_copy[n3 - 1] += 1;
            adj[adj_copy[n1 - 1] - 1] = n6;
            adj_copy[n1 - 1] += 1;
            adj[adj_copy[n6 - 1] - 1] = n1;
            adj_copy[n6 - 1] += 1;
            adj[adj_copy[n3 - 1] - 1] = n6;
            adj_copy[n3 - 1] += 1;
            adj[adj_copy[n6 - 1] - 1] = n3;
            adj_copy[n6 - 1] += 1;
        }

        //
        //  Ascending sort the entries for each node.
        //
        for (node = 1; node <= node_num; node++)
        {
            int k1 = adj_col[node - 1];
            int k2 = adj_col[node] - 1;
            typeMethods.i4vec_sort_heap_a(k2 + 1 - k1, ref adj, aIndex: + k1 - 1);
        }

        return adj;
    }
}