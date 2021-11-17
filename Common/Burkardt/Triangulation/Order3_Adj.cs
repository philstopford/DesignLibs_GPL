using Burkardt.Types;

namespace Burkardt.TriangulationNS;

public static partial class Adjacency
{
    public static int triangulation_order3_adj_count(int node_num, int triangle_num,
            int[] triangle_node, int[] triangle_neighbor, int[] adj_col )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULATION_ORDER3_ADJ_COUNT counts adjacencies in a triangulation.
        //
        //  Discussion:
        //
        //    This routine is called to count the adjacencies, so that the
        //    appropriate amount of memory can be set aside for storage when
        //    the adjacency structure is created.
        //
        //    The triangulation is assumed to involve 3-node triangles.
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
        //    e  |   .  side 2
        //       |    .
        //    3  |     .
        //       |      .
        //       1-------2
        //
        //         side 1
        //
        //    The local node numbering
        //
        //
        //   21-22-23-24-25
        //    |. |. |. |. |
        //    | .| .| .| .|
        //   16-17-18-19-20
        //    |. |. |. |. |
        //    | .| .| .| .|
        //   11-12-13-14-15
        //    |. |. |. |. |
        //    | .| .| .| .|
        //    6--7--8--9-10
        //    |. |. |. |. |
        //    | .| .| .| .|
        //    1--2--3--4--5
        //
        //    A sample grid.
        //
        //
        //    Below, we have a chart that summarizes the adjacency relationships
        //    in the sample grid.  On the left, we list the node, and its neighbors,
        //    with an asterisk to indicate the adjacency of the node to itself
        //    (in some cases, you want to count this self adjacency and in some
        //    you don't).  On the right, we list the number of adjancencies to
        //    lower-indexed nodes, to the node itself, to higher-indexed nodes,
        //    the total number of adjacencies for this node, and the location
        //    of the first and last entries required to list this set of adjacencies
        //    in a single list of all the adjacencies.
        //
        //    N   Adjacencies                Below  Self   Above   Total First  Last
        //
        //   --  -- -- -- -- -- -- --           --    --      --      --   ---     0
        //    1:  *  2  6                        0     1       2       3     1     3
        //    2:  1  *  3  6  7                  1     1       3       5     4     8
        //    3:  2  *  4  7  8                  1     1       3       5     9    13
        //    4:  3  *  5  8  9                  1     1       3       5    14    18
        //    5:  4  *  9 10                     1     1       2       4    19    22
        //    6:  1  2  *  7 11                  2     1       2       5    23    27
        //    7:  2  3  6  *  8 11 12            3     1       3       7    28    34
        //    8:  3  4  7  *  9 12 13            3     1       3       7    35    41
        //    9:  4  5  8  * 10 13 14            3     1       3       7    42    48
        //   10:  5  9  * 14 15                  2     1       2       5    49    53
        //   11:  6  7  * 12 16                  2     1       2       5    54    58
        //   12:  7  8 11  * 13 16 17            3     1       3       7    59    65
        //   13:  8  9 12  * 14 17 18            3     1       3       7    66    72
        //   14:  9 10 13  * 15 18 19            3     1       3       7    73    79
        //   15: 10 14  * 19 20                  2     1       2       5    80    84
        //   16: 11 12  * 17 21                  2     1       2       5    85    89
        //   17: 12 13 16  * 18 21 22            3     1       3       7    90    96
        //   18: 13 14 17  * 19 22 23            3     1       3       7    97   103
        //   19: 14 15 18  * 20 23 24            3     1       3       7   104   110
        //   20: 15 19  * 24 25                  2     1       2       5   111   115
        //   21: 16 17  * 22                     2     1       1       4   116   119
        //   22: 17 18 21  * 23                  3     1       1       5   120   124
        //   23: 18 19 22  * 24                  3     1       1       5   125   129
        //   24: 19 20 23  * 25                  3     1       1       5   130   134
        //   25: 20 24  *                        2     1       0       3   135   137
        //   --  -- -- -- -- -- -- --           --    --      --      --   138   ---
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
        //    Input, int TRIANGLE_NODE[3*TRIANGLE_NUM], lists the nodes that
        //    make up each triangle, in counterclockwise order.
        //
        //    Input, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], for each side of
        //    a triangle, lists the neighboring triangle, or -1 if there is
        //    no neighbor.
        //
        //    Output, TRIANGULATION_ORDER3_ADJ_COUNT, the number of adjacencies.
        //
        //    Output, int ADJ_COL[NODE_NUM+1].  Information about column J is stored
        //    in entries ADJ_COL(J) through ADJ_COL(J+1)-1 of ADJ.
        //
    {
        int adj_num;
        int i;
        int n1;
        int n2;
        int n3;
        int node;
        int triangle;
        int triangle_order = 3;
        int triangle2;

        adj_num = 0;
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
            n1 = triangle_node[0 + triangle * triangle_order];
            n2 = triangle_node[1 + triangle * triangle_order];
            n3 = triangle_node[2 + triangle * triangle_order];
            //
            //  Add edge (1,2) if this is the first occurrence,
            //  that is, if the edge (1,2) is on a boundary (TRIANGLE2 <= 0)
            //  or if this triangle is the first of the pair in which the edge
            //  occurs (TRIANGLE < TRIANGLE2).
            //
            triangle2 = triangle_neighbor[0 + triangle * 3];

            if (triangle2 < 0 || triangle < triangle2)
            {
                adj_col[n1 - 1] += 1;
                adj_col[n2 - 1] += 1;
            }

            //
            //  Add edge (2,3).
            //
            triangle2 = triangle_neighbor[1 + triangle * 3];

            if (triangle2 < 0 || triangle < triangle2)
            {
                adj_col[n2 - 1] += 1;
                adj_col[n3 - 1] += 1;
            }

            //
            //  Add edge (3,1).
            //
            triangle2 = triangle_neighbor[2 + triangle * 3];

            if (triangle2 < 0 || triangle < triangle2)
            {
                adj_col[n1 - 1] += 1;
                adj_col[n3 - 1] += 1;
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

        adj_num = adj_col[node_num] - 1;

        return adj_num;
    }

    public static int[] triangulation_order3_adj_set(int node_num, int triangle_num,
            int[] triangle_node, int[] triangle_neighbor, int adj_num, int[] adj_col )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULATION_ORDER3_ADJ_SET sets adjacencies in a triangulation.
        //
        //  Discussion:
        //
        //    This routine is called to set the adjacencies, after the
        //    appropriate amount of memory has been set aside for storage.
        //
        //    The triangulation is assumed to involve 3-node triangles.
        //
        //    Two nodes are "adjacent" if they are both nodes in some triangle.
        //    Also, a node is considered to be adjacent to itself.
        //
        //    This routine can be used to create the compressed column storage
        //    for a linear triangle finite element discretization of
        //    Poisson's equation in two dimensions.
        //
        //  Diagram:
        //
        //       3
        //    s  |.
        //    i  | .
        //    d  |  .
        //    e  |   .  side 2
        //       |    .
        //    3  |     .
        //       |      .
        //       1-------2
        //
        //         side 1
        //
        //    The local node numbering
        //
        //
        //   21-22-23-24-25
        //    |. |. |. |. |
        //    | .| .| .| .|
        //   16-17-18-19-20
        //    |. |. |. |. |
        //    | .| .| .| .|
        //   11-12-13-14-15
        //    |. |. |. |. |
        //    | .| .| .| .|
        //    6--7--8--9-10
        //    |. |. |. |. |
        //    | .| .| .| .|
        //    1--2--3--4--5
        //
        //    A sample grid
        //
        //
        //    Below, we have a chart that summarizes the adjacency relationships
        //    in the sample grid.  On the left, we list the node, and its neighbors,
        //    with an asterisk to indicate the adjacency of the node to itself
        //    (in some cases, you want to count this self adjacency and in some
        //    you don't).  On the right, we list the number of adjancencies to
        //    lower-indexed nodes, to the node itself, to higher-indexed nodes,
        //    the total number of adjacencies for this node, and the location
        //    of the first and last entries required to list this set of adjacencies
        //    in a single list of all the adjacencies.
        //
        //    N   Adjacencies                Below  Self    Above  Total First  Last
        //
        //   --  -- -- -- -- -- -- --           --    --      --      --   ---     0
        //    1:  *  2  6                        0     1       2       3     1     3
        //    2:  1  *  3  6  7                  1     1       3       5     4     8
        //    3:  2  *  4  7  8                  1     1       3       5     9    13
        //    4:  3  *  5  8  9                  1     1       3       5    14    18
        //    5:  4  *  9 10                     1     1       2       4    19    22
        //    6:  1  2  *  7 11                  2     1       2       5    23    27
        //    7:  2  3  6  *  8 11 12            3     1       3       7    28    34
        //    8:  3  4  7  *  9 12 13            3     1       3       7    35    41
        //    9:  4  5  8  * 10 13 14            3     1       3       7    42    48
        //   10:  5  9  * 14 15                  2     1       2       5    49    53
        //   11:  6  7  * 12 16                  2     1       2       5    54    58
        //   12:  7  8 11  * 13 16 17            3     1       3       7    59    65
        //   13:  8  9 12  * 14 17 18            3     1       3       7    66    72
        //   14:  9 10 13  * 15 18 19            3     1       3       7    73    79
        //   15: 10 14  * 19 20                  2     1       2       5    80    84
        //   16: 11 12  * 17 21                  2     1       2       5    85    89
        //   17: 12 13 16  * 18 21 22            3     1       3       7    90    96
        //   18: 13 14 17  * 19 22 23            3     1       3       7    97   103
        //   19: 14 15 18  * 20 23 24            3     1       3       7   104   110
        //   20: 15 19  * 24 25                  2     1       2       5   111   115
        //   21: 16 17  * 22                     2     1       1       4   116   119
        //   22: 17 18 21  * 23                  3     1       1       5   120   124
        //   23: 18 19 22  * 24                  3     1       1       5   125   129
        //   24: 19 20 23  * 25                  3     1       1       5   130   134
        //   25: 20 24  *                        2     1       0       3   135   137
        //   --  -- -- -- -- -- -- --           --    --      --      --   138   ---
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
        //    Input, int TRIANGLE_NODE[3*TRIANGLE_NUM], lists the nodes that
        //    make up each triangle in counterclockwise order.
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
        //    Output, int TRIANGULATION_ORDER3_ADJ_SET[ADJ_NUM], the adjacency
        //    information.
        //
    {
        int[] adj;
        int[] adj_copy;
        int k;
        int k1;
        int k2;
        int n1;
        int n2;
        int n3;
        int node;
        int triangle;
        int triangle2;
        int triangle_order = 3;

        adj = new int[adj_num];
        for (k = 0; k < adj_num; k++)
        {
            adj[k] = -1;
        }

        adj_copy = new int[node_num];
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
            n1 = triangle_node[0 + triangle * triangle_order];
            n2 = triangle_node[1 + triangle * triangle_order];
            n3 = triangle_node[2 + triangle * triangle_order];
            //
            //  Add edge (1,2) if this is the first occurrence,
            //  that is, if the edge (1,2) is on a boundary (TRIANGLE2 <= 0)
            //  or if this triangle is the first of the pair in which the edge
            //  occurs (TRIANGLE < TRIANGLE2).
            //
            triangle2 = triangle_neighbor[0 + triangle * 3];

            if (triangle2 < 0 || triangle < triangle2)
            {
                adj[adj_copy[n1 - 1] - 1] = n2;
                adj_copy[n1 - 1] += 1;
                adj[adj_copy[n2 - 1] - 1] = n1;
                adj_copy[n2 - 1] += 1;
            }

            //
            //  Add edge (2,3).
            //
            triangle2 = triangle_neighbor[1 + triangle * 3];

            if (triangle2 < 0 || triangle < triangle2)
            {
                adj[adj_copy[n2 - 1] - 1] = n3;
                adj_copy[n2 - 1] += 1;
                adj[adj_copy[n3 - 1] - 1] = n2;
                adj_copy[n3 - 1] += 1;
            }

            //
            //  Add edge (3,1).
            //
            triangle2 = triangle_neighbor[2 + triangle * 3];

            if (triangle2 < 0 || triangle < triangle2)
            {
                adj[adj_copy[n1 - 1] - 1] = n3;
                adj_copy[n1 - 1] += 1;
                adj[adj_copy[n3 - 1] - 1] = n1;
                adj_copy[n3 - 1] += 1;
            }
        }

        //
        //  Ascending sort the entries for each node.
        //
        for (node = 1; node <= node_num; node++)
        {
            k1 = adj_col[node - 1];
            k2 = adj_col[node] - 1;
            typeMethods.i4vec_sort_heap_a(k2 + 1 - k1, ref adj, aIndex: + k1 - 1);
        }
            
        return adj;
    }

    public static void triangulation_order3_adj_set2(int node_num, int triangle_num,
            int[] triangle_node, int[] triangle_neighbor, int adj_num, int[] adj_col,
            int[] ia, int[] ja )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULATION_ORDER3_ADJ_SET2 sets adjacencies in a triangulation.
        //
        //  Discussion:
        //
        //    This routine is called to set up the arrays IA and JA that
        //    record which nodes are adjacent in a triangulation.
        //
        //    The triangulation is assumed to involve 3-node triangles.
        //
        //    Two nodes are "adjacent" if they are both nodes in some triangle.
        //    Also, a node is considered to be adjacent to itself.
        //
        //    This routine can be used to create the compressed column storage
        //    for a linear triangle finite element discretization of
        //    Poisson's equation in two dimensions.
        //
        //  Diagram:
        //
        //       3
        //    s  |.
        //    i  | .
        //    d  |  .
        //    e  |   .  side 2
        //       |    .
        //    3  |     .
        //       |      .
        //       1-------2
        //
        //         side 1
        //
        //    The local node numbering
        //
        //
        //   21-22-23-24-25
        //    |. |. |. |. |
        //    | .| .| .| .|
        //   16-17-18-19-20
        //    |. |. |. |. |
        //    | .| .| .| .|
        //   11-12-13-14-15
        //    |. |. |. |. |
        //    | .| .| .| .|
        //    6--7--8--9-10
        //    |. |. |. |. |
        //    | .| .| .| .|
        //    1--2--3--4--5
        //
        //    A sample grid
        //
        //
        //    Below, we have a chart that summarizes the adjacency relationships
        //    in the sample grid.  On the left, we list the node, and its neighbors,
        //    with an asterisk to indicate the adjacency of the node to itself
        //    (in some cases, you want to count this self adjacency and in some
        //    you don't).  On the right, we list the number of adjancencies to
        //    lower-indexed nodes, to the node itself, to higher-indexed nodes,
        //    the total number of adjacencies for this node, and the location
        //    of the first and last entries required to list this set of adjacencies
        //    in a single list of all the adjacencies.
        //
        //    N   Adjacencies                Below  Self    Above  Total First  Last
        //
        //   --  -- -- -- -- -- -- --           --    --      --      --   ---     0
        //    1:  *  2  6                        0     1       2       3     1     3
        //    2:  1  *  3  6  7                  1     1       3       5     4     8
        //    3:  2  *  4  7  8                  1     1       3       5     9    13
        //    4:  3  *  5  8  9                  1     1       3       5    14    18
        //    5:  4  *  9 10                     1     1       2       4    19    22
        //    6:  1  2  *  7 11                  2     1       2       5    23    27
        //    7:  2  3  6  *  8 11 12            3     1       3       7    28    34
        //    8:  3  4  7  *  9 12 13            3     1       3       7    35    41
        //    9:  4  5  8  * 10 13 14            3     1       3       7    42    48
        //   10:  5  9  * 14 15                  2     1       2       5    49    53
        //   11:  6  7  * 12 16                  2     1       2       5    54    58
        //   12:  7  8 11  * 13 16 17            3     1       3       7    59    65
        //   13:  8  9 12  * 14 17 18            3     1       3       7    66    72
        //   14:  9 10 13  * 15 18 19            3     1       3       7    73    79
        //   15: 10 14  * 19 20                  2     1       2       5    80    84
        //   16: 11 12  * 17 21                  2     1       2       5    85    89
        //   17: 12 13 16  * 18 21 22            3     1       3       7    90    96
        //   18: 13 14 17  * 19 22 23            3     1       3       7    97   103
        //   19: 14 15 18  * 20 23 24            3     1       3       7   104   110
        //   20: 15 19  * 24 25                  2     1       2       5   111   115
        //   21: 16 17  * 22                     2     1       1       4   116   119
        //   22: 17 18 21  * 23                  3     1       1       5   120   124
        //   23: 18 19 22  * 24                  3     1       1       5   125   129
        //   24: 19 20 23  * 25                  3     1       1       5   130   134
        //   25: 20 24  *                        2     1       0       3   135   137
        //   --  -- -- -- -- -- -- --           --    --      --      --   138   ---
        //
        //    For this example, the initial portion of the IA and JA arrays will be:
        //
        //      (1,1), (1,2), (1,6),
        //      (2,1), (2,2), (2,3), (2,6), (2,7),
        //      (3,2), (3,3), (3,4), (3,7), (3,8),
        //     ...
        //      (25,20), (25,24), (25,25)
        //
        //    for a total of 137 pairs of values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 July 2007
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
        //    Input, int TRIANGLE_NODE[3*TRIANGLE_NUM], lists the nodes that
        //    make up each triangle in counterclockwise order.
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
        //    Output, int IA[ADJ_NUM], JA[ADJ_NUM], the adjacency information.
        //
    {
        int adj;
        int[] adj_copy;
        int n1;
        int n2;
        int n3;
        int node;
        int triangle;
        int triangle2;
        int triangle_order = 3;

        for (adj = 0; adj < adj_num; adj++)
        {
            ia[adj] = -1;
        }

        for (adj = 0; adj < adj_num; adj++)
        {
            ja[adj] = -1;
        }

        adj_copy = new int[node_num];
        for (node = 0; node < node_num; node++)
        {
            adj_copy[node] = adj_col[node];
        }

        //
        //  Set every node to be adjacent to itself.
        //
        for (node = 1; node <= node_num; node++)
        {
            ia[adj_copy[node - 1] - 1] = node;
            ja[adj_copy[node - 1] - 1] = node;
            adj_copy[node - 1] += 1;
        }

        //
        //  Examine each triangle.
        //
        for (triangle = 0; triangle < triangle_num; triangle++)
        {
            n1 = triangle_node[0 + triangle * triangle_order];
            n2 = triangle_node[1 + triangle * triangle_order];
            n3 = triangle_node[2 + triangle * triangle_order];
            //
            //  Add edge (1,2) if this is the first occurrence,
            //  that is, if the edge (1,2) is on a boundary (TRIANGLE2 <= 0)
            //  or if this triangle is the first of the pair in which the edge
            //  occurs (TRIANGLE < TRIANGLE2).
            //
            triangle2 = triangle_neighbor[0 + triangle * 3];

            if (triangle2 < 0 || triangle < triangle2)
            {
                ia[adj_copy[n1 - 1] - 1] = n1;
                ja[adj_copy[n1 - 1] - 1] = n2;
                adj_copy[n1 - 1] += 1;

                ia[adj_copy[n2 - 1] - 1] = n2;
                ja[adj_copy[n2 - 1] - 1] = n1;
                adj_copy[n2 - 1] += 1;
            }

            //
            //  Add edge (2,3).
            //
            triangle2 = triangle_neighbor[1 + triangle * 3];

            if (triangle2 < 0 || triangle < triangle2)
            {
                ia[adj_copy[n2 - 1] - 1] = n2;
                ja[adj_copy[n2 - 1] - 1] = n3;
                adj_copy[n2 - 1] += 1;

                ia[adj_copy[n3 - 1] - 1] = n3;
                ja[adj_copy[n3 - 1] - 1] = n2;
                adj_copy[n3 - 1] += 1;
            }

            //
            //  Add edge (3,1).
            //
            triangle2 = triangle_neighbor[2 + triangle * 3];

            if (triangle2 < 0 || triangle < triangle2)
            {
                ia[adj_copy[n1 - 1] - 1] = n1;
                ja[adj_copy[n1 - 1] - 1] = n3;
                adj_copy[n1 - 1] += 1;

                ia[adj_copy[n3 - 1] - 1] = n3;
                ja[adj_copy[n3 - 1] - 1] = n1;
                adj_copy[n3 - 1] += 1;
            }
        }

        //
        //  Lexically sort the IA, JA values.
        //
        typeMethods.i4vec2_sort_a(adj_num, ref ia, ref ja);
    }

    public static int[] triangulation_order3_adjacency(int node_num, int element_num,
            int[] element_node)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULATION_ORDER3_ADJACENCY computes the full adjacency matrix
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 March 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NODE_NUM, the number of nodes in the
        //    triangulation.
        //
        //    Input, int ELEMENT_NUM, the number of triangles in
        //    the triangulation.
        //
        //    Input, int ELEMENT_NODE[3*ELEMENT_NUM],
        //    the nodes making up each triangle.
        //
        //    Output, int TRIANGULATION_ORDER3_ADJACENCY[NODE_NUM*NODE_NUM], the adjacency
        //    matrix.  ADJ(I,J) is 1 if nodes I and J are adjacent, that is,
        //    they are immediate neighbors on an edge of the triangulation.
        //
    {
        int[] adj;
        int element;
        int i;
        int j;
        int k;

        adj = new int[node_num * node_num];

        for (j = 0; j < node_num; j++)
        {
            for (i = 0; i < node_num; i++)
            {
                adj[i + j * node_num] = 0;
            }
        }

        for (element = 0; element < element_num; element++)
        {
            i = element_node[0 + element * 3];
            j = element_node[1 + element * 3];
            k = element_node[2 + element * 3];

            adj[i + j * node_num] = 1;
            adj[i + k * node_num] = 1;
            adj[j + i * node_num] = 1;
            adj[j + k * node_num] = 1;
            adj[k + i * node_num] = 1;
            adj[k + j * node_num] = 1;
        }

        return adj;
    }
}