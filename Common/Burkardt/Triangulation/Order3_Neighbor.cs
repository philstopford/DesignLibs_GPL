using System;
using System.Globalization;
using Burkardt.Types;

namespace Burkardt.TriangulationNS;

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
        int t;

        int n1 = triangle_node[s1 - 1 + (t1 - 1) * 3];
        int ss = typeMethods.i4_wrap(s1 + 1, 1, 3);
        int n2 = triangle_node[ss - 1 + (t1 - 1) * 3];

        for (t = 0; t < triangle_num; t++)
        {
            int s;
            for (s = 0; s < 3; s++)
            {
                if (triangle_node[s + t * 3] != n1)
                {
                    continue;
                }

                ss = typeMethods.i4_wrap(s - 1, 0, 2);
                if (triangle_node[ss + t * 3] != n2)
                {
                    continue;
                }

                t2 = t + 1;
                s2 = ss + 1;
                return;
            }
        }

        t2 = -1;
        s2 = -1;
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
        int nabe;
        int tri;

        nabes = new int[nabes_max];
        int[] nabes1 = new int[nabes_max];
        //
        //  Step 1.  From the triangle list (I,J,K)
        //  construct the neighbor relations: (I,J), (J,K), (K,I), (J,I), (K,J), (I,K).
        //
        int n = 0;

        for (tri = 0; tri < triangle_num; tri++)
        {
            i = triangle_node[0 + tri * 3];
            int j = triangle_node[1 + tri * 3];
            int k = triangle_node[2 + tri * 3];
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

            n += 6;
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

        int i_current = 0;

        for (nabe = 1; nabe <= n; nabe++)
        {
            i = nabes1[nabe - 1];
            if (i == i_current)
            {
                nabes_num[i - 1] += 1;
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

        Console.WriteLine("");
        Console.WriteLine("  Node Nabes Index  List");
        Console.WriteLine("");
        string cout = "";

        for (i = 0; i < node_num; i++)
        {
            cout = i.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                                           + nabes_num[i].ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                                           + nabes_first[i].ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  ";

            int k = 0;
            int j;
            for (j = nabes_first[i] - 1; j < nabes_first[i] + nabes_num[i]; j++)
            {
                switch (k)
                {
                    case 10:
                        Console.WriteLine(cout);
                        cout = "                  ";
                        k = 0;
                        break;
                }

                cout += nabes[j].ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  ";
                k += 1;
            }
        }

        Console.WriteLine(cout);
    }

    public static int[] triangulation_order3_neighbor_triangles(int triangle_num,
            int[] triangle_node)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULATION_ORDER3_NEIGHBOR_TRIANGLES determines triangle neighbors.
        //
        //  Discussion:
        //
        //    A triangulation of a set of nodes can be completely described by
        //    the coordinates of the nodes, and the list of nodes that make up
        //    each triangle.  However, in some cases, it is necessary to know
        //    triangle adjacency information, that is, which triangle, if any,
        //    is adjacent to a given triangle on a particular side.
        //
        //    This routine creates a data structure recording this information.
        //
        //    The primary amount of work occurs in sorting a list of 3 * TRIANGLE_NUM
        //    data items.
        //
        //    This routine was modified to work with columns rather than rows.
        //
        //  Example:
        //
        //    The input information from TRIANGLE_NODE:
        //
        //    Triangle   Nodes
        //    --------   ---------------
        //     1         3      4      1
        //     2         3      1      2
        //     3         3      2      8
        //     4         2      1      5
        //     5         8      2     13
        //     6         8     13      9
        //     7         3      8      9
        //     8        13      2      5
        //     9         9     13      7
        //    10         7     13      5
        //    11         6      7      5
        //    12         9      7      6
        //    13        10      9      6
        //    14         6      5     12
        //    15        11      6     12
        //    16        10      6     11
        //
        //    The output information in TRIANGLE_NEIGHBOR:
        //
        //    Triangle  Neighboring Triangles
        //    --------  ---------------------
        //
        //     1        -1     -1      2
        //     2         1      4      3
        //     3         2      5      7
        //     4         2     -1      8
        //     5         3      8      6
        //     6         5      9      7
        //     7         3      6     -1
        //     8         5      4     10
        //     9         6     10     12
        //    10         9      8     11
        //    11        12     10     14
        //    12         9     11     13
        //    13        -1     12     16
        //    14        11     -1     15
        //    15        16     14     -1
        //    16        13     15     -1
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int TRIANGLE_NUM, the number of triangles.
        //
        //    Input, int TRIANGLE_NODE[3*TRIANGLE_NUM], the nodes that make up each 
        //    triangle.
        //
        //    Output, int TRIANGLE_ORDER3_NEIGHBOR_TRIANGLES[3*TRIANGLE_NUM], 
        //    the three triangles 
        //    that are direct neighbors of a given triangle.  TRIANGLE_NEIGHBOR(1,I) 
        //    is the index of the triangle which touches side 1, defined by nodes 2 
        //    and 3, and so on.  TRIANGLE_NEIGHBOR(1,I) is negative if there is no 
        //    neighbor on that side.  In this case, that side of the triangle lies 
        //    on the boundary of the triangulation.
        //
    {
        int i;
        int j;
        int tri;
        const int triangle_order = 3;

        int[] triangle_neighbor = new int[3 * triangle_num];
        int[] col = new int[4 * 3 * triangle_num];
        //
        //  Step 1.
        //  From the list of nodes for triangle T, of the form: (I,J,K)
        //  construct the three neighbor relations:
        //
        //    (I,J,1,T) or (J,I,1,T),
        //    (J,K,2,T) or (K,J,2,T),
        //    (K,I,3,T) or (I,K,3,T)
        //
        //  where we choose (I,J,1,T) if I < J, or else (J,I,1,T)
        //
        for (tri = 0; tri < triangle_num; tri++)
        {
            i = triangle_node[0 + tri * triangle_order];
            j = triangle_node[1 + tri * triangle_order];
            int k = triangle_node[2 + tri * triangle_order];

            if (i < j)
            {
                col[0 + (3 * tri + 0) * 4] = i;
                col[1 + (3 * tri + 0) * 4] = j;
                col[2 + (3 * tri + 0) * 4] = 1;
                col[3 + (3 * tri + 0) * 4] = tri + 1;
            }
            else
            {
                col[0 + (3 * tri + 0) * 4] = j;
                col[1 + (3 * tri + 0) * 4] = i;
                col[2 + (3 * tri + 0) * 4] = 1;
                col[3 + (3 * tri + 0) * 4] = tri + 1;
            }

            if (j < k)
            {
                col[0 + (3 * tri + 1) * 4] = j;
                col[1 + (3 * tri + 1) * 4] = k;
                col[2 + (3 * tri + 1) * 4] = 2;
                col[3 + (3 * tri + 1) * 4] = tri + 1;
            }
            else
            {
                col[0 + (3 * tri + 1) * 4] = k;
                col[1 + (3 * tri + 1) * 4] = j;
                col[2 + (3 * tri + 1) * 4] = 2;
                col[3 + (3 * tri + 1) * 4] = tri + 1;
            }

            if (k < i)
            {
                col[0 + (3 * tri + 2) * 4] = k;
                col[1 + (3 * tri + 2) * 4] = i;
                col[2 + (3 * tri + 2) * 4] = 3;
                col[3 + (3 * tri + 2) * 4] = tri + 1;
            }
            else
            {
                col[0 + (3 * tri + 2) * 4] = i;
                col[1 + (3 * tri + 2) * 4] = k;
                col[2 + (3 * tri + 2) * 4] = 3;
                col[3 + (3 * tri + 2) * 4] = tri + 1;
            }
        }

        //
        //  Step 2. Perform an ascending dictionary sort on the neighbor relations.
        //  We only intend to sort on rows 1 and 2; the routine we call here
        //  sorts on rows 1 through 4 but that won't hurt us.
        //
        //  What we need is to find cases where two triangles share an edge.
        //  Say they share an edge defined by the nodes I and J.  Then there are
        //  two columns of COL that start out ( I, J, ?, ? ).  By sorting COL,
        //  we make sure that these two columns occur consecutively.  That will
        //  make it easy to notice that the triangles are neighbors.
        //
        typeMethods.i4col_sort_a(4, 3 * triangle_num, ref col);
        //
        //  Step 3. Neighboring triangles show up as consecutive columns with
        //  identical first two entries.  Whenever you spot this happening,
        //  make the appropriate entries in TRIANGLE_NEIGHBOR.
        //
        for (j = 0; j < triangle_num; j++)
        {
            for (i = 0; i < 3; i++)
            {
                triangle_neighbor[i + j * 3] = -1;
            }
        }

        int icol = 1;

        for (;;)
        {
            if (3 * triangle_num <= icol)
            {
                break;
            }

            if (col[0 + (icol - 1) * 4] != col[0 + icol * 4] ||
                col[1 + (icol - 1) * 4] != col[1 + icol * 4])
            {
                icol += 1;
                continue;
            }

            int side1 = col[2 + (icol - 1) * 4];
            int tri1 = col[3 + (icol - 1) * 4];
            int side2 = col[2 + icol * 4];
            int tri2 = col[3 + icol * 4];

            triangle_neighbor[side1 - 1 + (tri1 - 1) * 3] = tri2;
            triangle_neighbor[side2 - 1 + (tri2 - 1) * 3] = tri1;

            icol += 2;
        }

        return triangle_neighbor;
    }

    public static void triangulation_order3_neighbor_triangles_a(int triangle_num,
            int[] triangle_node, ref int[] triangle_neighbor)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULATION_ORDER3_NEIGHBOR_TRIANGLES determines triangle neighbors.
        //
        //  Discussion:
        //
        //    A triangulation of a set of nodes can be completely described by
        //    the coordinates of the nodes, and the list of nodes that make up
        //    each triangle.  However, in some cases, it is necessary to know
        //    triangle adjacency information, that is, which triangle, if any,
        //    is adjacent to a given triangle on a particular side.
        //
        //    This routine creates a data structure recording this information.
        //
        //    The primary amount of work occurs in sorting a list of 3 * TRIANGLE_NUM
        //    data items.
        //
        //  Example:
        //
        //    The input information from TRIANGLE_NODE:
        //
        //    Triangle   Nodes
        //    --------   ---------------
        //     1         3      4      1
        //     2         3      1      2
        //     3         3      2      8
        //     4         2      1      5
        //     5         8      2     13
        //     6         8     13      9
        //     7         3      8      9
        //     8        13      2      5
        //     9         9     13      7
        //    10         7     13      5
        //    11         6      7      5
        //    12         9      7      6
        //    13        10      9      6
        //    14         6      5     12
        //    15        11      6     12
        //    16        10      6     11
        //
        //    The output information in TRIANGLE_NEIGHBOR:
        //
        //    Triangle  Neighboring Triangles
        //    --------  ---------------------
        //
        //     1        -1     -1      2
        //     2         1      4      3
        //     3         2      5      7
        //     4         2     -1      8
        //     5         3      8      6
        //     6         5      9      7
        //     7         3      6     -1
        //     8         5      4     10
        //     9         6     10     12
        //    10         9      8     11
        //    11        12     10     14
        //    12         9     11     13
        //    13        -1     12     16
        //    14        11     -1     15
        //    15        16     14     -1
        //    16        13     15     -1
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 October 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int TRIANGLE_NUM, the number of triangles.
        //
        //    Input, int TRIANGLE_NODE[3*TRIANGLE_NUM], the nodes that make up each 
        //    triangle.
        //
        //    Output, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], the three triangles 
        //    that are direct neighbors of a given triangle.  TRIANGLE_NEIGHBOR(1,I) 
        //    is the index of the triangle which touches side 1, defined by nodes 2 
        //    and 3, and so on.  TRIANGLE_NEIGHBOR(1,I) is negative if there is no 
        //    neighbor on that side.  In this case, that side of the triangle lies 
        //    on the boundary of the triangulation.
        //
    {
        int i;
        int j;
        int tri;
        const int triangle_order = 3;

        triangle_neighbor = new int[triangle_num * 3];

        int[] row = new int[triangle_order * triangle_num * 4];
        //
        //  Step 1.
        //  From the list of nodes for triangle T, of the form: (I,J,K)
        //  construct the three neighbor relations:
        //
        //    (I,J,1,T) or (J,I,1,T),
        //    (J,K,2,T) or (K,J,2,T),
        //    (K,I,3,T) or (I,K,3,T)
        //
        //  where we choose (I,J,1,T) if I < J, or else (J,I,1,T)
        //
        for (tri = 0; tri < triangle_num; tri++)
        {
            i = triangle_node[0 + tri * triangle_order];
            j = triangle_node[1 + tri * triangle_order];
            int k = triangle_node[2 + tri * triangle_order];

            if (i < j)
            {
                row[3 * tri + 0 + 0 * 3 * triangle_num] = i;
                row[3 * tri + 0 + 1 * 3 * triangle_num] = j;
                row[3 * tri + 0 + 2 * 3 * triangle_num] = 1;
                row[3 * tri + 0 + 3 * 3 * triangle_num] = tri + 1;
            }
            else
            {
                row[3 * tri + 0 + 0 * 3 * triangle_num] = j;
                row[3 * tri + 0 + 1 * 3 * triangle_num] = i;
                row[3 * tri + 0 + 2 * 3 * triangle_num] = 1;
                row[3 * tri + 0 + 3 * 3 * triangle_num] = tri + 1;
            }

            if (j < k)
            {
                row[3 * tri + 1 + 0 * 3 * triangle_num] = j;
                row[3 * tri + 1 + 1 * 3 * triangle_num] = k;
                row[3 * tri + 1 + 2 * 3 * triangle_num] = 2;
                row[3 * tri + 1 + 3 * 3 * triangle_num] = tri + 1;
            }
            else
            {
                row[3 * tri + 1 + 0 * 3 * triangle_num] = k;
                row[3 * tri + 1 + 1 * 3 * triangle_num] = j;
                row[3 * tri + 1 + 2 * 3 * triangle_num] = 2;
                row[3 * tri + 1 + 3 * 3 * triangle_num] = tri + 1;
            }

            if (k < i)
            {
                row[3 * tri + 2 + 0 * 3 * triangle_num] = k;
                row[3 * tri + 2 + 1 * 3 * triangle_num] = i;
                row[3 * tri + 2 + 2 * 3 * triangle_num] = 3;
                row[3 * tri + 2 + 3 * 3 * triangle_num] = tri + 1;
            }
            else
            {
                row[3 * tri + 2 + 0 * 3 * triangle_num] = i;
                row[3 * tri + 2 + 1 * 3 * triangle_num] = k;
                row[3 * tri + 2 + 2 * 3 * triangle_num] = 3;
                row[3 * tri + 2 + 3 * 3 * triangle_num] = tri + 1;
            }
        }

        //
        //  Step 2. Perform an ascending dictionary sort on the neighbor relations.
        //  We only intend to sort on columns 1 and 2; the routine we call here
        //  sorts on columns 1 through 4 but that won't hurt us.
        //
        //  What we need is to find cases where two triangles share an edge.
        //  Say they share an edge defined by the nodes I and J.  Then there are
        //  two rows of ROW that start out ( I, J, ?, ? ).  By sorting ROW,
        //  we make sure that these two rows occur consecutively.  That will
        //  make it easy to notice that the triangles are neighbors.
        //
        typeMethods.i4row_sort_a(3 * triangle_num, 4, ref row);
        //
        //  Step 3. Neighboring triangles show up as consecutive rows with
        //  identical first two entries.  Whenever you spot this happening,
        //  make the appropriate entries in TRIANGLE_NEIGHBOR.
        //
        for (j = 0; j < triangle_num; j++)
        {
            for (i = 0; i < 3; i++)
            {
                triangle_neighbor[i + j * 3] = -1;
            }
        }

        int irow = 1;

        for (;;)
        {
            if (3 * triangle_num <= irow)
            {
                break;
            }

            if (row[irow - 1 + 0 * 3 * triangle_num] != row[irow + 0 * 3 * triangle_num] ||
                row[irow - 1 + 1 * 3 * triangle_num] != row[irow + 1 * 3 * triangle_num])
            {
                irow += 1;
                continue;
            }

            int side1 = row[irow - 1 + 2 * 3 * triangle_num];
            int tri1 = row[irow - 1 + 3 * 3 * triangle_num];
            int side2 = row[irow + 2 * 3 * triangle_num];
            int tri2 = row[irow + 3 * 3 * triangle_num];

            triangle_neighbor[side1 - 1 + (tri1 - 1) * 3] = tri2;
            triangle_neighbor[side2 - 1 + (tri2 - 1) * 3] = tri1;

            irow += 2;
        }
    }
}