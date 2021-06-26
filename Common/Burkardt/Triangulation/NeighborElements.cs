using Burkardt.Types;

namespace Burkardt.TriangulationNS
{
    public static class NeighborElements
    {
        public static int[] triangulation_neighbor_elements(int triangle_order, int triangle_num,
                int[] triangle_node)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGULATION_NEIGHBOR_ELEMENTS determines element neighbors.
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
            //    07 September 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int TRIANGLE_ORDER, the order of the triangles.
            //
            //    Input, int TRIANGLE_NUM, the number of triangles.
            //
            //    Input, int TRIANGLE_NODE[TRIANGLE_ORDER*TRIANGLE_NUM], the nodes that
            //    make up each triangle.
            //
            //    Output, int TRIANGLE_NEIGHBOR_TRIANGLES[3*TRIANGLE_NUM],
            //    the three triangles
            //    that are direct neighbors of a given triangle.  TRIANGLE_NEIGHBOR(1,I)
            //    is the index of the triangle which touches side 1, defined by nodes 2
            //    and 3, and so on.  TRIANGLE_NEIGHBOR(1,I) is negative if there is no
            //    neighbor on that side.  In this case, that side of the triangle lies
            //    on the boundary of the triangulation.
            //
        {
            int[] col;
            int i;
            int icol;
            int j;
            int k;
            int side1;
            int side2;
            int tri;
            int tri1;
            int tri2;
            int[] triangle_neighbor;

            triangle_neighbor = new int[3 * triangle_num];
            col = new int[4 * (3 * triangle_num)];
            //
            //  Step 1.
            //  From the list of nodes for triangle T, of the form: (I,J,K)
            //  construct the three neighbor relations:
            //
            //    (I,J,3,T) or (J,I,3,T),
            //    (J,K,1,T) or (K,J,1,T),
            //    (K,I,2,T) or (I,K,2,T)
            //
            //  where we choose (I,J,3,T) if I < J, or else (J,I,3,T)
            //
            for (tri = 0; tri < triangle_num; tri++)
            {
                i = triangle_node[0 + tri * triangle_order];
                j = triangle_node[1 + tri * triangle_order];
                k = triangle_node[2 + tri * triangle_order];

                if (i < j)
                {
                    col[0 + (3 * tri + 0) * 4] = i;
                    col[1 + (3 * tri + 0) * 4] = j;
                    col[2 + (3 * tri + 0) * 4] = 3;
                    col[3 + (3 * tri + 0) * 4] = tri + 1;
                }
                else
                {
                    col[0 + (3 * tri + 0) * 4] = j;
                    col[1 + (3 * tri + 0) * 4] = i;
                    col[2 + (3 * tri + 0) * 4] = 3;
                    col[3 + (3 * tri + 0) * 4] = tri + 1;
                }

                if (j < k)
                {
                    col[0 + (3 * tri + 1) * 4] = j;
                    col[1 + (3 * tri + 1) * 4] = k;
                    col[2 + (3 * tri + 1) * 4] = 1;
                    col[3 + (3 * tri + 1) * 4] = tri + 1;
                }
                else
                {
                    col[0 + (3 * tri + 1) * 4] = k;
                    col[1 + (3 * tri + 1) * 4] = j;
                    col[2 + (3 * tri + 1) * 4] = 1;
                    col[3 + (3 * tri + 1) * 4] = tri + 1;
                }

                if (k < i)
                {
                    col[0 + (3 * tri + 2) * 4] = k;
                    col[1 + (3 * tri + 2) * 4] = i;
                    col[2 + (3 * tri + 2) * 4] = 2;
                    col[3 + (3 * tri + 2) * 4] = tri + 1;
                }
                else
                {
                    col[0 + (3 * tri + 2) * 4] = i;
                    col[1 + (3 * tri + 2) * 4] = k;
                    col[2 + (3 * tri + 2) * 4] = 2;
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

            icol = 1;

            for (;;)
            {
                if (3 * triangle_num <= icol)
                {
                    break;
                }

                if (col[0 + (icol - 1) * 4] != col[0 + icol * 4] ||
                    col[1 + (icol - 1) * 4] != col[1 + icol * 4])
                {
                    icol = icol + 1;
                    continue;
                }

                side1 = col[2 + (icol - 1) * 4];
                tri1 = col[3 + (icol - 1) * 4];
                side2 = col[2 + icol * 4];
                tri2 = col[3 + icol * 4];

                triangle_neighbor[side1 - 1 + (tri1 - 1) * 3] = tri2;
                triangle_neighbor[side2 - 1 + (tri2 - 1) * 3] = tri1;

                icol = icol + 2;
            }
            
            return triangle_neighbor;
        }
    }
}