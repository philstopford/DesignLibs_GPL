﻿using Burkardt.Types;

namespace Burkardt.TriangulationNS
{
    public static partial class Neighbor
    {
        public static void triangulation_order6_neighbor_triangles(int triangle_num,
            int[] triangle_node, ref int[] triangle_neighbor )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULATION_ORDER6_NEIGHBOR_TRIANGLES determines triangle neighbors.
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
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 June 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int TRIANGLE_NUM, the number of triangles.
        //
        //    Input, int TRIANGLE_NODE[6*TRIANGLE_NUM], the nodes that make up each triangle.
        //
        //    Output, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], the three triangles that are direct
        //    neighbors of a given triangle.  TRIANGLE_NEIGHBOR(1,I) is the index of the triangle
        //    which touches side 1, defined by nodes 2 and 3, and so on.  TRIANGLE_NEIGHBOR(1,I)
        //    is negative if there is no neighbor on that side.  In this case, that
        //    side of the triangle lies on the boundary of the triangulation.
        //
        {
            int i;
            int irow;
            int j;
            int k;
            int[] row;
            int side1;
            int side2;
            int tri;
            int tri1;
            int tri2;

            row = new int[3 * triangle_num * 4];
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
                i = triangle_node[0 + tri * 6];
                j = triangle_node[1 + tri * 6];
                k = triangle_node[2 + tri * 6];

                if (i < j)
                {
                    row[3 * tri + 0 * 3 * triangle_num] = i;
                    row[3 * tri + 1 * 3 * triangle_num] = j;
                    row[3 * tri + 2 * 3 * triangle_num] = 1;
                    row[3 * tri + 3 * 3 * triangle_num] = tri + 1;
                }
                else
                {
                    row[3 * tri + 0 * 3 * triangle_num] = j;
                    row[3 * tri + 1 * 3 * triangle_num] = i;
                    row[3 * tri + 2 * 3 * triangle_num] = 1;
                    row[3 * tri + 3 * 3 * triangle_num] = tri + 1;
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

            irow = 1;

            for (;;)
            {
                if (3 * triangle_num <= irow)
                {
                    break;
                }

                if (row[irow - 1 + 0 * 3 * triangle_num] != row[irow + 0 * 3 * triangle_num] ||
                    row[irow - 1 + 1 * 3 * triangle_num] != row[irow + 1 * 3 * triangle_num])
                {
                    irow = irow + 1;
                    continue;
                }

                side1 = row[irow - 1 + 2 * 3 * triangle_num];
                tri1 = row[irow - 1 + 3 * 3 * triangle_num];
                side2 = row[irow + 2 * 3 * triangle_num];
                tri2 = row[irow + 3 * 3 * triangle_num];

                triangle_neighbor[side1 - 1 + (tri1 - 1) * 3] = tri2;
                triangle_neighbor[side2 - 1 + (tri2 - 1) * 3] = tri1;

                irow = irow + 2;
            }
        }
    }
}