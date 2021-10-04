using System;
using Burkardt.Types;

namespace Burkardt.QuadMesh
{
    public class Neighbors
    {
        public static int[] neighbor_elements_q4_mesh(int element_num, int[] element_node)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NEIGHBOR_ELEMENTS_Q4_MESH determines element neighbors in a Q4 mesh.
            //
            //  Discussion:
            //
            //    A quadrilateral mesh of a set of nodes can be completely described by
            //    the coordinates of the nodes, and the list of nodes that make up
            //    each element.  However, in some cases, it is necessary to know
            //    element adjacency information, that is, which element, if any,
            //    is adjacent to a given element on a particular side.
            //
            //    This routine creates a data structure recording this information.
            //
            //    The primary amount of work occurs in sorting a list of 4 * ELEMENT_NUM
            //    data items.
            //
            //    Note that COL is a work array allocated dynamically inside this
            //    routine.  It is possible, for very large values of ELEMENT_NUM,
            //    that the necessary amount of memory will not be accessible, and the
            //    routine will fail.  This is a limitation of the implementation of
            //    dynamic arrays in FORTRAN90.  One way to get around this would be
            //    to require the user to declare ROW in the calling routine
            //    as an allocatable array, get the necessary memory explicitly with
            //    an ALLOCATE statement, and then pass ROW into this routine.
            //
            //    Of course, the point of dynamic arrays was to make it easy to
            //    hide these sorts of temporary work arrays from the poor user!
            //
            //    This routine was revised to store the edge data in a column
            //    array rather than a row array.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 September 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int ELEMENT_NUM, the number of elements.
            //
            //    Input, int ELEMENT_NODE[4*ELEMENT_NUM], the nodes that make up each element.
            //
            //    Output, int NEIGHBOR_ELEMENTS_Q4_MESH[4*ELEMENT_NUM], the elements that are direct 
            //    neighbors of a given element, or -1 if there is no neighbor on that side.
            //
        {
            int[] col;
            int element;
            int[] element_neighbor;
            int element_order = 4;
            int element1;
            int element2;
            int i;
            int icol;
            int j;
            int k;
            int l;
            int side1;
            int side2;

            element_neighbor = new int[4 * element_num];
            col = new int[4 * 4 * element_num];
            //
            //  Step 1.
            //  From the list of nodes for element E, of the form: (I,J,K,L)
            //  construct the four neighbor relations:
            //
            //    (I,J,0,E) or (J,I,0,E),
            //    (J,K,1,E) or (K,J,1,E),
            //    (K,L,2,E) or (L,K,2,E)
            //    (L,I,3,E) or (I,L,3,E)
            //
            //  where we choose (I,J,0,E) if I < J, or else (J,I,0,E)
            //
            for (element = 0; element < element_num; element++)
            {
                i = element_node[0 + element * element_order];
                j = element_node[1 + element * element_order];
                k = element_node[2 + element * element_order];
                l = element_node[3 + element * element_order];

                col[0 + 0 * 4 + 16 * element] = Math.Min(i, j);
                col[1 + 0 * 4 + 16 * element] = Math.Max(i, j);
                col[2 + 0 * 4 + 16 * element] = 0;
                col[3 + 0 * 4 + 16 * element] = element;

                col[0 + 1 * 4 + 16 * element] = Math.Min(j, k);
                col[1 + 1 * 4 + 16 * element] = Math.Max(j, k);
                col[2 + 1 * 4 + 16 * element] = 1;
                col[3 + 1 * 4 + 16 * element] = element;

                col[0 + 2 * 4 + 16 * element] = Math.Min(k, l);
                col[1 + 2 * 4 + 16 * element] = Math.Max(k, l);
                col[2 + 2 * 4 + 16 * element] = 2;
                col[3 + 2 * 4 + 16 * element] = element;

                col[0 + 3 * 4 + 16 * element] = Math.Min(l, i);
                col[1 + 3 * 4 + 16 * element] = Math.Max(l, i);
                col[2 + 3 * 4 + 16 * element] = 3;
                col[3 + 3 * 4 + 16 * element] = element;
            }

            //
            //  Step 2. Perform an ascending dictionary sort on the neighbor relations.
            //  We only intend to sort on rows 1 and 2; the routine we call here
            //  sorts on rows 1 through 4 but that won't hurt us.
            //
            //  What we need is to find cases where two elements share an edge.
            //  Say they share an edge defined by the nodes I and J.  Then there are
            // two columns of COL that start out ( I, J, ?, ? ).  By sorting COL,
            //  we make sure that these two columns occur consecutively.  That will
            //  make it easy to notice that the elements are neighbors.
            //
            typeMethods.i4col_sort_a(4, 4 * element_num, ref col);
            //
            //  Step 3. Neighboring elements show up as consecutive columns with
            //  identical first two entries.  Whenever you spot this happening,
            //  make the appropriate entries in ELEMENT_NEIGHBOR.
            //
            for (j = 0; j < element_num; j++)
            {
                for (i = 0; i < 4; i++)
                {
                    element_neighbor[i + j * 4] = -1;
                }
            }

            icol = 0;

            for (;;)
            {
                if (4 * element_num <= icol)
                {
                    break;
                }

                if (col[(0 + icol * 4) % col.Length] != col[(0 + (icol + 1) * 4) % col.Length] ||
                    col[(1 + icol * 4) % col.Length] != col[(1 + (icol + 1) * 4) % col.Length])
                {
                    icol = icol + 1;
                    continue;
                }

                side1 = col[(2 + icol * 4) % col.Length];
                element1 = col[(3 + icol * 4) % col.Length];

                side2 = col[(2 + (icol + 1) * 4) % col.Length];
                element2 = col[(3 + (icol + 1) * 4) % col.Length];

                element_neighbor[side1 + element1 * 4] = element2;
                element_neighbor[side2 + element2 * 4] = element1;

                icol = icol + 2;
            }

            return element_neighbor;
        }
    }
}