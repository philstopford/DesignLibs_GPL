using Burkardt.Types;

namespace Burkardt.QuadMesh;

public static class Adjacency
{
    public static int[] adj_set_q4_mesh(int node_num, int element_num,
            int[] element_node, int[] element_neighbor, int adj_num, int[] adj_row)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ADJ_SET_Q4_MESH sets adjacencies in a triangulation.
        //
        //  Discussion:
        //
        //    This routine is called to set the adjacencies, after the
        //    appropriate amount of memory has been set aside for storage.
        //
        //    The mesh is assumed to involve 4-node quadrilaterals.
        //
        //    Two nodes are "adjacent" if they are both nodes in some element.
        //    Also, a node is considered to be adjacent to itself.
        //
        //    This routine can be used to create the compressed column storage
        //    for a linear element finite element discretization of
        //    Poisson's equation in two dimensions.
        //
        //  Diagram:
        //
        //         side 3
        //       4-------3
        //    s  |       |  s
        //    i  |       |  i
        //    d  |       |  d
        //    e  |       |  e
        //       |       |
        //    4  |       |  2
        //       |       |
        //       1-------2
        //
        //         side 1
        //
        //    The local node numbering
        //
        //
        //   20-21-22-23-24
        //    |  |  |  |  |
        //    |  |  |  |  |
        //   15-16-17-18-19
        //    |  |  |  |  |
        //    |  |  |  |  |
        //   10-11-12-13-14
        //    |  |  |  |  |
        //    |  |  |  |  |
        //    5--6--7--8--9
        //    |  |  |  |  |
        //    |  |  |  |  |
        //    0--1--2--3--4
        //
        //    A sample grid.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 September 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, int ELEMENT_NUM, the number of elements.
        //
        //    Input, int ELEMENT_NODE[4*ELEMENT_NUM], lists the nodes that
        //    make up each element in counterclockwise order.
        //
        //    Input, int ELEMENT_NEIGHBOR[4*ELEMENT_NUM], for each side of
        //    an element, lists the neighboring element, or -1 if there is
        //    no neighbor.
        //
        //    Input, int ADJ_NUM, the number of adjacencies.
        //
        //    Input, int ADJ_ROW[NODE_NUM+1].  Information about column J is stored
        //    in entries ADJ_ROW(J) through ADJ_ROW(J+1)-1 of ADJ.
        //
        //    Output, int ADJ_SET_Q4_MESH[ADJ_NUM], the adjacency information.
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
        int n4;
        int node;
        int element;
        int element2;
        int element_order = 4;

        adj = new int[adj_num];
        for (k = 0; k < adj_num; k++)
        {
            adj[k] = -1;
        }

        adj_copy = new int[node_num];
        for (node = 0; node < node_num; node++)
        {
            adj_copy[node] = adj_row[node];
        }

        //
        //  Set every node to be adjacent to itself.
        //
        for (node = 0; node < node_num; node++)
        {
            adj[adj_copy[node]] = node;
            adj_copy[node] += 1;
        }

        //
        //  Examine each element.
        //
        for (element = 0; element < element_num; element++)
        {
            n1 = element_node[0 + element * element_order];
            n2 = element_node[1 + element * element_order];
            n3 = element_node[2 + element * element_order];
            n4 = element_node[3 + element * element_order];
            //
            //  Add edges (1,3) and (2,4).  There is no need to check for redundancy,
            //  since this is the only case when these nodes can share an element.
            //
            adj[adj_copy[n1]] = n3;
            adj_copy[n1] += 1;
            adj[adj_copy[n3]] = n1;
            adj_copy[n3] += 1;

            adj[adj_copy[n2]] = n4;
            adj_copy[n2] += 1;
            adj[adj_copy[n4]] = n2;
            adj_copy[n4] += 1;
            //
            //  Add edge (1,2) if this is the first occurrence,
            //  that is, if the edge (1,2) is on a boundary (ELEMENT2 <= 0)
            //  or if this element is the first of the pair in which the edge
            //  occurs (ELEMENT < ELEMENT2).
            //
            element2 = element_neighbor[0 + element * 4];

            if (element2 < 0 || element < element2)
            {
                adj[adj_copy[n1]] = n2;
                adj_copy[n1] += 1;
                adj[adj_copy[n2]] = n1;
                adj_copy[n2] += 1;
            }

            //
            //  Add edge (2,3).
            //
            element2 = element_neighbor[1 + element * 4];

            if (element2 < 0 || element < element2)
            {
                adj[adj_copy[n2]] = n3;
                adj_copy[n2] += 1;
                adj[adj_copy[n3]] = n2;
                adj_copy[n3] += 1;
            }

            //
            //  Add edge (3,4).
            //
            element2 = element_neighbor[2 + element * 4];

            if (element2 < 0 || element < element2)
            {
                adj[adj_copy[n4]] = n3;
                adj_copy[n4] += 1;
                adj[adj_copy[n3]] = n4;
                adj_copy[n3] += 1;
            }

            //
            //  Add edge (4,1).
            //
            element2 = element_neighbor[3 + element * 4];

            if (element2 < 0 || element < element2)
            {
                adj[adj_copy[n1]] = n4;
                adj_copy[n1] += 1;
                adj[adj_copy[n4]] = n1;
                adj_copy[n4] += 1;
            }
        }

        //
        //  Ascending sort the entries for each node.
        //
        for (node = 0; node < node_num; node++)
        {
            k1 = adj_row[node];
            k2 = adj_row[node + 1] - 1;
            typeMethods.i4vec_sort_heap_a(k2 + 1 - k1, ref adj, aIndex: +k1);
        }

        return adj;
    }

    public static int adj_size_q4_mesh(int node_num, int element_num, int[] element_node,
            int[] element_neighbor, ref int[] adj_row)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ADJ_SIZE_Q4_MESH counts adjacencies in a Q4 mesh.
        //
        //  Discussion:
        //
        //    This routine is called to count the adjacencies, so that the
        //    appropriate amount of memory can be set aside for storage when
        //    the adjacency structure is created.
        //
        //    The mesh is assumed to involve 4-node quadrilaterals.
        //
        //    Two nodes are "adjacent" if they are both nodes in some quadrilateral.
        //    Also, a node is considered to be adjacent to itself.
        //
        //  Diagram:
        //
        //         side 3
        //       4-------3
        //    s  |       |  s
        //    i  |       |  i
        //    d  |       |  d
        //    e  |       |  e
        //       |       |
        //    4  |       |  2
        //       |       |
        //       1-------2
        //
        //         side 1
        //
        //    The local node numbering
        //
        //
        //   20-21-22-23-24
        //    |  |  |  |  |
        //    |  |  |  |  |
        //   15-16-17-18-19
        //    |  |  |  |  |
        //    |  |  |  |  |
        //   10-11-12-13-14
        //    |  |  |  |  |
        //    |  |  |  |  |
        //    5--6--7--8--9
        //    |  |  |  |  |
        //    |  |  |  |  |
        //    0--1--2--3--4
        //
        //    A sample grid.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 September 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, int ELEMENT_NUM, the number of elements.
        //
        //    Input, int ELEMENT_NODE[4*ELEMENT_NUM], lists the nodes that
        //    make up each element, in counterclockwise order.
        //
        //    Input, int ELEMENT_NEIGHBOR[4*ELEMENT_NUM], for each side of
        //    a element, lists the neighboring elment, or -1 if there is
        //    no neighbor.
        //
        //    Output, int ADJ_ROW[NODE_NUM+1], Information about column J is stored
        //    in entries ADJ_ROW[J] through ADJ_ROW[J+1]-1 of ADJ.
        //
        //    Output, int ADJ_SIZE_Q4_MESH, the number of adjacencies.
        //
    {
        int adj_num;
        int element;
        int element_order = 4;
        int element2;
        int i;
        int n1;
        int n2;
        int n3;
        int n4;
        int node;

        adj_num = 0;
        //
        //  Set every node to be adjacent to itself.
        //
        for (node = 0; node < node_num; node++)
        {
            adj_row[node] = 1;
        }

        //
        //  Examine each element.
        //
        for (element = 0; element < element_num; element++)
        {
            n1 = element_node[0 + element * element_order];
            n2 = element_node[1 + element * element_order];
            n3 = element_node[2 + element * element_order];
            n4 = element_node[3 + element * element_order];
            //
            //  Add edge (1,3).
            //
            adj_row[n1] += 1;
            adj_row[n3] += 1;
            //
            //  Add edge (2,4).
            //
            adj_row[n2] += 1;
            adj_row[n4] += 1;
            //
            //  Add edge (1,2) if this is the first occurrence,
            //  that is, if the edge (1,2) is on a boundary (ELEMENT2 <= 0)
            //  or if this element is the first of the pair in which the edge
            //  occurs (ELEMENT < ELEMENT2).
            //
            element2 = element_neighbor[0 + element * 4];

            if (element2 < 0 || element < element2)
            {
                adj_row[n1] += 1;
                adj_row[n2] += 1;
            }

            //
            //  Add edge (2,3).
            //
            element2 = element_neighbor[1 + element * 4];

            if (element2 < 0 || element < element2)
            {
                adj_row[n2] += 1;
                adj_row[n3] += 1;
            }

            //
            //  Add edge (3,4).
            //
            element2 = element_neighbor[2 + element * 4];

            if (element2 < 0 || element < element2)
            {
                adj_row[n3] += 1;
                adj_row[n4] += 1;
            }

            //
            //  Add edge (4,1).
            //
            element2 = element_neighbor[3 + element * 4];

            if (element2 < 0 || element < element2)
            {
                adj_row[n4] += 1;
                adj_row[n1] += 1;
            }
        }

        //
        //  We used ADJ_ROW to count the number of entries in each column.
        //  Convert it to pointers into the ADJ array.
        //
        for (node = node_num; 1 <= node; node--)
        {
            adj_row[node] = adj_row[node - 1];
        }

        adj_row[0] = 0;
        for (i = 1; i <= node_num; i++)
        {
            adj_row[i] += adj_row[i - 1];
        }

        //
        //  Finally, record the total number of adjacencies.
        //
        adj_num = adj_row[node_num];

        return adj_num;
    }
}