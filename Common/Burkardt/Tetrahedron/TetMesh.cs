using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.TetrahedronNS
{
    public static class TetMesh
    {
        public static void tet_mesh_base_one(int node_num, int element_order, int element_num,
                ref int[] element_node)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TET_MESH_BASE_ONE ensures that the element definition is 1-based.
            //
            //  Discussion:
            //
            //    The ELEMENT_NODE array contains nodes indices that form elements.
            //    The convention for node indexing might start at 0 or at 1.
            //
            //    This function attempts to detect 0-based node indexing and correct it.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    27 September 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NODE_NUM, the number of nodes.
            //
            //    Input, int ELEMENT_ORDER, the order of the elements.
            //
            //    Input, int ELEMENT_NUM, the number of elements.
            //
            //    Input/output, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], the element
            //    definitions.
            //
        {
            int element;
            int node_max;
            int node_min;
            int order;

            node_min = typeMethods.i4mat_min(element_order, element_num, element_node);
            node_max = typeMethods.i4mat_max(element_order, element_num, element_node);

            if (node_min == 1 && node_max == node_num)
            {
                Console.WriteLine("");
                Console.WriteLine("TET_MESH_BASE_ONE:");
                Console.WriteLine("  The element indexing appears to be 1-based!");
                Console.WriteLine("  No conversion is necessary.");
            }
            else if (node_min == 0 && node_max == node_num - 1)
            {
                Console.WriteLine("");
                Console.WriteLine("TET_MESH_BASE_ONE:");
                Console.WriteLine("  The element indexing appears to be 0-based!");
                Console.WriteLine("  This will be converted to 1-based.");
                for (element = 0; element < element_num; element++)
                {
                    for (order = 0; order < element_order; order++)
                    {
                        element_node[order + element * element_order] =
                            element_node[order + element * element_order] + 1;
                    }
                }
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("TET_MESH_BASE_ONE - Warning!");
                Console.WriteLine("  The element indexing is not of a recognized type.");
            }
        }

        public static int tet_mesh_base_zero(int node_num, int element_order, int element_num,
                ref int[] element_node)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TET_MESH_BASE_ZERO ensures that the element definition is zero-based.
            //
            //  Discussion:
            //
            //    The ELEMENT_NODE array contains nodes indices that form elements.
            //    The convention for node indexing might start at 0 or at 1.
            //    Since a C++ program will naturally assume a 0-based indexing, it is
            //    necessary to check a given element definition and, if it is actually
            //    1-based, to convert it.
            //
            //    This function attempts to detect 1-based node indexing and correct it.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    27 September 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NODE_NUM, the number of nodes.
            //
            //    Input, int ELEMENT_ORDER, the order of the elements.
            //
            //    Input, int ELEMENT_NUM, the number of elements.
            //
            //    Input/output, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], the element
            //    definitions.
            //
        {
            int base_;
            int element;
            int node;
            int node_max;
            int node_min;
            int order;
            //
            //  If the element information is 1-based, make it 0-based.
            //
            node_min = node_num + 1;
            node_max = -1;
            for (element = 0; element < element_num; element++)
            {
                for (order = 0; order < element_order; order++)
                {
                    node = element_node[order + element * element_order];
                    node_min = Math.Min(node_min, node);
                    node_max = Math.Max(node_max, node);
                }
            }

            if (node_min == 1 && node_max == node_num)
            {
                base_ = 1;
                Console.WriteLine("");
                Console.WriteLine("TET_MESH_BASE_ZERO:");
                Console.WriteLine("  The element indexing appears to be 1-based!");
                Console.WriteLine("  This will be converted to 0-based.");
                for (element = 0; element < element_num; element++)
                {
                    for (order = 0; order < element_order; order++)
                    {
                        element_node[order + element * element_order] =
                            element_node[order + element * element_order] - 1;
                    }
                }
            }
            else if (node_min == 0 && node_max == node_num - 1)
            {
                base_ = 0;
                Console.WriteLine("");
                Console.WriteLine("TET_MESH_BASE_ZERO:");
                Console.WriteLine("  The element indexing appears to be 0-based!");
                Console.WriteLine("  No conversion is necessary.");
            }
            else
            {
                base_ = -1;
                Console.WriteLine("");
                Console.WriteLine("TET_MESH_BASE_ZERO - Warning!");
                Console.WriteLine("  The element indexing is not of a recognized type.");
            }

            return base_;
        }

        public static int[] tet_mesh_neighbor_tets(int tetra_order, int tetra_num,
                int[] tetra_node)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TET_MESH_NEIGHBOR_TETS determines tetrahedron neighbors.
            //
            //  Discussion:
            //
            //    A tet mesh of a set of nodes can be completely described by
            //    the coordinates of the nodes, and the list of nodes that make up
            //    each tetrahedron.  In the most common case, four nodes are used.
            //    There is also a 10 node case, where nodes are also placed on
            //    the midsides of the tetrahedral edges.
            //
            //    This routine can handle 4 or 10-node tetrahedral meshes.  The
            //    10-node case is handled simply by ignoring the six midside nodes,
            //    which are presumed to be listed after the vertices.
            //
            //    The tetrahedron adjacency information records which tetrahedron
            //    is adjacent to a given tetrahedron on a particular face.
            //
            //    This routine creates a data structure recording this information.
            //
            //    The primary amount of work occurs in sorting a list of 4 * TETRA_NUM
            //    data items.
            //
            //    The neighbor tetrahedrons are indexed by the face they share with
            //    the tetrahedron.
            //
            //    Each face of the tetrahedron is indexed by the node which is NOT
            //    part of the face.  That is:
            //
            //    * Neighbor 1 shares face 1 defined by nodes 2, 3, 4.
            //    * Neighbor 2 shares face 2 defined by nodes 1, 3, 4;
            //    * Neighbor 3 shares face 3 defined by nodes 1, 2, 4;
            //    * Neighbor 4 shares face 4 defined by nodes 1, 2, 3.
            //
            //    For instance, if the (transposed) TETRA_NODE array was:
            //
            //    Row       1      2      3      4
            //    Col
            //
            //      1       4      3      5      1
            //      2       4      2      5      1
            //      3       4      7      3      5
            //      4       4      7      8      5
            //      5       4      6      2      5
            //      6       4      6      8      5
            //
            //    then the (transposed) TETRA_NEIGHBOR array should be:
            //
            //    Row       1      2      3      4
            //    Col
            //
            //      1      -1      2     -1      3
            //      2      -1      1     -1      5
            //      3      -1      1      4     -1
            //      4      -1      6      3     -1
            //      5      -1      2      6     -1
            //      6      -1      4      5     -1
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int TETRA_ORDER, the order of the tetrahedrons.
            //
            //    Input, int TETRA_NUM, the number of tetrahedrons.
            //
            //    Input, int TETRA_NODE[TETRA_ORDER*TETRA_NUM], the indices of the nodes.
            //
            //    Output, int TET_MESH_NEIGHBORS[4*TETRA_NUM], the four tetrahedrons that
            //    are direct neighbors of a given tetrahedron.  If there is no neighbor
            //    sharing a given face, the index is set to -1.
            //
        {
            int a = 0;
            int b = 0;
            int c = 0;
            int face;
            int face1;
            int face2;
            int[] faces;
            int i;
            int j;
            int k;
            int l;
            int tetra;
            int[] tetra_neighbor;
            int tetra1;
            int tetra2;

            faces = new int[5 * (4 * tetra_num)];
            tetra_neighbor = new int[4 * tetra_num];
            //
            //  Step 1.
            //  From the list of nodes for tetrahedron T, of the form: (I,J,K,L)
            //  construct the four face relations:
            //
            //    (J,K,L,1,T)
            //    (I,K,L,2,T)
            //    (I,J,L,3,T)
            //    (I,J,K,4,T)
            //
            //  In order to make matching easier, we reorder each triple of nodes
            //  into ascending order.
            //
            for (tetra = 0; tetra < tetra_num; tetra++)
            {
                i = tetra_node[0 + tetra * tetra_order];
                j = tetra_node[1 + tetra * tetra_order];
                k = tetra_node[2 + tetra * tetra_order];
                l = tetra_node[3 + tetra * tetra_order];

                typeMethods.i4i4i4_sort_a(j, k, l, ref a, ref b, ref c);

                faces[0 + 0 * 5 + tetra * 5 * 4] = a;
                faces[1 + 0 * 5 + tetra * 5 * 4] = b;
                faces[2 + 0 * 5 + tetra * 5 * 4] = c;
                faces[3 + 0 * 5 + tetra * 5 * 4] = 0;
                faces[4 + 0 * 5 + tetra * 5 * 4] = tetra;

                typeMethods.i4i4i4_sort_a(i, k, l, ref a, ref b, ref c);

                faces[0 + 1 * 5 + tetra * 5 * 4] = a;
                faces[1 + 1 * 5 + tetra * 5 * 4] = b;
                faces[2 + 1 * 5 + tetra * 5 * 4] = c;
                faces[3 + 1 * 5 + tetra * 5 * 4] = 1;
                faces[4 + 1 * 5 + tetra * 5 * 4] = tetra;

                typeMethods.i4i4i4_sort_a(i, j, l, ref a, ref b, ref c);

                faces[0 + 2 * 5 + tetra * 5 * 4] = a;
                faces[1 + 2 * 5 + tetra * 5 * 4] = b;
                faces[2 + 2 * 5 + tetra * 5 * 4] = c;
                faces[3 + 2 * 5 + tetra * 5 * 4] = 2;
                faces[4 + 2 * 5 + tetra * 5 * 4] = tetra;

                typeMethods.i4i4i4_sort_a(i, j, k, ref a, ref b, ref c);

                faces[0 + 3 * 5 + tetra * 5 * 4] = a;
                faces[1 + 3 * 5 + tetra * 5 * 4] = b;
                faces[2 + 3 * 5 + tetra * 5 * 4] = c;
                faces[3 + 3 * 5 + tetra * 5 * 4] = 3;
                faces[4 + 3 * 5 + tetra * 5 * 4] = tetra;
            }

            //
            //  Step 2. Perform an ascending dictionary sort on the neighbor relations.
            //  We only intend to sort on rows 1:3; the routine we call here
            //  sorts on rows 1 through 5 but that won't hurt us.
            //
            //  What we need is to find cases where two tetrahedrons share a face.
            //  By sorting the columns of the FACES array, we will put shared faces
            //  next to each other.
            //
            typeMethods.i4col_sort_a(5, 4 * tetra_num, ref faces);
            //
            //  Step 3. Neighboring tetrahedrons show up as consecutive columns with
            //  identical first three entries.  Whenever you spot this happening,
            //  make the appropriate entries in TETRA_NEIGHBOR.
            //
            for (j = 0; j < tetra_num; j++)
            {
                for (i = 0; i < 4; i++)
                {
                    tetra_neighbor[i + j * 4] = -1;
                }
            }

            face = 0;

            for (;;)
            {
                if (4 * tetra_num - 1 <= face)
                {
                    break;
                }

                if (faces[0 + face * 5] == faces[0 + (face + 1) * 5] &&
                    faces[1 + face * 5] == faces[1 + (face + 1) * 5] &&
                    faces[2 + face * 5] == faces[2 + (face + 1) * 5])
                {
                    face1 = faces[3 + face * 5];
                    tetra1 = faces[4 + face * 5];
                    face2 = faces[3 + (face + 1) * 5];
                    tetra2 = faces[4 + (face + 1) * 5];
                    tetra_neighbor[face1 + tetra1 * 4] = tetra2;
                    tetra_neighbor[face2 + tetra2 * 4] = tetra1;
                    face = face + 2;
                }
                else
                {
                    face = face + 1;
                }
            }

            return tetra_neighbor;
        }

        public static int[] tet_mesh_node_order(int tetra_order, int tetra_num, int[] tetra_node,
                int node_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TET_MESH_NODE_ORDER: determines the order of nodes.
            //
            //  Discussion:
            //
            //    The order of a node is the number of tetrahedrons that use that node
            //    as a vertex.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    27 October 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int TETRA_ORDER, the order of the tetrahedrons.
            //
            //    Input, int TETRA_NUM, the number of tetrahedrons.
            //
            //    Input, int TETRA_NODE[TETRA_ORDER*TETRA_NUM], the indices of the nodes.
            //
            //    Input, int NODE_NUM, the number of nodes.
            //
            //    Output, int TET_MESH_NODE_ORDER[NODE_NUM], the order of each node.
            //
        {
            int i;
            int node;
            int[] node_order;
            int tetra;

            node_order = new int[node_num];

            typeMethods.i4vec_zero(node_num, ref node_order);

            for (tetra = 0; tetra < tetra_num; tetra++)
            {
                for (i = 0; i < tetra_order; i++)
                {
                    node = tetra_node[i + tetra * tetra_order];
                    if (node < 0 || node_num <= node)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("TET_MESH_NODE_ORDER - Fatal error!");
                        Console.WriteLine("  Illegal entry in TETRA_NODE.");
                        return null;
                    }
                    else
                    {
                        node_order[node] = node_order[node] + 1;
                    }
                }
            }

            return node_order;
        }

        public static void tet_mesh_order4_adj_count(int node_num, int tetra_num,
                int[] tetra_node, ref int adj_num, ref int[] adj_row)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TET_MESH_ORDER4_ADJ_COUNT counts the number of nodal adjacencies.
            //
            //  Discussion:
            //
            //    Assuming that the tet mesh is to be used in a finite element
            //    computation, we declare that two distinct nodes are "adjacent" if and
            //    only if they are both included in some tetrahedron.
            //
            //    It is the purpose of this routine to determine the number of
            //    such adjacency relationships.
            //
            //    The initial count gets only the (I,J) relationships, for which
            //    node I is strictly less than node J.  This value is doubled
            //    to account for symmetry.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    13 January 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NODE_NUM, the number of nodes.
            //
            //    Input, int TETRA_NUM, the number of tetrahedrons.
            //
            //    Input, int TETRA_NODE[4*TETRA_NUM], the indices of the nodes.
            //
            //    Output, int *ADJ_NUM, the total number of adjacency relationships,
            //
            //    Output, int ADJ_ROW[NODE_NUM+1], the ADJ pointer array.
            //
        {
            int i = 0;
            int j;
            int k;
            int node;
            int[] pair;
            int pair_num;
            int pair_unique_num;
            int tetra;
            //
            //  Each order 4 tetrahedron defines 6 adjacency pairs.
            //
            pair = new int[2 * 6 * tetra_num];

            for (tetra = 0; tetra < tetra_num; tetra++)
            {
                pair[0 + tetra * 2] = tetra_node[0 + tetra * 4];
                pair[1 + tetra * 2] = tetra_node[1 + tetra * 4];

                pair[0 + (tetra_num + tetra) * 2] = tetra_node[0 + tetra * 4];
                pair[1 + (tetra_num + tetra) * 2] = tetra_node[2 + tetra * 4];

                pair[0 + (2 * tetra_num + tetra) * 2] = tetra_node[0 + tetra * 4];
                pair[1 + (2 * tetra_num + tetra) * 2] = tetra_node[3 + tetra * 4];

                pair[0 + (3 * tetra_num + tetra) * 2] = tetra_node[1 + tetra * 4];
                pair[1 + (3 * tetra_num + tetra) * 2] = tetra_node[2 + tetra * 4];

                pair[0 + (4 * tetra_num + tetra) * 2] = tetra_node[1 + tetra * 4];
                pair[1 + (4 * tetra_num + tetra) * 2] = tetra_node[3 + tetra * 4];

                pair[0 + (5 * tetra_num + tetra) * 2] = tetra_node[2 + tetra * 4];
                pair[1 + (5 * tetra_num + tetra) * 2] = tetra_node[3 + tetra * 4];
            }

            pair_num = 6 * tetra_num;
            //
            //  Force the nodes of each pair to be listed in ascending order.
            //
            typeMethods.i4mat_transpose_print_some(2, pair_num, pair, 1, 1, 2, pair_num,
                "DEBUG: PAIR before first sort");

            typeMethods.i4col_sort2_a(2, pair_num, ref pair);

            typeMethods.i4mat_transpose_print_some(2, pair_num, pair, 1, 1, 2, pair_num,
                "DEBUG: PAIR after first sort");
            //
            //  Rearrange the columns in ascending order.
            //
            typeMethods.i4col_sort_a(2, pair_num, ref pair);
            //
            //  Get the number of unique columns.
            //
            pair_unique_num = typeMethods.i4col_sorted_unique_count(2, pair_num, pair);
            //
            //  The number of adjacencies is TWICE this value, plus the number of nodes.
            //
            adj_num = 2 * pair_unique_num;
            //
            //  Now set up the ADJ_ROW counts.
            //
            for (node = 0; node < node_num; node++)
            {
                adj_row[node] = 0;
            }

            for (k = 0; k < pair_num; k++)
            {
                if (0 < k)
                {
                    if (pair[0 + (k - 1) * 2] == pair[0 + k * 2] &&
                        pair[1 + (k - 1) * 2] == pair[1 + k * 2])
                    {
                        continue;
                    }
                }

                i = pair[0 + k * 2];
                j = pair[1 + k * 2];

                adj_row[i - 1] = adj_row[i - 1] + 1;
                adj_row[j - 1] = adj_row[j - 1] + 1;
            }

            //
            //  We used ADJ_ROW to count the number of entries in each row.
            //  Convert it to pointers into the ADJ array.
            //
            for (node = node_num - 1; 0 <= node; node--)
            {
                adj_row[node] = adj_row[node + 1];
            }

            adj_row[0] = 1;
            for (node = 1; node <= node_num; node++)
            {
                adj_row[node] = adj_row[node - 1] + adj_row[i];
            }
        }

        public static int[] tet_mesh_order4_adj_set(int node_num, int element_num,
                int[] element_node, int adj_num, int[] adj_row)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TET_MESH_ORDER4_ADJ_SET sets the nodal adjacency matrix.
            //
            //  Discussion:
            //
            //    A compressed format is used for the nodal adjacency matrix.
            //
            //    It is assumed that we know ADJ_NUM, the number of adjacency entries
            //    and the ADJ_ROW array, which keeps track of the list of slots
            //    in ADJ where we can store adjacency information for each row.
            //
            //    We essentially repeat the work of TET_MESH_ORDER4_ADJ_COUNT, but
            //    now we have a place to store the adjacency information.
            //
            //    A copy of the ADJ_ROW array is useful, as we can use it to keep track
            //    of the next available entry in ADJ for adjacencies associated with
            //    a given row.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    13 January 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NODE_NUM, the number of nodes.
            //
            //    Input, int TETRA_NUM, the number of tetrahedrons.
            //
            //    Input, int TETRA_NODE[4*TETRA_NUM], the indices of the nodes.
            //
            //    Input, int ADJ_NUM, the total number of adjacency relationships,
            //
            //    Input, int ADJ_ROW[NODE_NUM+1], the ADJ pointer array.
            //
            //    Output, int TET_MESH_ORDER4_ADJ_SET[ADJ_NUM], 
            //    the adjacency information.
            //
        {
            int[] adj;
            int[] adj_row_copy;
            int i;
            int j;
            int k;
            int node;
            int[] pair;
            int pair_num;
            int tetra;
            //
            //  Each order 4 tetrahedron defines 6 adjacency pairs.
            //
            pair = new int[2 * 6 * element_num];

            for (tetra = 0; tetra < element_num; tetra++)
            {
                pair[0 + tetra * 2] = element_node[0 + tetra * 4];
                pair[1 + tetra * 2] = element_node[1 + tetra * 4];

                pair[0 + (element_num + tetra) * 2] = element_node[0 + tetra * 4];
                pair[1 + (element_num + tetra) * 2] = element_node[2 + tetra * 4];

                pair[0 + (2 * element_num + tetra) * 2] = element_node[0 + tetra * 4];
                pair[1 + (2 * element_num + tetra) * 2] = element_node[3 + tetra * 4];

                pair[0 + (3 * element_num + tetra) * 2] = element_node[1 + tetra * 4];
                pair[1 + (3 * element_num + tetra) * 2] = element_node[2 + tetra * 4];

                pair[0 + (4 * element_num + tetra) * 2] = element_node[1 + tetra * 4];
                pair[1 + (4 * element_num + tetra) * 2] = element_node[3 + tetra * 4];

                pair[0 + (5 * element_num + tetra) * 2] = element_node[2 + tetra * 4];
                pair[1 + (5 * element_num + tetra) * 2] = element_node[3 + tetra * 4];
            }

            pair_num = 6 * element_num;
            //
            //  Force the nodes of each pair to be listed in ascending order.
            //
            typeMethods.i4col_sort2_a(2, pair_num, ref pair);
            //
            //  Rearrange the columns in ascending order.
            //
            typeMethods.i4col_sort_a(2, pair_num, ref pair);
            //
            //  Mark all entries of ADJ so we will know later if we missed one.
            //
            adj = new int[adj_num];

            for (i = 0; i < adj_num; i++)
            {
                adj[i] = -1;
            }

            //
            //  Copy the ADJ_ROW array and use it to keep track of the next
            //  free entry for each row.
            //
            adj_row_copy = new int[node_num];

            for (node = 0; node < node_num; node++)
            {
                adj_row_copy[node] = adj_row[node];
            }

            //
            //  Now set up the ADJ_ROW counts.
            //
            for (k = 0; k < pair_num; k++)
            {
                if (0 < k)
                {
                    if (pair[0 + (k - 1) * 2] == pair[0 + k * 2] &&
                        pair[1 + (k - 1) * 2] == pair[1 + k * 2])
                    {
                        continue;
                    }
                }

                i = pair[0 + k * 2];
                j = pair[1 + k * 2];

                adj[adj_row_copy[i]] = j;
                adj_row_copy[i] = adj_row_copy[i] + 1;
                adj[adj_row_copy[j]] = i;
                adj_row_copy[j] = adj_row_copy[j] + 1;
            }

            return adj;
        }

        public static int tet_mesh_order4_boundary_face_count(int tetra_num, int[] tetra_node)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TET_MESH_ORDER4_BOUNDARY_FACE_COUNT counts the number of boundary faces.
            //
            //  Discussion:
            //
            //    This routine is given a tet mesh, an abstract list of 
            //    quadruples of nodes.  It is assumed that the nodes forming each 
            //    face of each tetrahedron are listed in a counterclockwise order, 
            //    although the routine should work if the nodes are consistently 
            //    listed in a clockwise order as well.
            //
            //    It is assumed that each face of the tet mesh is either
            //    * an INTERIOR face, which is listed twice, once with positive
            //      orientation and once with negative orientation, or;
            //    * a BOUNDARY face, which will occur only once.
            //
            //    This routine should work even if the region has holes.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    11 January 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int TETRA_NUM, the number of tetrahedrons.
            //
            //    Input, int TETRA_NODE[4*TETRA_NUM], the indices of the nodes.
            //
            //    Output, int TET_MESH_ORDER4_BOUNDARY_FACE_COUNT, the number of 
            //    boundary faces.
            //
        {
            int boundary_face_num;
            int[] face;
            int face_num;
            int interior_face_num;
            int m;
            int tet;
            int unique_face_num;

            face = new int[3 * 4 * tetra_num];

            m = 3;
            face_num = 4 * tetra_num;
            //
            //  Set up the face array:
            //  (Omit node 1)
            //  (Omit node 2)
            //  (Omit node 3)
            //  (Omit node 4)
            //
            for (tet = 0; tet < tetra_num; tet++)
            {
                face[0 + (tet) * 3] = tetra_node[1 + tet * 4];
                face[1 + (tet) * 3] = tetra_node[2 + tet * 4];
                face[2 + (tet) * 3] = tetra_node[3 + tet * 4];

                face[0 + (tetra_num + tet) * 3] = tetra_node[0 + tet * 4];
                face[1 + (tetra_num + tet) * 3] = tetra_node[2 + tet * 4];
                face[2 + (tetra_num + tet) * 3] = tetra_node[3 + tet * 4];

                face[0 + (2 * tetra_num + tet) * 3] = tetra_node[0 + tet * 4];
                face[1 + (2 * tetra_num + tet) * 3] = tetra_node[1 + tet * 4];
                face[2 + (2 * tetra_num + tet) * 3] = tetra_node[3 + tet * 4];

                face[0 + (3 * tetra_num + tet) * 3] = tetra_node[0 + tet * 4];
                face[1 + (3 * tetra_num + tet) * 3] = tetra_node[1 + tet * 4];
                face[2 + (3 * tetra_num + tet) * 3] = tetra_node[2 + tet * 4];
            }

            //
            //  Force the nodes of each face to be listed in ascending order.
            //
            typeMethods.i4col_sort2_a(m, face_num, ref face);
            //
            //  Ascending sort the columns.
            //
            typeMethods.i4col_sort_a(m, face_num, ref face);
            //
            //  Get the number of unique columns.
            //
            unique_face_num = typeMethods.i4col_sorted_unique_count(m, face_num, face);
            //
            //  Determine the number of interior and boundary faces.
            //
            interior_face_num = 4 * tetra_num - unique_face_num;

            boundary_face_num = 4 * tetra_num - 2 * interior_face_num;

            return boundary_face_num;
        }

        public static int tet_mesh_order4_edge_count(int tetra_num, int[] tetra_node)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TET_MESH_ORDER4_EDGE_COUNT counts the number of edges.
            //
            //  Discussion:
            //
            //    This routine is given a tet mesh, an abstract list of
            //    quadruples of nodes.  Each tetrahedron defines 6 edges; however,
            //    assuming that tetrahedrons are touching each other, most edges
            //    will be used more than once.  This routine determines the actual
            //    number of "geometric" edges associated with the tet mesh.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    11 January 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int TETRA_NUM, the number of tetrahedrons.
            //
            //    Input, int TETRA_NODE[4*TETRA_NUM], the indices of the nodes.
            //
            //    Output, int TET_MESH_ORDER4_EDGE_COUNT, the number of edges.
            //
        {
            int[] edge;
            int edge_num;
            int edge_num_raw;
            int m;
            int tet;

            edge = new int[2 * 6 * tetra_num];

            m = 3;
            edge_num_raw = 6 * tetra_num;
            //
            //  Set up the raw edge array:
            //
            for (tet = 0; tet < tetra_num; tet++)
            {
                edge[0 + tet * 2] = tetra_node[0 + tet * 4];
                edge[1 + tet * 2] = tetra_node[1 + tet * 4];

                edge[0 + (tetra_num + tet) * 2] = tetra_node[0 + tet * 4];
                edge[1 + (tetra_num + tet) * 2] = tetra_node[2 + tet * 4];

                edge[0 + (2 * tetra_num + tet) * 2] = tetra_node[0 + tet * 4];
                edge[1 + (2 * tetra_num + tet) * 2] = tetra_node[3 + tet * 4];

                edge[0 + (3 * tetra_num + tet) * 2] = tetra_node[1 + tet * 4];
                edge[1 + (3 * tetra_num + tet) * 2] = tetra_node[2 + tet * 4];

                edge[0 + (4 * tetra_num + tet) * 2] = tetra_node[1 + tet * 4];
                edge[1 + (4 * tetra_num + tet) * 2] = tetra_node[3 + tet * 4];

                edge[0 + (5 * tetra_num + tet) * 2] = tetra_node[2 + tet * 4];
                edge[1 + (5 * tetra_num + tet) * 2] = tetra_node[3 + tet * 4];
            }

            //
            //  Force the nodes of each face to be listed in ascending order.
            //
            typeMethods.i4col_sort2_a(m, edge_num_raw, ref edge);
            //
            //  Ascending sort the columns.
            //
            typeMethods.i4col_sort_a(m, edge_num_raw, ref edge);
            //
            //  Get the number of unique columns.
            //
            edge_num = typeMethods.i4col_sorted_unique_count(m, edge_num_raw, edge);

            return edge_num;
        }

        public static void tet_mesh_order4_example_set(int node_num, int tetra_num,
                ref double[] node_xyz, ref int[] tetra_node)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TET_MESH_ORDER4_EXAMPLE_SET sets an example linear tet mesh.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 August 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NODE_NUM, the number of nodes.
            //
            //    Input, int TETRA_NUM, the number of tetrahedrons.
            //
            //    Output, double NODE_XYZ[3*NODE_NUM], the node coordinates.
            //
            //    Output, int TETRA_NODE[4*TETRA_NUM], the nodes forming each tet.
            //
        {
            int i;
            int j;
            double[] node_xyz_save =
                {
                    0.0, 0.0, 0.0,
                    0.0, 0.0, 0.5,
                    0.0, 0.0, 1.0,
                    0.0, 0.5, 0.0,
                    0.0, 0.5, 0.5,
                    0.0, 0.5, 1.0,
                    0.0, 1.0, 0.0,
                    0.0, 1.0, 0.5,
                    0.0, 1.0, 1.0,
                    0.5, 0.0, 0.0,
                    0.5, 0.0, 0.5,
                    0.5, 0.0, 1.0,
                    0.5, 0.5, 0.0,
                    0.5, 0.5, 0.5,
                    0.5, 0.5, 1.0,
                    0.5, 1.0, 0.0,
                    0.5, 1.0, 0.5,
                    0.5, 1.0, 1.0,
                    1.0, 0.0, 0.0,
                    1.0, 0.0, 0.5,
                    1.0, 0.0, 1.0,
                    1.0, 0.5, 0.0,
                    1.0, 0.5, 0.5,
                    1.0, 0.5, 1.0,
                    1.0, 1.0, 0.0,
                    1.0, 1.0, 0.5,
                    1.0, 1.0, 1.0,
                    1.5, 0.0, 0.0,
                    1.5, 0.0, 0.5,
                    1.5, 0.0, 1.0,
                    1.5, 0.5, 0.0,
                    1.5, 0.5, 0.5,
                    1.5, 0.5, 1.0,
                    1.5, 1.0, 0.0,
                    1.5, 1.0, 0.5,
                    1.5, 1.0, 1.0,
                    2.0, 0.0, 0.0,
                    2.0, 0.0, 0.5,
                    2.0, 0.0, 1.0,
                    2.0, 0.5, 0.0,
                    2.0, 0.5, 0.5,
                    2.0, 0.5, 1.0,
                    2.0, 1.0, 0.0,
                    2.0, 1.0, 0.5,
                    2.0, 1.0, 1.0,
                    2.5, 0.0, 0.0,
                    2.5, 0.0, 0.5,
                    2.5, 0.0, 1.0,
                    2.5, 0.5, 0.0,
                    2.5, 0.5, 0.5,
                    2.5, 0.5, 1.0,
                    2.5, 1.0, 0.0,
                    2.5, 1.0, 0.5,
                    2.5, 1.0, 1.0,
                    3.0, 0.0, 0.0,
                    3.0, 0.0, 0.5,
                    3.0, 0.0, 1.0,
                    3.0, 0.5, 0.0,
                    3.0, 0.5, 0.5,
                    3.0, 0.5, 1.0,
                    3.0, 1.0, 0.0,
                    3.0, 1.0, 0.5,
                    3.0, 1.0, 1.0
                }
                ;
            int[] tetra_node_save =
                {
                    1, 2, 4, 10,
                    2, 4, 5, 10,
                    2, 5, 10, 11,
                    2, 3, 5, 11,
                    4, 5, 10, 13,
                    3, 5, 6, 11,
                    5, 10, 11, 13,
                    4, 5, 7, 13,
                    5, 6, 8, 14,
                    5, 7, 8, 13,
                    6, 8, 9, 14,
                    11, 13, 14, 19,
                    12, 14, 15, 20,
                    3, 6, 11, 12,
                    5, 6, 11, 14,
                    6, 9, 14, 15,
                    6, 11, 12, 14,
                    6, 12, 14, 15,
                    7, 8, 13, 16,
                    5, 8, 13, 14,
                    10, 11, 13, 19,
                    8, 9, 14, 17,
                    11, 12, 14, 20,
                    5, 11, 13, 14,
                    8, 13, 14, 16,
                    9, 14, 15, 17,
                    13, 14, 16, 22,
                    8, 14, 16, 17,
                    14, 15, 17, 23,
                    14, 16, 17, 22,
                    9, 15, 17, 18,
                    15, 17, 18, 23,
                    14, 17, 22, 23,
                    13, 14, 19, 22,
                    11, 14, 19, 20,
                    14, 15, 20, 23,
                    15, 20, 21, 23,
                    21, 23, 24, 29,
                    20, 22, 23, 28,
                    14, 19, 20, 22,
                    15, 18, 23, 24,
                    12, 15, 20, 21,
                    15, 21, 23, 24,
                    16, 17, 22, 25,
                    19, 20, 22, 28,
                    17, 18, 23, 26,
                    20, 21, 23, 29,
                    14, 20, 22, 23,
                    17, 22, 23, 25,
                    18, 23, 24, 26,
                    22, 23, 25, 31,
                    17, 23, 25, 26,
                    23, 24, 26, 32,
                    23, 25, 26, 31,
                    18, 24, 26, 27,
                    24, 26, 27, 32,
                    23, 26, 31, 32,
                    22, 23, 28, 31,
                    20, 23, 28, 29,
                    23, 24, 29, 32,
                    24, 29, 30, 32,
                    30, 32, 33, 38,
                    29, 31, 32, 37,
                    23, 28, 29, 31,
                    24, 27, 32, 33,
                    21, 24, 29, 30,
                    24, 30, 32, 33,
                    25, 26, 31, 34,
                    28, 29, 31, 37,
                    26, 27, 32, 35,
                    29, 30, 32, 38,
                    23, 29, 31, 32,
                    26, 31, 32, 34,
                    27, 32, 33, 35,
                    31, 32, 34, 40,
                    26, 32, 34, 35,
                    32, 33, 35, 41,
                    32, 34, 35, 40,
                    27, 33, 35, 36,
                    33, 35, 36, 41,
                    32, 35, 40, 41,
                    31, 32, 37, 40,
                    29, 32, 37, 38,
                    32, 33, 38, 41,
                    33, 38, 39, 41,
                    39, 41, 42, 47,
                    38, 40, 41, 46,
                    32, 37, 38, 40,
                    33, 36, 41, 42,
                    30, 33, 38, 39,
                    33, 39, 41, 42,
                    34, 35, 40, 43,
                    37, 38, 40, 46,
                    35, 36, 41, 44,
                    38, 39, 41, 47,
                    32, 38, 40, 41,
                    35, 40, 41, 43,
                    36, 41, 42, 44,
                    40, 41, 43, 49,
                    35, 41, 43, 44,
                    41, 42, 44, 50,
                    41, 43, 44, 49,
                    36, 42, 44, 45,
                    42, 44, 45, 50,
                    41, 44, 49, 50,
                    40, 41, 46, 49,
                    38, 41, 46, 47,
                    41, 42, 47, 50,
                    42, 47, 48, 50,
                    48, 50, 51, 56,
                    47, 49, 50, 55,
                    41, 46, 47, 49,
                    42, 45, 50, 51,
                    39, 42, 47, 48,
                    42, 48, 50, 51,
                    43, 44, 49, 52,
                    46, 47, 49, 55,
                    44, 45, 50, 53,
                    47, 48, 50, 56,
                    41, 47, 49, 50,
                    44, 49, 50, 52,
                    45, 50, 51, 53,
                    49, 50, 52, 58,
                    44, 50, 52, 53,
                    50, 51, 53, 59,
                    50, 52, 53, 58,
                    45, 51, 53, 54,
                    51, 53, 54, 59,
                    50, 53, 58, 59,
                    49, 50, 55, 58,
                    47, 50, 55, 56,
                    50, 51, 56, 59,
                    51, 56, 57, 59,
                    50, 55, 56, 58,
                    51, 54, 59, 60,
                    48, 51, 56, 57,
                    51, 57, 59, 60,
                    52, 53, 58, 61,
                    53, 54, 59, 62,
                    50, 56, 58, 59,
                    53, 58, 59, 61,
                    54, 59, 60, 62,
                    53, 59, 61, 62,
                    54, 60, 62, 63
                }
                ;

            for (j = 0; j < node_num; j++)
            {
                for (i = 0; i < 3; i++)
                {
                    node_xyz[i + j * 3] = node_xyz_save[i + j * 3];
                }
            }

            for (j = 0; j < tetra_num; j++)
            {
                for (i = 0; i < 4; i++)
                {
                    tetra_node[i + j * 4] = tetra_node_save[i + j * 4] - 1;
                }
            }
        }

        public static void tet_mesh_order4_example_size(ref int node_num, ref int tetra_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TET_MESH_ORDER4_EXAMPLE_SIZE sizes an example linear tet mesh.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    03 August 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, int *NODE_NUM, the number of nodes.
            //
            //    Output, int *TETRA_NUM, the number of tetrahedrons.
            //
        {
            node_num = 63;
            tetra_num = 144;

        }

        public static void tet_mesh_order4_refine_compute(int node_num1, int tetra_num1,
                double[] node_xyz1, int[] tetra_node1, int node_num2, int tetra_num2,
                int[] edge_data, ref double[] node_xyz2, ref int[] tetra_node2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TET_MESH_ORDER4_REFINE_COMPUTE computes a refined order 4 tet mesh
            //
            //  Discussion:
            //
            //    A refined 4-node tet mesh can be derived from a given
            //    4-node tet mesh by interpolating nodes at the midpoint of
            //    every edge of the mesh.
            //
            //    The mesh is described indirectly, as the sum of individual
            //    tetrahedrons.  A single physical edge may be a logical edge of
            //    any number of tetrahedrons.  It is important, however, that a
            //    new node be created exactly once for each edge, assigned an index,
            //    and associated with every tetrahedron that shares this edge. 
            //
            //    This routine handles that problem.
            //
            //    The primary amount of work occurs in sorting a list of 6 * TETRA_NUM
            //    data items, one item for every edge of every tetrahedron.  Each
            //    data item records, for a given tetrahedron edge, the global indices
            //    of the two endpoints, the local indices of the two endpoints,
            //    and the index of the tetrahedron.
            //
            //    Through careful sorting, it is possible to arrange this data in
            //    a way that allows the proper generation of the interpolated nodes.
            //
            //    Let us add the new nodes and temporarily assign them local indices
            //    5 through X, based on the following ordering:
            //
            //      1, 2, 3, 4, (1+2), (1+3), (1+4), (2+3), (2+4), (3+4).
            //
            //    Then let us assign these nodes to eight subtetrahedrons as follows:
            //
            //      1, 5, 6, 7
            //      2, 5, 8, 9
            //      3, 6, 8, 9
            //      4, 7, 9, X
            //      5, 6, 7, 9
            //      5, 6, 8, 9
            //      6, 7, 9, X
            //      6, 8, 9, X
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    25 January 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Anwei Liu, Barry Joe,
            //    Quality Local Refinement of Tetrahedral Meshes Based
            //    on 8-Subtetrahedron Subdivision,
            //    Mathematics of Computation,
            //    Volume 65, Number 215, July 1996, pages 1183-1200.
            //
            //  Parameters:
            //
            //    Input, int NODE_NUM1, the number of nodes in the input mesh.
            //
            //    Input, int TETRA_NUM1, the number of tetrahedrons in the
            //    input mesh.
            //
            //    Input, double NODE_XYZ1[3*NODE_NUM1], the coordinates of
            //    the nodes that make up the input mesh.
            //
            //    Input, int TETRA_NODE1[4*TETRA_NUM], the indices of the nodes
            //    in the input mesh.
            //
            //    Input, int NODE_NUM2, the number of nodes in the refined mesh.
            //
            //    Input, int TETRA_NUM2, the number of tetrahedrons in the
            //    refined mesh.
            //
            //    Input, int EDGE_DATA[5*(6*TETRA_NUM1)], edge data.
            //
            //    Output, double NODE_XYZ2[3*NODE_NUM2], the coordinates of
            //    the nodes that make up the refined mesh.
            //
            //    Output, int TETRA_NODE2[4*TETRA_NUM2], the indices of the nodes 
            //    in the refined mesh.
            //
        {
            int dim_num = 3;
            int edge;
            int i;
            int j;
            int n1;
            int n1_old;
            int n2;
            int n2_old;
            int node;
            int tetra_order = 4;
            int tetra1;
            int v1;
            int v2;
            //
            //  Generate the index and coordinates of the new midside nodes, 
            //  and update the tetradehron-node data.
            //
            for (j = 0; j < node_num1; j++)
            {
                for (i = 0; i < dim_num; i++)
                {
                    node_xyz2[i + j * dim_num] = node_xyz1[i + j * dim_num];
                }
            }

            for (j = 0; j < tetra_num2; j++)
            {
                for (i = 0; i < tetra_order; i++)
                {
                    tetra_node2[i + j * tetra_order] = -1;
                }
            }

            //
            //  The vertices of the input tetrahedron can be assigned now.
            //
            for (tetra1 = 0; tetra1 < tetra_num1; tetra1++)
            {
                tetra_node2[0 + (tetra1 * 8 + 0) * tetra_order] = tetra_node1[0 + tetra1 * tetra_order];
                tetra_node2[0 + (tetra1 * 8 + 1) * tetra_order] = tetra_node1[1 + tetra1 * tetra_order];
                tetra_node2[0 + (tetra1 * 8 + 2) * tetra_order] = tetra_node1[2 + tetra1 * tetra_order];
                tetra_node2[0 + (tetra1 * 8 + 3) * tetra_order] = tetra_node1[3 + tetra1 * tetra_order];
            }

            node = node_num1;

            n1_old = -1;
            n2_old = -1;

            for (edge = 0; edge < 6 * tetra_num1; edge++)
            {
                //
                //  Read the data defining the edge.
                //
                n1 = edge_data[0 + edge * 5];
                n2 = edge_data[1 + edge * 5];
                //
                //  If this edge is new, create the coordinates and index.
                //
                if (n1 != n1_old || n2 != n2_old)
                {
                    if (node_num2 <= node)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("TET_MESH_ORDER4_REFINE_COMPUTE - Fatal error!");
                        Console.WriteLine("  Node index exceeds NODE_NUM2.");
                        return;
                    }

                    for (i = 0; i < dim_num; i++)
                    {
                        node_xyz2[i + node * dim_num] =
                            (node_xyz2[i + (n1 - 1) * dim_num] + node_xyz2[i + (n2 - 1) * dim_num]) / 2.0;
                    }

                    node = node + 1;
                    n1_old = n1;
                    n2_old = n2;
                }

                //
                //  Assign the node to the tetrahedron.
                //
                v1 = edge_data[2 + edge * 5];
                v2 = edge_data[3 + edge * 5];
                tetra1 = edge_data[4 + edge * 5];
                //
                //  We know the two vertices that bracket this new node.
                //  This tells us whether it is new node number 5, 6, 7, 8, 9 or 10.
                //  This tells us which of the new subtetrahedrons it belongs to,
                //  and what position it occupies.
                //
                if (v1 == 1 && v2 == 2)
                {
                    tetra_node2[1 + (tetra1 * 8 + 0) * tetra_order] = node;
                    tetra_node2[1 + (tetra1 * 8 + 1) * tetra_order] = node;
                    tetra_node2[0 + (tetra1 * 8 + 4) * tetra_order] = node;
                    tetra_node2[0 + (tetra1 * 8 + 5) * tetra_order] = node;
                }
                else if (v1 == 1 && v2 == 3)
                {
                    tetra_node2[2 + (tetra1 * 8 + 0) * tetra_order] = node;
                    tetra_node2[1 + (tetra1 * 8 + 2) * tetra_order] = node;
                    tetra_node2[1 + (tetra1 * 8 + 4) * tetra_order] = node;
                    tetra_node2[1 + (tetra1 * 8 + 5) * tetra_order] = node;
                    tetra_node2[0 + (tetra1 * 8 + 6) * tetra_order] = node;
                    tetra_node2[0 + (tetra1 * 8 + 7) * tetra_order] = node;
                }
                else if (v1 == 1 && v2 == 4)
                {
                    tetra_node2[3 + (tetra1 * 8 + 0) * tetra_order] = node;
                    tetra_node2[1 + (tetra1 * 8 + 3) * tetra_order] = node;
                    tetra_node2[2 + (tetra1 * 8 + 4) * tetra_order] = node;
                    tetra_node2[1 + (tetra1 * 8 + 6) * tetra_order] = node;
                }
                else if (v1 == 2 && v2 == 3)
                {
                    tetra_node2[2 + (tetra1 * 8 + 1) * tetra_order] = node;
                    tetra_node2[2 + (tetra1 * 8 + 2) * tetra_order] = node;
                    tetra_node2[2 + (tetra1 * 8 + 5) * tetra_order] = node;
                    tetra_node2[1 + (tetra1 * 8 + 7) * tetra_order] = node;
                }
                else if (v1 == 2 && v2 == 4)
                {
                    tetra_node2[3 + (tetra1 * 8 + 1) * tetra_order] = node;
                    tetra_node2[3 + (tetra1 * 8 + 2) * tetra_order] = node;
                    tetra_node2[2 + (tetra1 * 8 + 3) * tetra_order] = node;
                    tetra_node2[3 + (tetra1 * 8 + 4) * tetra_order] = node;
                    tetra_node2[3 + (tetra1 * 8 + 5) * tetra_order] = node;
                    tetra_node2[2 + (tetra1 * 8 + 6) * tetra_order] = node;
                    tetra_node2[2 + (tetra1 * 8 + 7) * tetra_order] = node;
                }
                else if (v1 == 3 && v2 == 4)
                {
                    tetra_node2[3 + (tetra1 * 8 + 3) * tetra_order] = node;
                    tetra_node2[3 + (tetra1 * 8 + 6) * tetra_order] = node;
                    tetra_node2[3 + (tetra1 * 8 + 7) * tetra_order] = node;
                }
            }

            return;
        }

        public static void tet_mesh_order4_refine_size(int node_num1, int tetra_num1,
                int[] tetra_node1, ref int node_num2, ref int tetra_num2, ref int[] edge_data)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TET_MESH_ORDER4_REFINE_SIZE sizes a refined order 4 tet mesh.
            //
            //  Discussion:
            //
            //    A refined tet mesh can be derived from an existing one by interpolating 
            //    nodes at the midpoint of every edge of the mesh.
            //
            //    The mesh is described indirectly, as the sum of individual
            //    tetrahedrons.  A single physical edge may be a logical edge of
            //    any number of tetrahedrons.  It is important, however, that a
            //    new node be created exactly once for each edge, assigned an index,
            //    and associated with every tetrahedron that shares this edge. 
            //
            //    This routine handles that problem.
            //
            //    The primary amount of work occurs in sorting a list of 6 * TETRA_NUM
            //    data items, one item for every edge of every tetrahedron.  Each
            //    data item records, for a given tetrahedron edge, the global indices
            //    of the two endpoints, the local indices of the two endpoints,
            //    and the index of the tetrahedron.
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
            //    25 January 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NODE_NUM1, the number of nodes in the original mesh.
            //
            //    Input, int TETRA_NUM1, the number of tetrahedrons in the
            //    original mesh.
            //
            //    Input, int TETRA_NODE1[4*TETRA_NUM1], the indices of the nodes
            //    in the original mesh.
            //
            //    Output, int *NODE_NUM2, the number of nodes in the refined mesh.
            //
            //    Output, int *TETRA_NUM2, the number of tetrahedrons in the refined mesh.
            //
            //    Output, int EDGE_DATA[5*(6*TETRA_NUM1)], edge data.
            //
        {
            int a = 0;
            int b = 0;
            int edge;
            int i;
            int j;
            int k;
            int l;
            int n1;
            int n1_old;
            int n2;
            int n2_old;
            int tetra;
            int tetra_order = 4;
            //
            //  Step 1.
            //  From the list of nodes for tetrahedron T, of the form: (I,J,K,L)
            //  construct the six edge relations:
            //
            //    (I,J,1,2,T)
            //    (I,K,1,3,T)
            //    (I,L,1,4,T)
            //    (J,K,2,3,T)
            //    (J,L,2,4,T)
            //    (K,L,3,4,T)
            //
            //  In order to make matching easier, we reorder each pair of nodes
            //  into ascending order.
            //
            for (tetra = 0; tetra < tetra_num1; tetra++)
            {
                i = tetra_node1[0 + tetra * tetra_order];
                j = tetra_node1[1 + tetra * tetra_order];
                k = tetra_node1[2 + tetra * tetra_order];
                l = tetra_node1[3 + tetra * tetra_order];

                typeMethods.i4i4_sort_a(i, j, ref a, ref b);

                edge_data[0 + (6 * tetra) * 5] = a;
                edge_data[1 + (6 * tetra) * 5] = b;
                edge_data[2 + (6 * tetra) * 5] = 1;
                edge_data[3 + (6 * tetra) * 5] = 2;
                edge_data[4 + (6 * tetra) * 5] = tetra;

                typeMethods.i4i4_sort_a(i, k, ref a, ref b);

                edge_data[0 + (6 * tetra + 1) * 5] = a;
                edge_data[1 + (6 * tetra + 1) * 5] = b;
                edge_data[2 + (6 * tetra + 1) * 5] = 1;
                edge_data[3 + (6 * tetra + 1) * 5] = 3;
                edge_data[4 + (6 * tetra + 1) * 5] = tetra;

                typeMethods.i4i4_sort_a(i, l, ref a, ref b);

                edge_data[0 + (6 * tetra + 2) * 5] = a;
                edge_data[1 + (6 * tetra + 2) * 5] = b;
                edge_data[2 + (6 * tetra + 2) * 5] = 1;
                edge_data[3 + (6 * tetra + 2) * 5] = 4;
                edge_data[4 + (6 * tetra + 2) * 5] = tetra;

                typeMethods.i4i4_sort_a(j, k, ref a, ref b);

                edge_data[0 + (6 * tetra + 3) * 5] = a;
                edge_data[1 + (6 * tetra + 3) * 5] = b;
                edge_data[2 + (6 * tetra + 3) * 5] = 2;
                edge_data[3 + (6 * tetra + 3) * 5] = 3;
                edge_data[4 + (6 * tetra + 3) * 5] = tetra;

                typeMethods.i4i4_sort_a(j, l, ref a, ref b);

                edge_data[0 + (6 * tetra + 4) * 5] = a;
                edge_data[1 + (6 * tetra + 4) * 5] = b;
                edge_data[2 + (6 * tetra + 4) * 5] = 2;
                edge_data[3 + (6 * tetra + 4) * 5] = 4;
                edge_data[4 + (6 * tetra + 4) * 5] = tetra;

                typeMethods.i4i4_sort_a(k, l, ref a, ref b);

                edge_data[0 + (6 * tetra + 5) * 5] = a;
                edge_data[1 + (6 * tetra + 5) * 5] = b;
                edge_data[2 + (6 * tetra + 5) * 5] = 3;
                edge_data[3 + (6 * tetra + 5) * 5] = 4;
                edge_data[4 + (6 * tetra + 5) * 5] = tetra;
            }

            //
            //  Step 2. Perform an ascending dictionary sort on the neighbor relations.
            //  We only intend to sort on rows 1:2; the routine we call here
            //  sorts on the full column but that won't hurt us.
            //
            //  What we need is to find all cases where tetrahedrons share an edge.
            //  By sorting the columns of the EDGE_DATA array, we will put shared edges
            //  next to each other.
            //
            typeMethods.i4col_sort_a(5, 6 * tetra_num1, ref edge_data);
            //
            //  Step 3. All the tetrahedrons which share an edge show up as consecutive
            //  columns with identical first two entries.  Figure out how many new
            //  nodes there are, and allocate space for their coordinates.
            //
            node_num2 = node_num1;

            n1_old = -1;
            n2_old = -1;

            for (edge = 0; edge < 6 * tetra_num1; edge++)
            {
                n1 = edge_data[0 + edge * 5];
                n2 = edge_data[1 + edge * 5];
                if (n1 != n1_old || n2 != n2_old)
                {
                    node_num2 = node_num2 + 1;
                    n1_old = n1;
                    n2_old = n2;
                }
            }

            tetra_num2 = 8 * tetra_num1;
        }

        public static void tet_mesh_order4_to_order10_compute(int tetra_num, int[] tetra_node1,
                int node_num1, double[] node_xyz1, int[] edge_data, ref int[] tetra_node2,
                int node_num2, ref double[] node_xyz2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TET_MESH_ORDER4_TO_ORDER10_COMPUTE computes a quadratic tet mesh from a linear one.
            //
            //  Discussion:
            //
            //    A quadratic (10 node) tet mesh can be derived from a linear
            //    (4 node) tet mesh by interpolating nodes at the midpoint of
            //    every edge of the mesh.
            //
            //    The mesh is described indirectly, as the sum of individual
            //    tetrahedrons.  A single physical edge may be a logical edge of
            //    any number of tetrahedrons.  It is important, however, that a
            //    new node be created exactly once for each edge, assigned an index,
            //    and associated with every tetrahedron that shares this edge. 
            //
            //    This routine handles that problem.
            //
            //    The primary amount of work occurs in sorting a list of 6 * TETRA_NUM
            //    data items, one item for every edge of every tetrahedron.  Each
            //    data item records, for a given tetrahedron edge, the global indices
            //    of the two endpoints, the local indices of the two endpoints,
            //    and the index of the tetrahedron.
            //
            //    Through careful sorting, it is possible to arrange this data in
            //    a way that allows the proper generation of the interpolated nodes.
            //
            //    The node ordering for the quadratic tetrahedron is somewhat
            //    arbitrary.  In the current scheme, the vertices are listed
            //    first, followed by the 6 midside nodes.  Each midside node
            //    may be identified by the two vertices that bracket it.  Thus,
            //    the node ordering may be suggested by:
            //
            //      1  2  3  4 (1+2) (1+3) (1+4) (2+3) (2+4) (3+4)
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 January 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int TETRA_NUM, the number of tetrahedrons in the
            //    linear mesh.
            //
            //    Input, int TETRA_NODE1[4*TETRA_NUM], the indices of the nodes
            //    in the linear mesh.
            //
            //    Input, int NODE_NUM1, the number of nodes for the linear mesh.
            //
            //    Input, double NODE_XYZ1[3*NODE_NUM1], the coordinates of
            //    the nodes that make up the linear mesh.
            //
            //    Input, int EDGE_DATA[5*(6*TETRA_NUM)], edge data.
            //
            //    Output, int TETRA_NODE2[10*TETRA_NUM], the indices of the nodes 
            //    in the quadratic mesh.
            //
            //    Input, int NODE_NUM2, the number of nodes for the quadratic mesh.
            //
            //    Output, double NODE_XYZ2[3*NODE_NUM2], the coordinates of
            //    the nodes that make up the quadratic mesh.
            //
        {
            int dim_num = 3;
            int edge;
            int i;
            int j;
            int n1;
            int n1_old;
            int n2;
            int n2_old;
            int node;
            int tetra;
            int tetra_order1 = 4;
            int tetra_order2 = 10;
            int v = 0;
            int v1;
            int v2;
            //
            //  Generate the index and coordinates of the new midside nodes, 
            //  and update the tetradehron-node data.
            //
            for (j = 0; j < node_num1; j++)
            {
                for (i = 0; i < dim_num; i++)
                {
                    node_xyz2[i + j * dim_num] = node_xyz1[i + j * dim_num];
                }
            }

            for (j = 0; j < tetra_num; j++)
            {
                for (i = 0; i < tetra_order1; i++)
                {
                    tetra_node2[i + j * tetra_order2] = tetra_node1[i + j * tetra_order1];
                }
            }

            node = node_num1;

            n1_old = -1;
            n2_old = -1;

            for (edge = 0; edge < 6 * tetra_num; edge++)
            {
                //
                //  Read the data defining the edge.
                //
                n1 = edge_data[0 + edge * 5];
                n2 = edge_data[1 + edge * 5];
                //
                //  If this edge is new, create the coordinates and index.
                //
                if (n1 != n1_old || n2 != n2_old)
                {
                    if (node_num2 <= node)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("TET_MESH_ORDER4_TO_ORDER10_COMPUTE - Fatal error!");
                        Console.WriteLine("  Node index exceeds NODE_NUM2.");
                        return;
                    }

                    for (i = 0; i < dim_num; i++)
                    {
                        node_xyz2[i + node * dim_num] =
                            (node_xyz2[i + (n1 - 1) * dim_num] + node_xyz2[i + (n2 - 1) * dim_num]) / 2.0;
                    }

                    node = node + 1;
                    n1_old = n1;
                    n2_old = n2;
                }

                //
                //  Assign the node to the tetrahedron.
                //
                v1 = edge_data[2 + edge * 5];
                v2 = edge_data[3 + edge * 5];
                //
                //  Here is where the local ordering of the nodes is effected:
                //
                if (v1 == 1 && v2 == 2)
                {
                    v = 5;
                }
                else if (v1 == 1 && v2 == 3)
                {
                    v = 6;
                }
                else if (v1 == 1 && v2 == 4)
                {
                    v = 7;
                }
                else if (v1 == 2 && v2 == 3)
                {
                    v = 8;
                }
                else if (v1 == 2 && v2 == 4)
                {
                    v = 9;
                }
                else if (v1 == 3 && v2 == 4)
                {
                    v = 10;
                }

                tetra = edge_data[4 + edge * 5];

                tetra_node2[v - 1 + tetra * tetra_order2] = node;
            }

        }

        public static void tet_mesh_order4_to_order10_size(int tetra_num, int[] tetra_node1,
                int node_num1, ref int[] edge_data, ref int node_num2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TET_MESH_ORDER4_TO_ORDER10_SIZE sizes a quadratic tet mesh from a linear one.
            //
            //  Discussion:
            //
            //    A quadratic (10 node) tet mesh can be derived from a linear
            //    (4 node) tet mesh by interpolating nodes at the midpoint of
            //    every edge of the mesh.
            //
            //    The mesh is described indirectly, as the sum of individual
            //    tetrahedrons.  A single physical edge may be a logical edge of
            //    any number of tetrahedrons.  It is important, however, that a
            //    new node be created exactly once for each edge, assigned an index,
            //    and associated with every tetrahedron that shares this edge. 
            //
            //    This routine handles that problem.
            //
            //    The primary amount of work occurs in sorting a list of 6 * TETRA_NUM
            //    data items, one item for every edge of every tetrahedron.  Each
            //    data item records, for a given tetrahedron edge, the global indices
            //    of the two endpoints, the local indices of the two endpoints,
            //    and the index of the tetrahedron.
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
            //    09 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int TETRA_NUM, the number of tetrahedrons in the
            //    linear mesh.
            //
            //    Input, int TETRA_NODE1[4*TETRA_NUM], the indices of the nodes
            //    in the linear mesh.
            //
            //    Input, int NODE_NUM1, the number of nodes for the linear mesh.
            //
            //    Output, int EDGE_DATA[5*(6*TETRA_NUM)], edge data.
            //
            //    Output, int *NODE_NUM2, the number of nodes for the quadratic mesh.
            //
        {
            int a = 0;
            int b = 0;
            int edge;
            int i;
            int j;
            int k;
            int l;
            int n1;
            int n1_old;
            int n2;
            int n2_old;
            int tetra;
            int tetra_order1 = 4;
            //
            //  Step 1.
            //  From the list of nodes for tetrahedron T, of the form: (I,J,K,L)
            //  construct the six edge relations:
            //
            //    (I,J,1,2,T)
            //    (I,K,1,3,T)
            //    (I,L,1,4,T)
            //    (J,K,2,3,T)
            //    (J,L,2,4,T)
            //    (K,L,3,4,T)
            //
            //  In order to make matching easier, we reorder each pair of nodes
            //  into ascending order.
            //
            for (tetra = 0; tetra < tetra_num; tetra++)
            {
                i = tetra_node1[0 + tetra * tetra_order1];
                j = tetra_node1[1 + tetra * tetra_order1];
                k = tetra_node1[2 + tetra * tetra_order1];
                l = tetra_node1[3 + tetra * tetra_order1];

                typeMethods.i4i4_sort_a(i, j, ref a, ref b);

                edge_data[0 + (6 * tetra) * 5] = a;
                edge_data[1 + (6 * tetra) * 5] = b;
                edge_data[2 + (6 * tetra) * 5] = 1;
                edge_data[3 + (6 * tetra) * 5] = 2;
                edge_data[4 + (6 * tetra) * 5] = tetra;

                typeMethods.i4i4_sort_a(i, k, ref a, ref b);

                edge_data[0 + (6 * tetra + 1) * 5] = a;
                edge_data[1 + (6 * tetra + 1) * 5] = b;
                edge_data[2 + (6 * tetra + 1) * 5] = 1;
                edge_data[3 + (6 * tetra + 1) * 5] = 3;
                edge_data[4 + (6 * tetra + 1) * 5] = tetra;

                typeMethods.i4i4_sort_a(i, l, ref a, ref b);

                edge_data[0 + (6 * tetra + 2) * 5] = a;
                edge_data[1 + (6 * tetra + 2) * 5] = b;
                edge_data[2 + (6 * tetra + 2) * 5] = 1;
                edge_data[3 + (6 * tetra + 2) * 5] = 4;
                edge_data[4 + (6 * tetra + 2) * 5] = tetra;

                typeMethods.i4i4_sort_a(j, k, ref a, ref b);

                edge_data[0 + (6 * tetra + 3) * 5] = a;
                edge_data[1 + (6 * tetra + 3) * 5] = b;
                edge_data[2 + (6 * tetra + 3) * 5] = 2;
                edge_data[3 + (6 * tetra + 3) * 5] = 3;
                edge_data[4 + (6 * tetra + 3) * 5] = tetra;

                typeMethods.i4i4_sort_a(j, l, ref a, ref b);

                edge_data[0 + (6 * tetra + 4) * 5] = a;
                edge_data[1 + (6 * tetra + 4) * 5] = b;
                edge_data[2 + (6 * tetra + 4) * 5] = 2;
                edge_data[3 + (6 * tetra + 4) * 5] = 4;
                edge_data[4 + (6 * tetra + 4) * 5] = tetra;

                typeMethods.i4i4_sort_a(k, l, ref a, ref b);

                edge_data[0 + (6 * tetra + 5) * 5] = a;
                edge_data[1 + (6 * tetra + 5) * 5] = b;
                edge_data[2 + (6 * tetra + 5) * 5] = 3;
                edge_data[3 + (6 * tetra + 5) * 5] = 4;
                edge_data[4 + (6 * tetra + 5) * 5] = tetra;
            }

            //
            //  Step 2. Perform an ascending dictionary sort on the neighbor relations.
            //  We only intend to sort on rows 1:2; the routine we call here
            //  sorts on the full column but that won't hurt us.
            //
            //  What we need is to find all cases where tetrahedrons share an edge.
            //  By sorting the columns of the EDGE_DATA array, we will put shared edges
            //  next to each other.
            //
            typeMethods.i4col_sort_a(5, 6 * tetra_num, ref edge_data);
            //
            //  Step 3. All the tetrahedrons which share an edge show up as consecutive
            //  columns with identical first two entries.  Figure out how many new
            //  nodes there are, and allocate space for their coordinates.
            //
            node_num2 = node_num1;

            n1_old = -1;
            n2_old = -1;

            for (edge = 0; edge < 6 * tetra_num; edge++)
            {
                n1 = edge_data[0 + edge * 5];
                n2 = edge_data[1 + edge * 5];
                if (n1 != n1_old || n2 != n2_old)
                {
                    node_num2 = node_num2 + 1;
                    n1_old = n1;
                    n2_old = n2;
                }
            }

        }

        public static void tet_mesh_order10_adj_count(int node_num, int tet_num,
                int[] tet_node, ref int adj_num, ref int[] adj_row)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TET_MESH_ORDER10_ADJ_COUNT counts the number of nodal adjacencies.
            //
            //  Discussion:
            //
            //    Assuming that the tet mesh is to be used in a finite element
            //    computation, we declare that two distinct nodes are "adjacent" if and
            //    only if they are both included in some tetrahedron.
            //
            //    It is the purpose of this routine to determine the number of
            //    such adjacency relationships.
            //
            //    The initial count gets only the (I,J) relationships, for which
            //    node I is strictly less than node J.  This value is doubled
            //    to account for symmetry.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 March 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NODE_NUM, the number of nodes.
            //
            //    Input, int TET_NUM, the number of tetrahedrons.
            //
            //    Input, int TET_NODE[10*TET_NUM], the indices of the nodes.
            //
            //    Output, int *ADJ_NUM, the total number of adjacency relationships,
            //
            //    Output, int ADJ_ROW[NODE_NUM+1], the ADJ pointer array.
            //
        {
            int i;
            int j;
            int k;
            int l;
            int node;
            int[] pair;
            int pair_num;
            int pair_unique_num;
            //
            //  Each order 10 tetrahedron defines 45 adjacency pairs.
            //
            pair = new int[2 * 45 * tet_num];

            k = 0;
            for (i = 0; i < 9; i++)
            {
                for (j = i + 1; j < 10; j++)
                {
                    for (l = 0; l < tet_num; l++)
                    {
                        pair[0 + (k * tet_num + l) * 2] = tet_node[i + l * 10];
                        pair[1 + (k * tet_num + l) * 2] = tet_node[j + l * 10];
                    }

                    k = k + 1;
                }
            }

            //
            //  Force the nodes of each pair to be listed in ascending order.
            //
            pair_num = 45 * tet_num;

            typeMethods.i4col_sort2_a(2, pair_num, ref pair);
            //
            //  Rearrange the columns in ascending order.
            //
            typeMethods.i4col_sort_a(2, pair_num, ref pair);
            //
            //  Get the number of unique columns.
            //
            pair_unique_num = typeMethods.i4col_sorted_unique_count(2, pair_num, pair);
            //
            //  The number of adjacencies is TWICE this value, plus the number of nodes.
            //
            adj_num = 2 * pair_unique_num;
            //
            //  Now set up the ADJ_ROW counts.
            //
            for (node = 0; node < node_num; node++)
            {
                adj_row[node] = 0;
            }

            for (k = 0; k < pair_num; k++)
            {
                if (0 < k)
                {
                    if (pair[0 + (k - 1) * 2] == pair[0 + k * 2] &&
                        pair[1 + (k - 1) * 2] == pair[1 + k * 2])
                    {
                        continue;
                    }
                }

                i = pair[0 + k * 2];
                j = pair[1 + k * 2];

                adj_row[i - 1] = adj_row[i - 1] + 1;
                adj_row[j - 1] = adj_row[j - 1] + 1;
            }

            //
            //  We used ADJ_ROW to count the number of entries in each row.
            //  Convert it to pointers into the ADJ array.
            //
            for (node = node_num - 1; 0 <= node; node--)
            {
                adj_row[node] = adj_row[node + 1];
            }

            adj_row[0] = 1;
            for (node = 1; node <= node_num; node++)
            {
                adj_row[node] = adj_row[node - 1] + adj_row[i];
            }
        }

        public static int[] tet_mesh_order10_adj_set(int node_num, int tet_num,
                int[] tet_node, ref int adj_num, ref int[] adj_row)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TET_MESH_ORDER10_ADJ_SET sets the nodal adjacency matrix.
            //
            //  Discussion:
            //
            //    A compressed format is used for the nodal adjacency matrix.
            //
            //    It is assumed that we know ADJ_NUM, the number of adjacency entries
            //    and the ADJ_ROW array, which keeps track of the list of slots
            //    in ADJ where we can store adjacency information for each row.
            //
            //    We essentially repeat the work of TET_MESH_ORDER4_ADJ_COUNT, but
            //    now we have a place to store the adjacency information.
            //
            //    A copy of the ADJ_ROW array is useful, as we can use it to keep track
            //    of the next available entry in ADJ for adjacencies associated with
            //    a given row.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 March 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NODE_NUM, the number of nodes.
            //
            //    Input, int TET_NUM, the number of tetrahedrons.
            //
            //    Input, int TET_NODE[10*TET_NUM], the indices of the nodes.
            //
            //    Input, int ADJ_NUM, the total number of adjacency relationships,
            //
            //    Input, int ADJ_ROW[NODE_NUM+1], the ADJ pointer array.
            //
            //    Output, int TET_MESH_ORDER4_ADJ_SET[ADJ_NUM], 
            //    the adjacency information.
            //
        {
            int[] adj;
            int[] adj_row_copy;
            int i;
            int j;
            int k;
            int l;
            int node;
            int[] pair;
            int pair_num;
            //
            //  Each order 10 tetrahedron defines 45 adjacency pairs.
            //
            pair = new int[2 * 45 * tet_num];

            k = 0;
            for (i = 0; i < 9; i++)
            {
                for (j = i + 1; j < 10; j++)
                {
                    for (l = 0; l < tet_num; l++)
                    {
                        pair[0 + (k * tet_num + l) * 2] = tet_node[i + l * 10];
                        pair[1 + (k * tet_num + l) * 2] = tet_node[j + l * 10];
                    }

                    k = k + 1;
                }
            }

            //
            //  Force the nodes of each pair to be listed in ascending order.
            //
            pair_num = 45 * tet_num;

            typeMethods.i4col_sort2_a(2, pair_num, ref pair);
            //
            //  Rearrange the columns in ascending order.
            //
            typeMethods.i4col_sort_a(2, pair_num, ref pair);
            //
            //  Mark all entries of ADJ so we will know later if we missed one.
            //
            adj = new int[adj_num];

            for (i = 0; i < adj_num; i++)
            {
                adj[i] = -1;
            }

            //
            //  Copy the ADJ_ROW array and use it to keep track of the next
            //  free entry for each row.
            //
            adj_row_copy = new int[node_num];

            for (node = 0; node < node_num; node++)
            {
                adj_row_copy[node] = adj_row[node];
            }

            //
            //  Now set up the ADJ_ROW counts.
            //
            for (k = 0; k < pair_num; k++)
            {
                if (0 < k)
                {
                    if (pair[0 + (k - 1) * 2] == pair[0 + k * 2] &&
                        pair[1 + (k - 1) * 2] == pair[1 + k * 2])
                    {
                        continue;
                    }
                }

                i = pair[0 + k * 2];
                j = pair[1 + k * 2];

                adj[adj_row_copy[i]] = j;
                adj_row_copy[i] = adj_row_copy[i] + 1;
                adj[adj_row_copy[j]] = i;
                adj_row_copy[j] = adj_row_copy[j] + 1;
            }

            return adj;
        }

        public static void tet_mesh_order10_example_set(int node_num, int tetra_num,
                ref double[] node_xyz, ref int[] tetra_node)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TET_MESH_ORDER10_EXAMPLE_SET sets an example quadratic tet mesh.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    03 August 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NODE_NUM, the number of nodes.
            //
            //    Input, int TETRA_NUM, the number of tetrahedrons.
            //
            //    Output, double NODE_XYZ[3*NODE_NUM], the node coordinates.
            //
            //    Output, int TETRA_NODE[10*TETRA_NUM], the nodes forming each tet.
            //
        {
            int i;
            int j;
            double[] node_xyz_save =
                {
                    0.0, 0.0, 0.0,
                    0.0, 0.0, 1.0,
                    0.0, 1.0, 0.0,
                    0.0, 1.0, 1.0,
                    1.0, 0.0, 0.0,
                    1.0, 0.0, 1.0,
                    1.0, 1.0, 0.0,
                    1.0, 1.0, 1.0,
                    0.0, 0.0, 0.5,
                    0.0, 0.5, 0.0,
                    0.0, 0.5, 0.5,
                    0.5, 0.0, 0.0,
                    0.0, 0.5, 1.0,
                    0.5, 0.0, 0.5,
                    0.5, 0.0, 1.0,
                    0.0, 1.0, 0.5,
                    0.5, 0.5, 0.0,
                    0.5, 1.0, 0.0,
                    0.5, 0.5, 0.5,
                    0.5, 0.5, 1.0,
                    0.5, 1.0, 0.5,
                    0.5, 1.0, 1.0,
                    1.0, 0.0, 0.5,
                    1.0, 0.5, 0.0,
                    1.0, 0.5, 0.5,
                    1.0, 0.5, 1.0,
                    1.0, 1.0, 0.5
                }
                ;
            int[] tetra_node_save =
                {
                    4, 3, 5, 1, 16, 19, 17, 11, 10, 12,
                    4, 2, 5, 1, 13, 19, 14, 11, 9, 12,
                    4, 7, 3, 5, 21, 16, 18, 19, 24, 17,
                    4, 7, 8, 5, 21, 22, 27, 19, 24, 25,
                    4, 6, 2, 5, 20, 13, 15, 19, 23, 14,
                    4, 6, 8, 5, 20, 22, 26, 19, 23, 25
                }
                ;

            for (j = 0; j < node_num; j++)
            {
                for (i = 0; i < 3; i++)
                {
                    node_xyz[i + j * 3] = node_xyz_save[i + j * 3];
                }
            }

            for (j = 0; j < tetra_num; j++)
            {
                for (i = 0; i < 10; i++)
                {
                    tetra_node[i + j * 10] = tetra_node_save[i + j * 10] - 1;
                }
            }
        }

        public static void tet_mesh_order10_example_size(ref int node_num, ref int tetra_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TET_MESH_ORDER10_EXAMPLE_SIZE sizes an example quadratic tet mesh.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    03 August 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, int *NODE_NUM, the number of nodes.
            //
            //    Output, int *TETRA_NUM, the number of tetrahedrons.
            //
        {
            node_num = 27;
            tetra_num = 6;

        }

        public static void tet_mesh_order10_to_order4_compute(int tetra_num1, int[] tetra_node1,
                int tetra_num2, ref int[] tetra_node2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TET_MESH_ORDER10_TO_ORDER4_COMPUTE linearizes a quadratic tet mesh.
            //
            //  Discussion:
            //
            //    A quadratic tet mesh is assumed to consist of 10-node
            //    tetrahedrons.
            //
            //    This routine rearranges the information so as to define a 4-node
            //    tet mesh.
            //
            //    The same nodes are used, but there are 8 times as many
            //    tetrahedrons.
            //
            //    The node ordering for the quadratic tetrahedron is somewhat
            //    arbitrary.  In the current scheme, the vertices are listed
            //    first, followed by the 6 midside nodes.  Each midside node
            //    may be identified by the two vertices that bracket it.  Thus,
            //    the node ordering may be suggested by:
            //
            //      1  2  3  4 (1+2) (1+3) (1+4) (2+3) (2+4) (3+4)
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    22 January 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Anwei Liu, Barry Joe,
            //    Quality Local Refinement of Tetrahedral Meshes Based
            //    on 8-Subtetrahedron Subdivision,
            //    Mathematics of Computation,
            //    Volume 65, Number 215, July 1996, pages 1183-1200.
            //
            //  Parameters:
            //
            //    Input, int TETRA_NUM1, the number of tetrahedrons in the quadratic
            //    tet mesh.
            //
            //    Input, int TETRA_NODE1[10*TETRA_NUM1], the indices of the nodes
            //    that made up the quadratic mesh.
            //
            //    Input, int TETRA_NUM2, the number of tetrahedrons in the linear
            //    tet mesh.  TETRA_NUM2 = 8 * TETRA_NUM1.
            //
            //    Output, int TETRA_NODE2[4*TETRA_NUM2], the indices of the nodes
            //    that make up the linear mesh.
            //
        {
            int n1;
            int n2;
            int n3;
            int n4;
            int n5;
            int n6;
            int n7;
            int n8;
            int n9;
            int nx;
            int tetra1;
            int tetra2;

            tetra2 = 0;

            for (tetra1 = 0; tetra1 < tetra_num1; tetra1++)
            {
                n1 = tetra_node1[0 + tetra1 * 10];
                n2 = tetra_node1[1 + tetra1 * 10];
                n3 = tetra_node1[2 + tetra1 * 10];
                n4 = tetra_node1[3 + tetra1 * 10];
                n5 = tetra_node1[4 + tetra1 * 10];
                n6 = tetra_node1[5 + tetra1 * 10];
                n7 = tetra_node1[6 + tetra1 * 10];
                n8 = tetra_node1[7 + tetra1 * 10];
                n9 = tetra_node1[8 + tetra1 * 10];
                nx = tetra_node1[9 + tetra1 * 10];

                tetra_node2[0 + tetra2 * 4] = n1;
                tetra_node2[1 + tetra2 * 4] = n5;
                tetra_node2[2 + tetra2 * 4] = n6;
                tetra_node2[3 + tetra2 * 4] = n7;
                tetra2 = tetra2 + 1;

                tetra_node2[0 + tetra2 * 4] = n2;
                tetra_node2[1 + tetra2 * 4] = n5;
                tetra_node2[2 + tetra2 * 4] = n8;
                tetra_node2[3 + tetra2 * 4] = n9;
                tetra2 = tetra2 + 1;

                tetra_node2[0 + tetra2 * 4] = n3;
                tetra_node2[1 + tetra2 * 4] = n6;
                tetra_node2[2 + tetra2 * 4] = n8;
                tetra_node2[3 + tetra2 * 4] = n9;
                tetra2 = tetra2 + 1;

                tetra_node2[0 + tetra2 * 4] = n4;
                tetra_node2[1 + tetra2 * 4] = n7;
                tetra_node2[2 + tetra2 * 4] = n9;
                tetra_node2[3 + tetra2 * 4] = nx;
                tetra2 = tetra2 + 1;

                tetra_node2[0 + tetra2 * 4] = n5;
                tetra_node2[1 + tetra2 * 4] = n6;
                tetra_node2[2 + tetra2 * 4] = n7;
                tetra_node2[3 + tetra2 * 4] = n9;
                tetra2 = tetra2 + 1;

                tetra_node2[0 + tetra2 * 4] = n5;
                tetra_node2[1 + tetra2 * 4] = n6;
                tetra_node2[2 + tetra2 * 4] = n8;
                tetra_node2[3 + tetra2 * 4] = n9;
                tetra2 = tetra2 + 1;

                tetra_node2[0 + tetra2 * 4] = n6;
                tetra_node2[1 + tetra2 * 4] = n7;
                tetra_node2[2 + tetra2 * 4] = n9;
                tetra_node2[3 + tetra2 * 4] = nx;
                tetra2 = tetra2 + 1;

                tetra_node2[0 + tetra2 * 4] = n6;
                tetra_node2[1 + tetra2 * 4] = n8;
                tetra_node2[2 + tetra2 * 4] = n9;
                tetra_node2[3 + tetra2 * 4] = nx;
                tetra2 = tetra2 + 1;
            }
        }

        public static void tet_mesh_order10_to_order4_size(int node_num1, int tetra_num1,
                ref int node_num2, ref int tetra_num2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TET_MESH_ORDER10_TO_ORDER4_SIZE sizes a linear tet mesh from a quadratic one.
            //
            //  Discussion:
            //
            //    A linear (4 node) tet mesh can be derived from a quadratic
            //    (10 node) tet mesh using the same set of nodes, but reassigning
            //    the nodes of each quadratic tet among 8 linear subtets.
            //
            //    This routine returns the number of nodes and tetrahedra in the
            //    linear mesh.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    06 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Anwei Liu, Barry Joe,
            //    Quality Local Refinement of Tetrahedral Meshes Based
            //    on 8-Subtetrahedron Subdivision,
            //    Mathematics of Computation,
            //    Volume 65, Number 215, July 1996, pages 1183-1200.
            //
            //  Parameters:
            //
            //    Input, int NODE_NUM1, the number of nodes in the quadratic mesh.
            //
            //    Input, int TETRA_NUM1, the number of tetrahedrons in the
            //    quadratic mesh.
            //
            //    Output, int *NODE_NUM2, the number of nodes for the linear mesh.
            //
            //    Output, int *TETRA_NUM2, the number of tetrahedrons in the
            //    linear mesh.
            //
        {
            node_num2 = node_num1;
            tetra_num2 = 8 * tetra_num1;
        }

        public static void tet_mesh_quad(int node_num, double[] node_xyz, int tetra_order,
                int tetra_num, int[] tetra_node,
                Func<int, double[], double[], double[]> quad_fun,
                int quad_num, double[] quad_xyz, double[] quad_w, ref double quad_value,
                ref double region_volume)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TET_MESH_QUAD approximates an integral over a tet mesh.
            //
            //  Discussion:
            //
            //    The routine will accept tetrahedral meshes of order higher than 4.
            //    However, only the first four nodes (the vertices) of each
            //    tetrahedron will be used.  This will still produce correct results
            //    for higher order tet meshes, as long as the sides of each
            //    tetrahedron are flat (linear).
            //
            //    We assume that the vertices of each tetrahedron are listed first
            //    in the description of higher order tetrahedrons.
            //
            //    The approximation of the integral is made using a quadrature rule 
            //    defined on the unit tetrahedron, and supplied by the user.  
            //
            //    The user also supplies the name of a subroutine, here called "QUAD_FUN", 
            //    which evaluates the integrand at a set of points.  The form is:
            //
            //      void quad_fun ( int n, double xyz_vec[3*n], double f_vec[n] )
            //
            //    and it returns in each entry F_VEC(1:N), the value of the integrand
            //    at XYZ_VEC(1:3,1:N).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    06 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NODE_NUM, the number of nodes in the tet mesh.
            //
            //    Input, double NODE_XYZ[3*NODE_NUM], the coordinates of the nodes.
            //
            //    Input, int TETRA_ORDER, the order of tetrahedrons in the tet mesh.
            //
            //    Input, int TETRA_NUM, the number of tetrahedrons in the tet mesh.
            //
            //    Input, int TETRA_NODE[TETRA_ORDER*TETRA_NUM], indices of the nodes.
            //
            //    Input, void QUAD_FUN ( int N, double XYZ_VEC[3*N], F_VEC[N] ), the name 
            //    of the routine that evaluates the integrand.
            //
            //    Input, int QUAD_NUM, the order of the quadrature rule.
            //
            //    Input, double QUAD_XYZ[3*QUAD_NUM], the abscissas of the 
            //    quadrature rule, in the unit tetrahedron.
            //
            //    Input, double QUAD_W[QUAD_NUM], the weights of the 
            //    quadrature rule.
            //
            //    Output, double *QUAD_VALUE, the estimate of the integral
            //    of F(X,Y) over the region covered by the tet mesh.
            //
            //    Output, double *REGION_VOLUME, the volume of the region.
            //
        {
            int i;
            int j;
            int quad;
            double[] quad_f = new double[quad_num];
            double[] quad2_xyz = new double[3 * quad_num];
            double temp;
            int tet;
            double tetra_volume;
            double[] tetra_xyz = new double[3 * 4];

            quad_value = 0.0;
            region_volume = 0.0;

            for (tet = 0; tet < tetra_num; tet++)
            {
                for (j = 0; j < 4; j++)
                {
                    for (i = 0; i < 3; i++)
                    {
                        tetra_xyz[i + j * 3] = node_xyz[i + (tetra_node[j + tet * 4] - 1) * 3];
                    }
                }

                tetra_volume = Tetrahedron.tetrahedron_volume(tetra_xyz);

                Tetrahedron.tetrahedron_order4_reference_to_physical(tetra_xyz, quad_num,
                    quad_xyz, ref quad2_xyz);

                quad_f = quad_fun(quad_num, quad2_xyz, quad_f);

                temp = 0.0;
                for (quad = 0; quad < quad_num; quad++)
                {
                    temp = temp + quad_w[quad] * quad_f[quad];
                }

                quad_value = quad_value + tetra_volume * temp;

                region_volume = region_volume + tetra_volume;
            }

        }

        public static void tet_mesh_quality1(int node_num, double[] node_xyz,
                int tetra_order, int tetra_num, int[] tetra_node, ref double value_min,
                ref double value_mean, ref double value_max, ref double value_var)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TET_MESH_QUALITY1 returns a tet mesh quality factor.
            //
            //  Discussion:
            //
            //    The tet mesh quality measure is the minimum of the 
            //    corresponding tetrahedron quality measure, over all tetrahedrons in the 
            //    tet mesh.
            //
            //    This routine is designed for a 4-node tet mesh.  It can handle a 10-node
            //    tet mesh, but it simply ignores the extra nodes.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    27 October 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NODE_NUM, the number of nodes.
            //
            //    Input, double NODE_XYZ[3*NODE_NUM], the coordinates of the nodes.
            //
            //    Input, int TETRA_ORDER, the order of the tetrahedrons.
            //
            //    Input, int TETRA_NUM, the number of tetrahedrons.
            //
            //    Input, int TETRA_NODE[TETRA_ORDER*TETRA_NUM], the indices of the nodes.
            //
            //    Output, double *VALUE_MIN, *VALUE_MEAN, *VALUE_MAX, *VALUE_VAR,
            //    the minimum, mean, maximum and variance of the quality measure.
            //
        {
            int DIM_NUM = 3;

            int i;
            int j;
            int node;
            int tetra;
            double[] tetrahedron = new double[DIM_NUM * 4];
            double[] tetrahedron_quality;

            tetrahedron_quality = new double[tetra_num];

            for (tetra = 0; tetra < tetra_num; tetra++)
            {
                for (j = 0; j < 4; j++)
                {
                    node = tetra_node[j + tetra * tetra_order];
                    for (i = 0; i < DIM_NUM; i++)
                    {
                        tetrahedron[i + j * DIM_NUM] = node_xyz[i + (node - 1) * DIM_NUM];
                    }
                }

                tetrahedron_quality[tetra] = Properties.tetrahedron_quality1_3d(tetrahedron);
            }

            value_max = typeMethods.r8vec_max(tetra_num, tetrahedron_quality);
            value_min = typeMethods.r8vec_min(tetra_num, tetrahedron_quality);
            value_mean = typeMethods.r8vec_mean(tetra_num, tetrahedron_quality);
            value_var = typeMethods.r8vec_variance(tetra_num, tetrahedron_quality);

        }

        public static void tet_mesh_quality2(int node_num, double[] node_xyz, int tetra_order,
                int tetra_num, int[] tetra_node, ref double value_min, ref double value_mean,
                ref double value_max, ref double value_var)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TET_MESH_QUALITY2 returns a tet mesh quality factor.
            //
            //  Discussion:
            //
            //    The tet mesh quality measure is the minimum of the 
            //    corresponding tetrahedron quality measure, over all tetrahedrons in the 
            //    tet mesh.
            //
            //    This routine is designed for a 4-node tet mesh.  It can handle a 10-node
            //    tet mesh, but it simply ignores the extra nodes.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    27 October 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NODE_NUM, the number of nodes.
            //
            //    Input, double NODE_XYZ[3*NODE_NUM], the coordinates of the nodes.
            //
            //    Input, int TETRA_ORDER, the order of the tetrahedrons.
            //
            //    Input, int TETRA_NUM, the number of tetrahedrons.
            //
            //    Input, int TETRA_NODE[TETRA_ORDER*TETRA_NUM], the indices of the nodes.
            //
            //    Output, double *VALUE_MIN, *VALUE_MEAN, *VALUE_MAX, *VALUE_VAR,
            //    the minimum, mean, maximum and variance of the quality measure.
            //
        {
            int DIM_NUM = 3;

            int i;
            int j;
            int node;
            int tetra;
            double[] tetrahedron = new double[DIM_NUM * 4];
            double[] tetrahedron_quality;

            tetrahedron_quality = new double[tetra_num];

            for (tetra = 0; tetra < tetra_num; tetra++)
            {
                for (j = 0; j < 4; j++)
                {
                    node = tetra_node[j + tetra * tetra_order];
                    for (i = 0; i < DIM_NUM; i++)
                    {
                        tetrahedron[i + j * DIM_NUM] = node_xyz[i + (node - 1) * DIM_NUM];
                    }
                }

                tetrahedron_quality[tetra] = Properties.tetrahedron_quality2_3d(tetrahedron);
            }

            value_max = typeMethods.r8vec_max(tetra_num, tetrahedron_quality);
            value_min = typeMethods.r8vec_min(tetra_num, tetrahedron_quality);
            value_mean = typeMethods.r8vec_mean(tetra_num, tetrahedron_quality);
            value_var = typeMethods.r8vec_variance(tetra_num, tetrahedron_quality);
        }

        public static void tet_mesh_quality3(int node_num, double[] node_xyz, int tetra_order,
                int tetra_num, int[] tetra_node, ref double value_min, ref double value_mean,
                ref double value_max, ref double value_var)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TET_MESH_QUALITY3 returns a tet mesh quality factor.
            //
            //  Discussion:
            //
            //    The tet mesh quality measure is the minimum of the 
            //    corresponding tetrahedron quality measure, over all tetrahedrons in the 
            //    tet mesh.
            //
            //    This routine is designed for a 4-node tet mesh.  It can handle a 10-node
            //    tet mesh, but it simply ignores the extra nodes.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    27 October 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NODE_NUM, the number of nodes.
            //
            //    Input, double NODE_XYZ[3*NODE_NUM], the coordinates of the nodes.
            //
            //    Input, int TETRA_ORDER, the order of the tetrahedrons.
            //
            //    Input, int TETRA_NUM, the number of tetrahedrons.
            //
            //    Input, int TETRA_NODE[TETRA_ORDER*TETRA_NUM], the indices of the nodes.
            //
            //    Output, double *VALUE_MIN, *VALUE_MEAN, *VALUE_MAX, *VALUE_VAR,
            //    the minimum, mean, maximum and variance of the quality measure.
            //
        {
            int DIM_NUM = 3;

            int i;
            int j;
            int node;
            int tetra;
            double[] tetrahedron = new double[DIM_NUM * 4];
            double[] tetrahedron_quality;

            tetrahedron_quality = new double[tetra_num];

            for (tetra = 0; tetra < tetra_num; tetra++)
            {
                for (j = 0; j < 4; j++)
                {
                    node = tetra_node[j + tetra * tetra_order];
                    for (i = 0; i < DIM_NUM; i++)
                    {
                        tetrahedron[i + j * DIM_NUM] = node_xyz[i + (node - 1) * DIM_NUM];
                    }
                }

                tetrahedron_quality[tetra] = Properties.tetrahedron_quality3_3d(tetrahedron);
            }

            value_max = typeMethods.r8vec_max(tetra_num, tetrahedron_quality);
            value_min = typeMethods.r8vec_min(tetra_num, tetrahedron_quality);
            value_mean = typeMethods.r8vec_mean(tetra_num, tetrahedron_quality);
            value_var = typeMethods.r8vec_variance(tetra_num, tetrahedron_quality);
        }

        public static void tet_mesh_quality4(int node_num, double[] node_xyz, int tetra_order,
                int tetra_num, int[] tetra_node, ref double value_min, ref double value_mean,
                ref double value_max, ref double value_var)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TET_MESH_QUALITY4 returns a tet mesh quality factor.
            //
            //  Discussion:
            //
            //    The tet mesh quality measure is the minimum of the 
            //    corresponding tetrahedron quality measure, over all tetrahedrons in the 
            //    tet mesh.
            //
            //    This routine is designed for a 4-node tet mesh.  It can handle a 10-node
            //    tet mesh, but it simply ignores the extra nodes.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    27 October 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NODE_NUM, the number of nodes.
            //
            //    Input, double NODE_XYZ[3*NODE_NUM], the coordinates of the nodes.
            //
            //    Input, int TETRA_ORDER, the order of the tetrahedrons.
            //
            //    Input, int TETRA_NUM, the number of tetrahedrons.
            //
            //    Input, int TETRA_NODE[TETRA_ORDER*TETRA_NUM], the indices of the nodes.
            //
            //    Output, double *VALUE_MIN, *VALUE_MEAN, *VALUE_MAX, *VALUE_VAR,
            //    the minimum, mean, maximum and variance of the quality measure.
            //
        {
            int DIM_NUM = 3;

            int i;
            int j;
            int node;
            int tetra;
            double[] tetrahedron = new double[DIM_NUM * 4];
            double[] tetrahedron_quality;

            tetrahedron_quality = new double[tetra_num];

            for (tetra = 0; tetra < tetra_num; tetra++)
            {
                for (j = 0; j < 4; j++)
                {
                    node = tetra_node[j + tetra * tetra_order];
                    for (i = 0; i < DIM_NUM; i++)
                    {
                        tetrahedron[i + j * DIM_NUM] = node_xyz[i + (node - 1) * DIM_NUM];
                    }
                }

                tetrahedron_quality[tetra] = Properties.tetrahedron_quality4_3d(tetrahedron);
            }

            value_max = typeMethods.r8vec_max(tetra_num, tetrahedron_quality);
            value_min = typeMethods.r8vec_min(tetra_num, tetrahedron_quality);
            value_mean = typeMethods.r8vec_mean(tetra_num, tetrahedron_quality);
            value_var = typeMethods.r8vec_variance(tetra_num, tetrahedron_quality);
        }

        public static void tet_mesh_quality5(int node_num, double[] node_xyz, int tetra_order,
                int tetra_num, int[] tetra_node, ref double value_min, ref double value_mean,
                ref double value_max, ref double value_var)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TET_MESH_QUALITY5 returns a tet mesh quality factor.
            //
            //  Discussion:
            //
            //    The tet mesh quality measure is the ratio of the minimum
            //    tetrahedron volume to the maximum tetrahedron volume.
            //
            //    This routine is designed for a 4-node tet mesh.  It can handle a 10-node
            //    tet mesh, but it simply ignores the extra nodes.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    27 October 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NODE_NUM, the number of nodes.
            //
            //    Input, double NODE_XYZ[3*NODE_NUM], the coordinates of the nodes.
            //
            //    Input, int TETRA_ORDER, the order of the tetrahedrons.
            //
            //    Input, int TETRA_NUM, the number of tetrahedrons.
            //
            //    Input, int TETRA_NODE[TETRA_ORDER*TETRA_NUM], the indices of the nodes.
            //
            //    Output, double *VALUE_MIN, *VALUE_MEAN, *VALUE_MAX, *VALUE_VAR,
            //    the minimum, mean, maximum and variance of the quality measure.
            //
        {
            int DIM_NUM = 3;

            int i;
            int j;
            int node;
            int tetra;
            double[] tetrahedron = new double[DIM_NUM * 4];
            double[] tetrahedron_quality;
            double volume_max;

            tetrahedron_quality = new double[tetra_num];

            for (tetra = 0; tetra < tetra_num; tetra++)
            {
                for (j = 0; j < 4; j++)
                {
                    node = tetra_node[j + tetra * tetra_order];
                    for (i = 0; i < DIM_NUM; i++)
                    {
                        tetrahedron[i + j * DIM_NUM] = node_xyz[i + (node - 1) * DIM_NUM];
                    }
                }

                tetrahedron_quality[tetra] = Tetrahedron.tetrahedron_volume(tetrahedron);
            }

            volume_max = typeMethods.r8vec_max(tetra_num, tetrahedron_quality);

            for (tetra = 0; tetra < tetra_num; tetra++)
            {
                tetrahedron_quality[tetra] = tetrahedron_quality[tetra] / volume_max;
            }

            value_max = typeMethods.r8vec_max(tetra_num, tetrahedron_quality);
            value_min = typeMethods.r8vec_min(tetra_num, tetrahedron_quality);
            value_mean = typeMethods.r8vec_mean(tetra_num, tetrahedron_quality);
            value_var = typeMethods.r8vec_variance(tetra_num, tetrahedron_quality);
        }
        
        public static int tet_mesh_search_delaunay(int node_num, double[] node_xyz, int tet_order,
                int tet_num, int[] tet_node, int[] tet_neighbor, double[] p, ref int face,
                ref int step_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TET_MESH_SEARCH_DELAUNAY searches a Delaunay tet mesh for a point.
            //
            //  Discussion:
            //
            //    The algorithm "walks" from one tetrahedron to its neighboring tetrahedron,
            //    and so on, until a tetrahedron is found containing point P, or P is found
            //    to be outside the convex hull.
            //
            //    The algorithm computes the barycentric coordinates of the point with
            //    respect to the current tetrahedron.  If all 4 quantities are positive,
            //    the point is contained in the tetrahedron.  If the I-th coordinate is
            //    negative, then P lies on the far side of edge I, which is opposite
            //    from vertex I.  This gives a hint as to where to search next.
            //
            //    For a Delaunay tet mesh, the search is guaranteed to terminate.
            //    For other meshes, a continue may occur.
            //
            //    Note the surprising fact that, even for a Delaunay tet mesh of
            //    a set of nodes, the nearest node to P need not be one of the
            //    vertices of the tetrahedron containing P.
            //
            //    The code can be called for tet meshes of any order, but only
            //    the first 4 nodes in each tetrahedron are considered.  Thus, if
            //    higher order tetrahedrons are used, and the extra nodes are intended
            //    to give the tetrahedron a polygonal shape, these will have no effect,
            //    and the results obtained here might be misleading.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    19 August 2009
            //
            //  Author:
            //
            //    John Burkardt.
            //
            //  Reference:
            //
            //    Barry Joe,
            //    GEOMPACK - a software package for the generation of meshes
            //    using geometric algorithms,
            //    Advances in Engineering Software,
            //    Volume 13, pages 325-331, 1991.
            //
            //  Parameters:
            //
            //    Input, int NODE_NUM, the number of nodes.
            //
            //    Input, double NODE_XYZ[3*NODE_NUM], the coordinates of 
            //    the nodes.
            //
            //    Input, int TET_ORDER, the order of the tetrahedrons.
            //
            //    Input, int TET_NUM, the number of tetrahedrons.
            //
            //    Input, int TET_NODE[TET_ORDER*TET_NUM],
            //    the nodes that make up each tetrahedron.
            //
            //    Input, int TET_NEIGHBOR[4*TET_NUM], the 
            //    tetrahedron neighbor list.
            //
            //    Input, double P[3], the coordinates of a point.
            //
            //    Output, int *FACE, indicates the position of the point P in
            //    face TET_INDEX:
            //    0, the interior or boundary of the tetrahedron;
            //    -1, outside the convex hull of the tet mesh, past face 1;
            //    -2, outside the convex hull of the tet mesh, past face 2;
            //    -3, outside the convex hull of the tet mesh, past face 3.
            //    -4, outside the convex hull of the tet mesh, past face 4.
            //
            //    Output, int *STEP_NUM, the number of steps taken.
            //
            //    Output, int TET_MESH_SEARCH_DELAUNAY, the index of the tetrahedron 
            //    where the search ended.  If a cycle occurred, then -1 is returned.
            //
        {
            double[] alpha;
            int i;
            int j;
            int k;
            int tet_index;
            double[] tet_xyz = new double[3 * 4];
            int tet_index_save = -1;
            //
            //  If possible, start with the previous successful value of TET_INDEX.
            //
            if (tet_index_save < 1 || tet_num < tet_index_save)
            {
                tet_index = (tet_num + 1) / 2;
            }
            else
            {
                tet_index = tet_index_save;
            }

            step_num = -1;
            face = 0;

            for (;;)
            {
                step_num = step_num + 1;

                if (tet_num < step_num)
                {
                    Console.WriteLine("");
                    Console.WriteLine("TET_MESH_SEARCH_DELAUNAY - Fatal error!");
                    Console.WriteLine("  The algorithm seems to be cycling.");
                    tet_index = -1;
                    face = -1;
                    return 1;
                }

                for (j = 0; j < 4; j++)
                {
                    k = tet_node[j + tet_index * 4];
                    for (i = 0; i < 3; i++)
                    {
                        tet_xyz[i + j * 3] = node_xyz[i + k * 3];
                    }
                }

                alpha = Tetrahedron.tetrahedron_barycentric(tet_xyz, p);
                //
                //  If the barycentric coordinates are all positive, then the point
                //  is inside the tetrahedron and we're done.
                //
                if (0.0 <= alpha[0] && 0.0 <= alpha[1] && 0.0 <= alpha[2] && 0.0 <= alpha[3])
                {
                    break;
                }

                //
                //  At least one barycentric coordinate is negative.
                //
                //  If there is a negative barycentric coordinate for which there exists an
                //  opposing tetrahedron neighbor closer to the point, move to that tetrahedron.
                //
                if (alpha[0] < 0.0 && 0 < tet_neighbor[0 + tet_index * 4])
                {
                    tet_index = tet_neighbor[0 + tet_index * 4];
                    continue;
                }
                else if (alpha[1] < 0.0 && 0 < tet_neighbor[1 + tet_index * 4])
                {
                    tet_index = tet_neighbor[1 + tet_index * 4];
                    continue;
                }
                else if (alpha[2] < 0.0 && 0 < tet_neighbor[2 + tet_index * 4])
                {
                    tet_index = tet_neighbor[2 + tet_index * 4];
                    continue;
                }
                else if (alpha[3] < 0.0 && 0 < tet_neighbor[3 + tet_index * 4])
                {
                    tet_index = tet_neighbor[3 + tet_index * 4];
                    continue;
                }

                //
                //  All negative barycentric coordinates correspond to vertices opposite
                //  faces on the convex hull.
                //
                //  Note the face and exit.
                //
                if (alpha[0] < 0.0)
                {
                    face = -1;
                    break;
                }
                else if (alpha[1] < 0.0)
                {
                    face = -2;
                    break;
                }
                else if (alpha[2] < 0.0)
                {
                    face = -3;
                    break;
                }
                else if (alpha[3] < 0.0)
                {
                    face = -4;
                    break;
                }
            }

            tet_index_save = tet_index;

            return tet_index;
        }

        public static int tet_mesh_search_naive(int node_num, double[] node_xyz,
                int tet_order, int tet_num, int[] tet_node, double[] p, ref int step_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TET_MESH_SEARCH_NAIVE naively searches a tet mesh.
            //
            //  Discussion:
            //
            //    The algorithm simply checks each tetrahedron to see if point P is
            //    contained in it.  
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    19 August 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NODE_NUM, the number of nodes.
            //
            //    Input, double NODE_XYZ[3*NODE_NUM], the coordinates 
            //    of the nodes.
            //
            //    Input, int TET_ORDER, the order of the tetrahedrons.
            //
            //    Input, int TET_NUM, the number of tetrahedrons in
            //    the mesh.
            //
            //    Input, int TET_NODE[TET_ORDER*TET_NUM], 
            //    the nodes that make up each tetrahedron.
            //
            //    Input, double P[3], the coordinates of a point.
            //
            //    Output, int TET_MESH_ORDER4_SEARCH_NAIE, the index of the tetrahedron
            //    where the search ended, or -1 if no tetrahedron was found containing
            //    the point.
            //
            //    Output, int *STEP_NUM, the number of tetrahedrons examined.
        {
            double[] alpha;
            int i;
            int j;
            int tet;
            int tet_index;
            double[] tet_xyz = new double[3 * 4];

            tet_index = -1;
            step_num = 0;

            for (tet = 0; tet < tet_num; tet++)
            {
                for (j = 0; j < 4; j++)
                {
                    for (i = 0; i < 3; i++)
                    {
                        tet_xyz[i + j * 3] = node_xyz[i + tet_node[j + tet * 4] * 3];
                    }
                }

                alpha = Tetrahedron.tetrahedron_barycentric(tet_xyz, p);

                if (typeMethods.r8vec_is_nonnegative(4, alpha))
                {
                    tet_index = tet;
                    step_num = tet;
                    return tet_index;
                }
            }

            return tet_index;
        }
    }
}