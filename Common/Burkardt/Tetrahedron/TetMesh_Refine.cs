using System;
using Burkardt.Types;

namespace Burkardt.TetrahedronNS
{
    public static class TetMesh_Refine
    {
        public static void tet_mesh_order4_refine_compute(int node_num1, int element_num1,
                double[] node_xyz1, int[] element_node1, int node_num2, int element_num2,
                int[] edge_data, ref double[] node_xyz2, ref int[] element_node2)

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
            int element_order = 4;
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

            for (j = 0; j < element_num2; j++)
            {
                for (i = 0; i < element_order; i++)
                {
                    element_node2[i + j * element_order] = -1;
                }
            }

            //
            //  The vertices of the input tetrahedron can be assigned now.
            //
            for (tetra1 = 0; tetra1 < element_num1; tetra1++)
            {
                element_node2[0 + (tetra1 * 8 + 0) * element_order] = element_node1[0 + tetra1 * element_order];
                element_node2[0 + (tetra1 * 8 + 1) * element_order] = element_node1[1 + tetra1 * element_order];
                element_node2[0 + (tetra1 * 8 + 2) * element_order] = element_node1[2 + tetra1 * element_order];
                element_node2[0 + (tetra1 * 8 + 3) * element_order] = element_node1[3 + tetra1 * element_order];
            }

            node = node_num1;

            n1_old = -1;
            n2_old = -1;

            for (edge = 0; edge < 6 * element_num1; edge++)
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
                    element_node2[1 + (tetra1 * 8 + 0) * element_order] = node;
                    element_node2[1 + (tetra1 * 8 + 1) * element_order] = node;
                    element_node2[0 + (tetra1 * 8 + 4) * element_order] = node;
                    element_node2[0 + (tetra1 * 8 + 5) * element_order] = node;
                }
                else if (v1 == 1 && v2 == 3)
                {
                    element_node2[2 + (tetra1 * 8 + 0) * element_order] = node;
                    element_node2[1 + (tetra1 * 8 + 2) * element_order] = node;
                    element_node2[1 + (tetra1 * 8 + 4) * element_order] = node;
                    element_node2[1 + (tetra1 * 8 + 5) * element_order] = node;
                    element_node2[0 + (tetra1 * 8 + 6) * element_order] = node;
                    element_node2[0 + (tetra1 * 8 + 7) * element_order] = node;
                }
                else if (v1 == 1 && v2 == 4)
                {
                    element_node2[3 + (tetra1 * 8 + 0) * element_order] = node;
                    element_node2[1 + (tetra1 * 8 + 3) * element_order] = node;
                    element_node2[2 + (tetra1 * 8 + 4) * element_order] = node;
                    element_node2[1 + (tetra1 * 8 + 6) * element_order] = node;
                }
                else if (v1 == 2 && v2 == 3)
                {
                    element_node2[2 + (tetra1 * 8 + 1) * element_order] = node;
                    element_node2[2 + (tetra1 * 8 + 2) * element_order] = node;
                    element_node2[2 + (tetra1 * 8 + 5) * element_order] = node;
                    element_node2[1 + (tetra1 * 8 + 7) * element_order] = node;
                }
                else if (v1 == 2 && v2 == 4)
                {
                    element_node2[3 + (tetra1 * 8 + 1) * element_order] = node;
                    element_node2[3 + (tetra1 * 8 + 2) * element_order] = node;
                    element_node2[2 + (tetra1 * 8 + 3) * element_order] = node;
                    element_node2[3 + (tetra1 * 8 + 4) * element_order] = node;
                    element_node2[3 + (tetra1 * 8 + 5) * element_order] = node;
                    element_node2[2 + (tetra1 * 8 + 6) * element_order] = node;
                    element_node2[2 + (tetra1 * 8 + 7) * element_order] = node;
                }
                else if (v1 == 3 && v2 == 4)
                {
                    element_node2[3 + (tetra1 * 8 + 3) * element_order] = node;
                    element_node2[3 + (tetra1 * 8 + 6) * element_order] = node;
                    element_node2[3 + (tetra1 * 8 + 7) * element_order] = node;
                }
            }
        }

        public static void tet_mesh_order4_refine_size(int node_num1, int element_num1,
                int[] element_node1, ref int node_num2, ref int element_num2, ref int[] edge_data)

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
            int element_order = 4;
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
            for (tetra = 0; tetra < element_num1; tetra++)
            {
                i = element_node1[0 + tetra * element_order];
                j = element_node1[1 + tetra * element_order];
                k = element_node1[2 + tetra * element_order];
                l = element_node1[3 + tetra * element_order];

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
            typeMethods.i4col_sort_a(5, 6 * element_num1, ref edge_data);
            //
            //  Step 3. All the tetrahedrons which share an edge show up as consecutive
            //  columns with identical first two entries.  Figure out how many new
            //  nodes there are, and allocate space for their coordinates.
            //
            node_num2 = node_num1;

            n1_old = -1;
            n2_old = -1;

            for (edge = 0; edge < 6 * element_num1; edge++)
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

            element_num2 = 8 * element_num1;

        }
    }
}