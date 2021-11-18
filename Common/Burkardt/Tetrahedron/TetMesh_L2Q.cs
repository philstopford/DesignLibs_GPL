using System;
using Burkardt.Types;

namespace Burkardt.TetrahedronNS;

public static class TetMesh_L2Q
{
    public static void tet_mesh_order4_to_order10_compute(int element_num, int[] element_node1,
            int node_num1, double[] node_xyz1, int[] edge_data, ref int[] element_node2,
            int node_num2, ref double[] node_xyz2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TET_MESH_ORDER4_TO_ORDER10_COMPUTE: quadratic tet mesh from a linear one.
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
        //    Thanks to kaushikkn for pointing out an indexing error in this
        //    function, 24 May 2017.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 May 2017
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
        const int dim_num = 3;
        int edge;
        int i;
        int j;
        const int element_order1 = 4;
        const int element_order2 = 10;
        int v = 0;
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

        for (j = 0; j < element_num; j++)
        {
            for (i = 0; i < element_order1; i++)
            {
                element_node2[i + j * element_order2] = element_node1[i + j * element_order1];
            }
        }

        int node = node_num1;

        int n1_old = -1;
        int n2_old = -1;

        for (edge = 0; edge < 6 * element_num; edge++)
        {
            //
            //  Read the data defining the edge.
            //
            int n1 = edge_data[0 + edge * 5];
            int n2 = edge_data[1 + edge * 5];
            //
            //  If this edge is new, we need to create a new node between
            //  the endpoints.
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
                        (node_xyz2[i + n1 * dim_num] + node_xyz2[i + n2 * dim_num]) / 2.0;
                }

                node += 1;
                n1_old = n1;
                n2_old = n2;
            }

            //
            //  Assign the node to the tetrahedron.
            //
            int v1 = edge_data[2 + edge * 5];
            int v2 = edge_data[3 + edge * 5];
            v = v1 switch
            {
                //
                //  Here is where the local ordering of the nodes is effected:
                //
                0 when v2 == 1 => 4,
                0 when v2 == 2 => 5,
                0 when v2 == 3 => 6,
                1 when v2 == 2 => 7,
                1 when v2 == 3 => 8,
                2 when v2 == 3 => 9,
                _ => v
            };

            int tetra = edge_data[4 + edge * 5];

            element_node2[v + tetra * element_order2] = node - 1;
        }

    }

    public static void tet_mesh_order4_to_order10_size(int element_num, int[] element_node1,
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
        int tetra;
        const int element_order1 = 4;
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
        for (tetra = 0; tetra < element_num; tetra++)
        {
            int i = element_node1[0 + tetra * element_order1];
            int j = element_node1[1 + tetra * element_order1];
            int k = element_node1[2 + tetra * element_order1];
            int l = element_node1[3 + tetra * element_order1];

            typeMethods.i4i4_sort_a(i, j, ref a, ref b);

            edge_data[0 + 6 * tetra * 5] = a;
            edge_data[1 + 6 * tetra * 5] = b;
            edge_data[2 + 6 * tetra * 5] = 0;
            edge_data[3 + 6 * tetra * 5] = 1;
            edge_data[4 + 6 * tetra * 5] = tetra;

            typeMethods.i4i4_sort_a(i, k, ref a, ref b);

            edge_data[0 + (6 * tetra + 1) * 5] = a;
            edge_data[1 + (6 * tetra + 1) * 5] = b;
            edge_data[2 + (6 * tetra + 1) * 5] = 0;
            edge_data[3 + (6 * tetra + 1) * 5] = 2;
            edge_data[4 + (6 * tetra + 1) * 5] = tetra;

            typeMethods.i4i4_sort_a(i, l, ref a, ref b);

            edge_data[0 + (6 * tetra + 2) * 5] = a;
            edge_data[1 + (6 * tetra + 2) * 5] = b;
            edge_data[2 + (6 * tetra + 2) * 5] = 0;
            edge_data[3 + (6 * tetra + 2) * 5] = 3;
            edge_data[4 + (6 * tetra + 2) * 5] = tetra;

            typeMethods.i4i4_sort_a(j, k, ref a, ref b);

            edge_data[0 + (6 * tetra + 3) * 5] = a;
            edge_data[1 + (6 * tetra + 3) * 5] = b;
            edge_data[2 + (6 * tetra + 3) * 5] = 1;
            edge_data[3 + (6 * tetra + 3) * 5] = 2;
            edge_data[4 + (6 * tetra + 3) * 5] = tetra;

            typeMethods.i4i4_sort_a(j, l, ref a, ref b);

            edge_data[0 + (6 * tetra + 4) * 5] = a;
            edge_data[1 + (6 * tetra + 4) * 5] = b;
            edge_data[2 + (6 * tetra + 4) * 5] = 1;
            edge_data[3 + (6 * tetra + 4) * 5] = 3;
            edge_data[4 + (6 * tetra + 4) * 5] = tetra;

            typeMethods.i4i4_sort_a(k, l, ref a, ref b);

            edge_data[0 + (6 * tetra + 5) * 5] = a;
            edge_data[1 + (6 * tetra + 5) * 5] = b;
            edge_data[2 + (6 * tetra + 5) * 5] = 2;
            edge_data[3 + (6 * tetra + 5) * 5] = 3;
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
        typeMethods.i4col_sort_a(5, 6 * element_num, ref edge_data);
        //
        //  Step 3. All the tetrahedrons which share an edge show up as consecutive
        //  columns with identical first two entries.  Figure out how many new
        //  nodes there are, and allocate space for their coordinates.
        //
        node_num2 = node_num1;

        int n1_old = -1;
        int n2_old = -1;

        for (edge = 0; edge < 6 * element_num; edge++)
        {
            int n1 = edge_data[0 + edge * 5];
            int n2 = edge_data[1 + edge * 5];
            if (n1 == n1_old && n2 == n2_old)
            {
                continue;
            }

            node_num2 += 1;
            n1_old = n1;
            n2_old = n2;
        }
    }
}