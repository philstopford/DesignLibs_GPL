namespace Burkardt.TetrahedronNS;

public static class TetMesh_Q2L
{
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
        //    The same nodes are used, but there are s8 times as many
        //    tetrahedrons.
        //
        //    The node ordering for the quadratic tetrahedron is somewhat
        //    arbitrary.  In the current scheme, the vertices are listed
        //    first, followed by the 6 midside nodes.  Each midside node
        //    may be identified by the two vertices that bracket it.  Thus,
        //    the node ordering may be suggested by:
        //
        //     1  2  3  4 (1+2) (1+3) (1+4) (2+3) (2+4) (3+4)
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
        int tetra1;

        int tetra2 = 0;

        for (tetra1 = 0; tetra1 < tetra_num1; tetra1++)
        {
            int n1 = tetra_node1[0 + tetra1 * 10];
            int n2 = tetra_node1[1 + tetra1 * 10];
            int n3 = tetra_node1[2 + tetra1 * 10];
            int n4 = tetra_node1[3 + tetra1 * 10];
            int n5 = tetra_node1[4 + tetra1 * 10];
            int n6 = tetra_node1[5 + tetra1 * 10];
            int n7 = tetra_node1[6 + tetra1 * 10];
            int n8 = tetra_node1[7 + tetra1 * 10];
            int n9 = tetra_node1[8 + tetra1 * 10];
            int nx = tetra_node1[9 + tetra1 * 10];

            tetra_node2[0 + tetra2 * 4] = n1;
            tetra_node2[1 + tetra2 * 4] = n5;
            tetra_node2[2 + tetra2 * 4] = n6;
            tetra_node2[3 + tetra2 * 4] = n7;
            tetra2 += 1;

            tetra_node2[0 + tetra2 * 4] = n2;
            tetra_node2[1 + tetra2 * 4] = n5;
            tetra_node2[2 + tetra2 * 4] = n8;
            tetra_node2[3 + tetra2 * 4] = n9;
            tetra2 += 1;

            tetra_node2[0 + tetra2 * 4] = n3;
            tetra_node2[1 + tetra2 * 4] = n6;
            tetra_node2[2 + tetra2 * 4] = n8;
            tetra_node2[3 + tetra2 * 4] = n9;
            tetra2 += 1;

            tetra_node2[0 + tetra2 * 4] = n4;
            tetra_node2[1 + tetra2 * 4] = n7;
            tetra_node2[2 + tetra2 * 4] = n9;
            tetra_node2[3 + tetra2 * 4] = nx;
            tetra2 += 1;

            tetra_node2[0 + tetra2 * 4] = n5;
            tetra_node2[1 + tetra2 * 4] = n6;
            tetra_node2[2 + tetra2 * 4] = n7;
            tetra_node2[3 + tetra2 * 4] = n9;
            tetra2 += 1;

            tetra_node2[0 + tetra2 * 4] = n5;
            tetra_node2[1 + tetra2 * 4] = n6;
            tetra_node2[2 + tetra2 * 4] = n8;
            tetra_node2[3 + tetra2 * 4] = n9;
            tetra2 += 1;

            tetra_node2[0 + tetra2 * 4] = n6;
            tetra_node2[1 + tetra2 * 4] = n7;
            tetra_node2[2 + tetra2 * 4] = n9;
            tetra_node2[3 + tetra2 * 4] = nx;
            tetra2 += 1;

            tetra_node2[0 + tetra2 * 4] = n6;
            tetra_node2[1 + tetra2 * 4] = n8;
            tetra_node2[2 + tetra2 * 4] = n9;
            tetra_node2[3 + tetra2 * 4] = nx;
            tetra2 += 1;
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
}