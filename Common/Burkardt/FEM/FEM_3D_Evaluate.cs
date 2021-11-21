using Burkardt.TetrahedronNS;

namespace Burkardt.FEM;

public static class FEM_3D_Evaluate
{
    public static double[] fem3d_evaluate(int fem_node_num, double[] fem_node_xyz,
            int fem_element_order, int fem_element_num, int[] fem_element_node,
            int[] fem_element_neighbor, int fem_value_dim, double[] fem_value,
            int sample_node_num, double[] sample_node_xyz)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FEM3D_EVALUATE samples an FEM function on a T4 tet mesh.
        //
        //  Discussion:
        //
        //    Note that the sample values returned are true values of the underlying
        //    finite element function.  They are NOT produced by constructing some
        //    other function that interpolates the data at the finite element nodes
        //    (something which MATLAB's griddata function can easily do.)  Instead, 
        //    each sampling node is located within one of the associated finite
        //    element tetrahedrons, and the finite element function is developed and 
        //    evaluated there.  
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int FEM_NODE_NUM, the number of nodes.
        //
        //    Input, double FEM_NODE_XYZ[3*FEM_NODE_NUM], the coordinates 
        //    of the nodes.
        //
        //    Input, int FEM_ELEMENT_ORDER, the order of the elements, 
        //    which should be 4.
        //
        //    Input, int FEM_ELEMENT_NUM, the number of elements.
        //
        //    Input, int FEM_ELEMENT_NODE[FEM_ELEMENT_ORDER*FEM_ELEMENT_NUM], the
        //    nodes that make up each element.
        //
        //    Input, int FEM_ELEMENT_NEIGHBOR[4*FEM_ELEMENT_NUM], the 
        //    index of the neighboring element on each side, or -1 if no neighbor there.
        //
        //    Input, int FEM_VALUE_DIM, the "dimension" of the values.
        //
        //    Input, double FEM_VALUE[FEM_VALUE_DIM*FEM_NODE_NUM], the 
        //    finite element coefficient values at each node.
        //
        //    Input, int SAMPLE_NODE_NUM, the number of sample nodes.
        //
        //    Input, double SAMPLE_NODE_XYZ[3*SAMPLE_NODE_NUM], the sample nodes.
        //
        //    Output, double) FEM3D_EVALUATE[FEM_VALUE_DIM*SAMPLE_NODE_NUM],
        //    the sampled values.
        //
    {
        int j;
        double[] p_xyz = new double[3];

        double[] b = new double[fem_element_order];
        double[] sample_value = new double[fem_value_dim * sample_node_num];
        int[] t_node = new int[fem_element_order];
        double[] t_xyz = new double[3 * fem_element_order];
        //
        //  For each sample point: find the element T that contains it,
        //  and evaluate the finite element function there.
        //
        for (j = 0; j < sample_node_num; j++)
        {
            p_xyz[0] = sample_node_xyz[0 + j * 3];
            p_xyz[1] = sample_node_xyz[1 + j * 3];
            p_xyz[2] = sample_node_xyz[2 + j * 3];
            //
            //  Find the element T that contains the point.
            //
            int step_num = 0;
            int t = TetMesh.tet_mesh_search_naive(fem_node_num, fem_node_xyz,
                fem_element_order, fem_element_num, fem_element_node,
                p_xyz, ref step_num);
            //
            //  Evaluate the finite element basis functions at the point in T.
            //
            int i;
            for (i = 0; i < fem_element_order; i++)
            {
                t_node[i] = fem_element_node[i + t * fem_element_order];
            }

            for (i = 0; i < fem_element_order; i++)
            {
                t_xyz[0 + i * 3] = fem_node_xyz[0 + t_node[i] * 3];
                t_xyz[1 + i * 3] = fem_node_xyz[1 + t_node[i] * 3];
                t_xyz[2 + i * 3] = fem_node_xyz[2 + t_node[i] * 3];
            }

            Basis_mn.basis_mn_tet4(t_xyz, 1, p_xyz, ref b);
            //
            //  Multiply by the finite element values to get the sample values.
            //
            for (i = 0; i < fem_value_dim; i++)
            {
                sample_value[i + j * fem_value_dim] = 0.0;
                int k;
                for (k = 0; k < fem_element_order; k++)
                {
                    sample_value[i + j * fem_value_dim] += fem_value[i + t_node[k] * fem_value_dim] * b[k];
                }
            }
        }

        return sample_value;
    }
}