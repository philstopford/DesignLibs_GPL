using System;
using Burkardt.TetrahedronNS;

namespace Burkardt.FEM
{
    public static class FEM_3D_Projection
    {
        public static double[] projection(int fem_node_num, double[] fem_node_xyz,
        int fem_element_order, int fem_element_num, int[] fem_element_node,
        int[] fem_element_neighbor, int fem_value_dim, double[] fem_value,
        int sample_node_num, double[] sample_node_xyz )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PROJECTION evaluates an FEM function on a T3 or T6 triangulation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 August 2009
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
        //    Input, int FEM_ELEMENT_ORDER, the order of the elements.
        //
        //    Input, int FEM_ELEMENT_NUM, the number of elements.
        //
        //    Input, int FEM_ELEMENT_NODE(FEM_ELEMENT_ORDER,FEM_ELEMENT_NUM), the
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
        //    Output, double PROJECTION[FEM_VALUE_DIM*SAMPLE_NODE_NUM],
        //    the sampled values.
        //
        {
            double[] b;
            double dot;
            int face = 0;
            int i;
            int j;
            int k;
            double[] p_xyz = new double[3];
            double[] sample_value;
            int step_num = 0;
            int t;
            int t_node;
            double[] t_xyz;

            b = new double[fem_element_order];
            sample_value = new double[fem_value_dim * sample_node_num];
            t_xyz = new double[3 * fem_element_order];
            //
            //  For each sample point: find the element T that contains it,
            //  and evaluate the finite element function there.
            //
            for (j = 0; j < sample_node_num; j++)
            {
                p_xyz[0] = sample_node_xyz[0 + 3 * j];
                p_xyz[1] = sample_node_xyz[1 + 3 * j];
                p_xyz[2] = sample_node_xyz[2 + 3 * j];
                //
                //  Find the triangle T that contains the point.
                //
                t = TetMesh.tet_mesh_search_delaunay(fem_node_num, fem_node_xyz,
                    fem_element_order, fem_element_num, fem_element_node,
                    fem_element_neighbor, p_xyz, ref face, ref step_num);

                if (t == -1)
                {
                    Console.WriteLine("");
                    Console.WriteLine("PROJECTION - Fatal error!");
                    Console.WriteLine("  Search failed.");
                    return null;
                }

                //
                //  Evaluate the finite element basis functions at the point in T.
                //
                for (i = 0; i < fem_element_order; i++)
                {
                    t_node = fem_element_node[i + t * fem_element_order];
                    t_xyz[0 + i * 3] = fem_node_xyz[0 + t_node * 3];
                    t_xyz[1 + i * 3] = fem_node_xyz[1 + t_node * 3];
                    t_xyz[2 + i * 3] = fem_node_xyz[2 + t_node * 3];
                }

                Basis_mn.basis_mn_tet4(t_xyz, 1, p_xyz, ref b);
                //
                //  Multiply by the finite element values to get the sample values.
                //
                for (i = 0; i < fem_value_dim; i++)
                {
                    dot = 0.0;
                    for (k = 0; k < fem_element_order; k++)
                    {
                        t_node = fem_element_node[k + t * fem_element_order];
                        dot = dot + fem_value[i + t_node * fem_value_dim] * b[k];
                    }

                    sample_value[i + j * fem_value_dim] = dot;
                }
            }

            return sample_value;
        }
    }
}