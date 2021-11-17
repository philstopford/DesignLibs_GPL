using System;
using Burkardt.TriangulationNS;

namespace Burkardt.FEM;

public static class FEM_2D_Projection
{
    public static double[] projection(int fem_node_num, double[] fem_node_xy,
            int fem_element_order, int fem_element_num, int[] fem_element_node,
            int[] fem_element_neighbor, int fem_value_dim, double[] fem_value,
            int sample_node_num, double[] sample_node_xy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PROJECTION evaluates an FEM function on a T3 or T6 triangulation.
        //
        //  Discussion:
        //
        //    Note that the sample values returned are true values of the underlying
        //    finite element function.  They are NOT produced by constructing some
        //    other function that interpolates the data at the finite element nodes
        //    (something which MATLAB's griddata function can easily do.)  Instead, 
        //    each sampling node is located within one of the associated finite
        //    element triangles, and the finite element function is developed and 
        //    evaluated there.  
        //
        //    MATLAB's scattered data interpolation is wonderful, but it cannot
        //    be guaranteed to reproduce the finite element function corresponding
        //    to nodal data.  This routine can (or at least tries to//).
        // 
        //    So if you are using finite elements, then using THIS routine
        //    (but not MATLAB's griddata function), what you see is what you have//
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int FEM_NODE_NUM, the number of nodes.
        //
        //    Input, double FEM_NODE_XY[2*FEM_NODE_NUM], the coordinates 
        //    of the nodes.
        //
        //    Input, int FEM_ELEMENT_ORDER, the order of the elements, 
        //    either 3 or 6.
        //
        //    Input, int FEM_ELEMENT_NUM, the number of triangles.
        //
        //    Input, int FEM_ELEMENT_NODE(FEM_ELEMENT_ORDER,FEM_ELEMENT_NUM), the
        //    nodes that make up each triangle.
        //
        //    Input, int FEM_ELEMENT_NEIGHBOR[3*FEM_ELEMENT_NUM], the 
        //    index of the neighboring triangle on each side, or -1 if no neighbor there.
        //
        //    Input, int FEM_VALUE_DIM, the "dimension" of the values.
        //
        //    Input, double FEM_VALUE[FEM_VALUE_DIM*FEM_NODE_NUM], the 
        //    finite element coefficient values at each node.
        //
        //    Input, int SAMPLE_NODE_NUM, the number of sample nodes.
        //
        //    Input, double SAMPLE_NODE_XY[2*SAMPLE_NODE_NUM], the sample nodes.
        //
        //    Output, double PROJECTION[FEM_VALUE_DIM*SAMPLE_NODE_NUM],
        //    the sampled values.
        //
    {
        double[] b;
        double[] dbdx;
        double[] dbdy;
        double dot;
        int edge = 0;
        int i;
        int j;
        int k;
        double[] p_xy = new double[2];
        double[] sample_value;
        int t = 0;
        int[] t_node;
        double[] t_xy;
        DelaunaySearchData data = new();
        double alpha = 0;
        double beta = 0;
        double gammaa = 0;
        int step = 0;

        b = new double[fem_element_order];
        dbdx = new double[fem_element_order];
        dbdy = new double[fem_element_order];
        sample_value = new double[fem_value_dim * sample_node_num];
        t_node = new int[fem_element_order];
        t_xy = new double[2 * fem_element_order];
        //
        //  For each sample point: find the triangle T that contains it,
        //  and evaluate the finite element function there.
        //
        for (j = 0; j < sample_node_num; j++)
        {
            p_xy[0] = sample_node_xy[0 + 2 * j];
            p_xy[1] = sample_node_xy[1 + 2 * j];
            //
            //  Find the triangle T that contains the point.
            //
            Search.triangulation_search_delaunay_a(ref data, fem_node_num, fem_node_xy,
                fem_element_order, fem_element_num, fem_element_node,
                fem_element_neighbor, p_xy, ref t, ref alpha, ref beta, ref gammaa, ref edge, ref step);

            switch (t)
            {
                case -1:
                    Console.WriteLine("");
                    Console.WriteLine("PROJECTION - Fatal error!");
                    Console.WriteLine("  Triangulation search failed.");
                    return null;
            }

            //
            //  Evaluate the finite element basis functions at the point in T.
            //
            for (i = 0; i < fem_element_order; i++)
            {
                t_node[i] = fem_element_node[(i + t * fem_element_order)  % fem_element_node.Length];
                t_xy[0 + i * 2] = fem_node_xy[0 + t_node[i] * 2];
                t_xy[1 + i * 2] = fem_node_xy[1 + t_node[i] * 2];
            }

            switch (fem_element_order)
            {
                case 3:
                    Basis_mn.basis_mn_t3(t_xy, 1, p_xy, ref b, ref dbdx, ref dbdy);
                    break;
                case 6:
                    Basis_mn.basis_mn_t6(t_xy, 1, p_xy, ref b, ref dbdx, ref dbdy);
                    break;
            }

            //
            //  Multiply by the finite element values to get the sample values.
            //
            for (i = 0; i < fem_value_dim; i++)
            {
                dot = 0.0;
                for (k = 0; k < fem_element_order; k++)
                {
                    dot += fem_value[i + t_node[k] * fem_value_dim] * b[k];
                }

                sample_value[i + j * fem_value_dim] = dot;
            }
        }

        return sample_value;
    }
}