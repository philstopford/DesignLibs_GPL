using System;
using Burkardt.Types;

namespace Burkardt.FEM
{
    public static class FEM_1D_Sample
    {
        public static double[] fem1d_evaluate(int node_num, double[] node_x, int element_order,
                int element_num, int value_dim, double[] value, int sample_node_num,
                double[] sample_node_x)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FEM1D_EVALUATE evaluates a 1D FEM function at sample points.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 April 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NODE_NUM, the number of FEM nodes.
            //
            //    Input, double NODE_X[NODE_NUM], the nodes.  
            //
            //    Input, int ELEMENT_ORDER, the element order.
            //
            //    Input, int ELEMENT_NUM, the number of elements.
            //
            //    Input, int VALUE_DIM, the value dimension.
            //
            //    Input, double VALUE[VALUE_DIM*NODE_NUM], the FEM values.
            //
            //    Input, int SAMPLE_NODE_NUM, the number of sample points.
            //
            //    INput, double SAMPLE_NODE_X[SAMPLE_NODE_NUM], the sample nodes.
            //
            //    Output, double FEM1D_EVALUATE[VALUE_DIM*SAMPLE_NODE_NUM],
            //    the interpolated FEM values at sample nodes.
            //
        {
            //
            //  For each sample point, find NODE_LEFT and NODE_RIGHT that bracket it.
            //
            int[] sample_left = new int[sample_node_num];

            typeMethods.r8vec_bracket4(node_num, node_x, sample_node_num, sample_node_x,
                ref sample_left);

            double[] sample_value = new double[value_dim * sample_node_num];

            if (element_order == 1)
            {
                for (int sample = 0; sample < sample_node_num; sample++)
                {
                    for (int i = 0; i < value_dim; i++)
                    {
                        sample_value[i + sample * value_dim] = value[i + sample_left[sample] * value_dim];
                    }
                }
            }
            else if (element_order == 2)
            {
                for (int sample = 0; sample < sample_node_num; sample++)
                {
                    int l = sample_left[sample];
                    int r = sample_left[sample] + 1;
                    for (int i = 0; i < value_dim; i++)
                    {
                        sample_value[i + sample * value_dim] =
                            (value[i + l * value_dim] * (node_x[r] - sample_node_x[sample])
                             + value[i + r * value_dim] * (sample_node_x[sample] - node_x[l]))
                            / (node_x[r] - node_x[l]);
                    }
                }
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("FEM_EVALUATE - Fatal error!");
                Console.WriteLine("  Cannot handle elements of this order.");
                return new double[1];
            }

            return sample_value;
        }
    }
}