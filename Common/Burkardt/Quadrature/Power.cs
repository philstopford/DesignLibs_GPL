using System;
using Burkardt.Types;

namespace Burkardt
{
    public static class PowerQuadrature
    {
        public static void power_rule_set(int point_num_1d, double[] x_1d, double[] w_1d,
                double[] r_1d, int dim_num, int point_num, ref double[] x,
                ref double[] w, ref double[] r)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POWER_RULE_SET sets up a power rule.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 May 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int POINT_NUM_1D, the order of the 1D rule.
            //
            //    Input, double X_1D[POINT_NUM_1D], the points of the 1D rule.
            //
            //    Input, double W_1D[POINT_NUM_1D], the weights of the
            //    1D rule.
            //
            //    Input, double R_1D[2], the extreme points that define
            //    the range of the 1D region.
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, int POINT_NUM, the number of points in the rule.
            //
            //    Output, double X[DIM_NUM*POINT_NUM], the points of the rule.
            //
            //    Output, double W[POINT_NUM], the weights of the rule.
            //
            //    Output, double R[DIM_NUM*2], the extreme points 
            //    that define the range of the product rule region.
            //
        {
            int dim;
            int[] indx;
            int k;

            indx = new int[dim_num];
            k = 0;

            for (;;)
            {
                BTuple.tuple_next(0, point_num_1d - 1, dim_num, ref k, ref indx);

                if (k == 0)
                {
                    break;
                }

                w[k - 1] = 1.0;
                for (dim = 0; dim < dim_num; dim++)
                {
                    w[k - 1] = w[k - 1] * w_1d[indx[dim]];
                }

                for (dim = 0; dim < dim_num; dim++)
                {
                    x[dim + (k - 1) * dim_num] = x_1d[indx[dim]];
                }
            }

            for (dim = 0; dim < dim_num; dim++)
            {
                r[dim + 0 * dim_num] = r_1d[0];
                r[dim + 1 * dim_num] = r_1d[1];
            }

        }

        public static int power_rule_size(int point_num_1d, int dim_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POWER_RULE_SIZE returns the size of a power rule.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 May 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int POINT_NUM_1D, the number of points in the 1D rule.
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Output, int POWER_RULE_SIZE, the number of points in the rule.
            //
        {
            int point_num;

            point_num = (int)Math.Pow(point_num_1d, dim_num);

            return point_num;
        }
    }
}