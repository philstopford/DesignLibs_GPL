using System;
using Burkardt.Types;

namespace Burkardt.IntegralNS
{
    public static class Box
    {
        public static double box_nd(Func<int, double[], double> func, int dim_num,
                int order, double[] xtab, double[] weight, ref int eval_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BOX_ND estimates a multidimensional integral using a product rule.
            //
            //  Discussion:
            //
            //    The routine creates a DIM_NUM-dimensional product rule from a 1D rule
            //    supplied by the user.  The routine is fairly inflexible.  If
            //    you supply a rule for integration from -1 to 1, then your product
            //    box must be a product of DIM_NUM copies of the interval [-1,1].
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    25 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Philip Davis, Philip Rabinowitz,
            //    Methods of Numerical Integration,
            //    Second Edition,
            //    Dover, 2007,
            //    ISBN: 0486453391,
            //    LC: QA299.3.D28.
            //
            //  Parameters:
            //
            //    Input, double FUNC ( int dim_num, double x[] ), evaluates
            //    the function to be integrated.
            //      
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, int ORDER, the number of points used in the 1D rule.
            //
            //    Input, double XTAB[ORDER], the abscissas of the 1D rule.
            //
            //    Input, double WEIGHT[ORDER], the weights of the 1D rule.
            //
            //    Output, int *EVAL_NUM, the number of function evaluations.
            //
            //    Output, double BOX_ND, the approximate value of the integral.
            //
        {
            int dim;
            int[] indx;
            int k;
            double result;
            double w;
            double[] x;

            eval_num = 0;

            if (dim_num < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("BOX_ND - Fatal error!");
                Console.WriteLine("  DIM_NUM < 1.");
                Console.WriteLine("  DIM_NUM = " + dim_num + "");
                return (1);
            }

            if (order < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("BOX_ND - Fatal error!");
                Console.WriteLine("  ORDER < 1.");
                Console.WriteLine("  ORDER = " + order + "");
                return (1);
            }

            k = 0;
            result = 0.0;

            indx = new int[dim_num];
            x = new double[dim_num];

            for (;;)
            {
                BTuple.tuple_next(1, order, dim_num, ref k, ref indx);

                if (k == 0)
                {
                    break;
                }

                w = 1.0;
                for (dim = 0; dim < dim_num; dim++)
                {
                    w = w * weight[indx[dim] - 1];
                }

                for (dim = 0; dim < dim_num; dim++)
                {
                    x[dim] = xtab[indx[dim] - 1];
                }

                result = result + w * func(dim_num, x);
                eval_num = eval_num + 1;
            }


            return result;
        }
    }
}