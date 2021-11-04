using System;
using Burkardt.Types;

namespace Burkardt.IntegralNS
{
    public static class Romberg
    {
        public static double romberg_nd(Func<int, double[], double> func, double[] a,
                double[] b, int dim_num, int[] sub_num, int it_max, double tol, ref int ind,
                ref int eval_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ROMBERG_ND estimates a multidimensional integral using Romberg integration.
            //
            //  Discussion:
            //
            //    The routine uses a Romberg method based on the midpoint rule.
            //
            //    In the reference, this routine is called "NDIMRI".
            //
            //    Thanks to Barak Bringoltz for pointing out problems in a previous
            //    FORTRAN90 implementation of this routine.
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
            //    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
            //    C++ version by John Burkardt
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
            //    Input, double A[DIM_NUM], B[DIM_NUM], the integration limits.
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, int SUB_NUM[DIM_NUM], the number of subintervals into
            //    which the I-th integration interval (A(I), B(I)) is
            //    initially subdivided.  SUB_NUM(I) must be greater than 0.
            //
            //    Input, int IT_MAX, the maximum number of iterations to
            //    be performed.  The number of function evaluations on
            //    iteration J is at least J**DIM_NUM, which grows very rapidly.
            //    IT_MAX should be small!
            //
            //    Input, double TOL, an error tolerance for the approximation
            //    of the integral.
            //
            //    Output, int *IND, error return flag.
            //    IND = -1 if the error tolerance could not be achieved.
            //    IND = 1 if the error tolerance was achieved.
            //
            //    Output, int eval_NUM, the number of function evaluations.
            //
            //    Output, double ROMBERG_ND, the approximate value of the integral.
            //
            //  Local Parameters:
            //
            //    Local, int IWORK[DIM_NUM], a pointer used to generate all the
            //    points X in the product region.
            //
            //    Local, int IWORK2[IT_MAX], a counter of the number of points
            //    used at each step of the Romberg iteration.
            //
            //    Local, int SUB_NUM2[DIM_NUM], the number of subintervals used
            //    in each direction, a refinement of the user's input SUB_NUM.
            //
            //    Local, double TABLE[IT_MAX], the difference table.
            //
            //    Local, double X[DIM_NUM], an evaluation point.
            //
        {
            int dim;
            double factor;
            int i;
            int it;
            int[] iwork;
            int[] iwork2;
            int kdim;
            int ll;
            double result;
            double result_old = 0;
            double rnderr;
            int[] sub_num2;
            double sum1;
            double weight;
            double[] table;
            double[] x;

            eval_num = 0;

            if (dim_num < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("ROMBERG_ND - Fatal error!");
                Console.WriteLine("  DIM_NUM is less than 1.  DIM_NUM = " + dim_num + "");
                return (1);
            }

            if (it_max < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("ROMBERG_ND - Fatal error!");
                Console.WriteLine("  IT_MAX is less than 1.  IT_MAX = " + it_max + "");
                return (1);
            }

            for (i = 0; i < dim_num; i++)
            {
                if (sub_num[i] <= 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("ROMBERG_ND - Fatal error!");
                    Console.WriteLine("  SUB_NUM(I) is less than 1.");
                    Console.WriteLine("  for I = " + i + "");
                    Console.WriteLine("  SUB_NUM(I) = " + sub_num[i] + "");
                    return (1);
                }
            }

            iwork = new int[dim_num];
            iwork2 = new int[it_max];
            sub_num2 = new int[dim_num];
            table = new double[it_max];
            x = new double[dim_num];

            ind = 0;
            rnderr = typeMethods.r8_epsilon();
            iwork2[0] = 1;

            for (dim = 0; dim < dim_num; dim++)
            {
                sub_num2[dim] = sub_num[dim];
            }

            if (1 < it_max)
            {
                iwork2[1] = 2;
            }

            it = 1;

            for (;;)
            {
                sum1 = 0.0;

                weight = 1.0;
                for (dim = 0; dim < dim_num; dim++)
                {
                    weight = weight * (b[dim] - a[dim]) / (double) sub_num2[dim];
                }

                //
                //  Generate every point X in the product region, and evaluate F(X).
                //
                for (dim = 0; dim < dim_num; dim++)
                {
                    iwork[dim] = 1;
                }

                for (;;)
                {
                    for (dim = 0; dim < dim_num; dim++)
                    {
                        x[dim] =
                            ((double) (2 * sub_num2[dim] - 2 * iwork[dim] + 1) * a[dim]
                             + (double) (+2 * iwork[dim] - 1) * b[dim])
                            / (double) (2 * sub_num2[dim]);
                    }

                    sum1 = sum1 + func(dim_num, x);
                    eval_num = eval_num + 1;

                    kdim = dim_num;

                    while (0 < kdim)
                    {
                        if (iwork[kdim - 1] < sub_num2[kdim - 1])
                        {
                            iwork[kdim - 1] = iwork[kdim - 1] + 1;
                            break;
                        }

                        iwork[kdim - 1] = 1;
                        kdim = kdim - 1;
                    }

                    if (kdim == 0)
                    {
                        break;
                    }
                }

                //
                //  Done with summing.
                //
                table[it - 1] = weight * sum1;

                if (it <= 1)
                {
                    result = table[0];
                    result_old = result;

                    if (it_max <= it)
                    {
                        ind = 1;
                        break;
                    }

                    it = it + 1;
                    for (dim = 0; dim < dim_num; dim++)
                    {
                        sub_num2[dim] = iwork2[it - 1] * sub_num2[dim];
                    }

                    continue;
                }

                //
                //  Compute the difference table for Richardson extrapolation.
                // 
                for (ll = 2; ll <= it; ll++)
                {
                    i = it + 1 - ll;
                    factor = (double) (iwork2[i - 1] * iwork2[i - 1])
                             / (double) (iwork2[it - 1] * iwork2[it - 1] - iwork2[i - 1] * iwork2[i - 1]);
                    table[i] = table[i] + (table[i] - table[i - 1]) * factor;
                }

                result = table[0];
                //
                //  Terminate successfully if the estimated error is acceptable.
                //
                if (Math.Abs(result - result_old) <= Math.Abs(result * (tol + rnderr)))
                {
                    ind = 1;
                    break;
                }

                //
                //  Terminate unsuccessfully if the iteration limit has been reached.
                //
                if (it_max <= it)
                {
                    ind = -1;
                    break;
                }

                //
                //  Prepare for another step.
                //
                result_old = result;

                it = it + 1;

                iwork2[it - 1] = (int) (1.5 * (double) (iwork2[it - 2]));

                for (dim = 0; dim < dim_num; dim++)
                {
                    sub_num2[dim] = (int) (1.5 * (double) (sub_num2[dim]));
                }
            }

            return result;
        }
    }
}