using System;

namespace Burkardt.IntegralNS;

public static class P5
{
    public static double p5_nd(Func<int, double[], double> func, int dim_num,
            double[] a, double[] b, ref int eval_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P5_ND estimates a multidimensional integral with a formula of exactness 5.
        //
        //  Discussion:
        //
        //    The routine uses a method which is exact for polynomials of total 
        //    degree 5 or less.
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
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, double A[DIM_NUM], B[DIM_NUM], the integration limits.
        //
        //    Output, int eval_num, the number of function evaluations.
        //
        //    Output, double P5_ND, the approximate value of the integral.
        //
    {
        int dim;
        int i;

        eval_num = 0;

        switch (dim_num)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("P5_ND - Fatal error!");
                Console.WriteLine("  DIM_NUM < 1, DIM_NUM = " + dim_num + "");
                return 1;
        }

        const double a2 = 25.0 / 324.0;
        double a3 = Math.Sqrt(0.6);
        double en = dim_num;
        double a0 = (25.0 * en * en - 115.0 * en + 162.0) / 162.0;
        double a1 = (70.0 - 25.0 * en) / 162.0;

        double volume = 1.0;
        for (dim = 0; dim < dim_num; dim++)
        {
            volume *= b[dim] - a[dim];
        }

        double[] work = new double[dim_num];
        for (dim = 0; dim < dim_num; dim++)
        {
            work[dim] = 0.5 * (a[dim] + b[dim]);
        }

        double result = 0.0;
        switch (volume)
        {
            case 0.0:
                Console.WriteLine("");
                Console.WriteLine("P5_ND - Warning!");
                Console.WriteLine("  Volume = 0, integral = 0.");
                return result;
        }

        double sum1 = a0 * func(dim_num, work);
        eval_num += 1;

        double sum2 = 0.0;
        double sum3 = 0.0;

        for (i = 0; i < dim_num; i++)
        {
            work[i] = 0.5 * (a[i] + b[i] + a3 * (b[i] - a[i]));
            sum2 += func(dim_num, work);
            eval_num += 1;

            work[i] = 0.5 * (a[i] + b[i] - a3 * (b[i] - a[i]));
            sum2 += func(dim_num, work);
            eval_num += 1;

            work[i] = 0.5 * (a[i] + b[i]);
        }

        switch (dim_num)
        {
            case > 1:
            {
                double a4 = a3;

                for (;;)
                {
                    for (i = 0; i < dim_num - 1; i++)
                    {
                        work[i] = 0.5 * (a[i] + b[i] + a4 * (b[i] - a[i]));
                        double a5 = a3;

                        for (;;)
                        {
                            int j;
                            for (j = i + 1; j < dim_num; j++)
                            {
                                work[j] = 0.5 * (a[j] + b[j] + a5 * (b[j] - a[j]));
                                sum3 += func(dim_num, work);
                                eval_num += 1;
                                work[j] = 0.5 * (a[j] + b[j]);
                            }

                            a5 = -a5;

                            if (0.0 <= a5)
                            {
                                break;
                            }
                        }

                        work[i] = 0.5 * (a[i] + b[i]);
                    }

                    a4 = -a4;

                    if (0.0 <= a4)
                    {
                        break;
                    }
                }

                break;
            }
        }

        result = volume * (sum1 + a1 * sum2 + a2 * sum3);

        return result;
    }
}