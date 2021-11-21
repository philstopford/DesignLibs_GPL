using System;

namespace Burkardt.IntegralNS;

public static class Sum2
{
    public static double sum2_nd(Func<int, double[], double> func, double[] xtab,
            double[] weight, int[] order, int dim_num, ref int eval_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SUM2_ND estimates a multidimensional integral using a product rule.
        //
        //  Discussion:
        //
        //    The routine uses a product rule supplied by the user.
        //
        //    The region may be a product of any combination of finite,
        //    semi-infinite, or infinite intervals.
        //
        //    For each factor in the region, it is assumed that an integration
        //    rule is given, and hence, the region is defined implicitly by
        //    the integration rule chosen.
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
        //    Input, double XTAB[DIM_NUM*ORDER_MAX].  XTAB(I,J) is the 
        //    I-th abscissa of the J-th rule.
        //
        //    Input, double WEIGHT[DIM_NUM*ORDER_MAX].  WEIGHT(I,J) is the 
        //    I-th weight for the J-th rule.
        //
        //    Input, int ORDER[DIM_NUM].  ORDER(I) is the number of
        //    abscissas to be used in the J-th rule.  ORDER(I) must be
        //    greater than 0 and less than or equal to ORDER_MAX.
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Output, int EVAL_NUM, the number of function evaluations.
        //
        //    Output, double SUM2_ND, the approximate value of the integral.
        //
    {
        int dim;
        int i;
        //
        //  Default values.
        //
        double result = 0.0;
        eval_num = 0;

        switch (dim_num)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("SUM2_ND - Fatal error!");
                Console.WriteLine("  DIM_NUM < 1");
                Console.WriteLine("  DIM_NUM = " + dim_num + "");
                return 1;
        }

        for (i = 0; i < dim_num; i++)
        {
            switch (order[i])
            {
                case < 1:
                    Console.WriteLine("");
                    Console.WriteLine("SUM2_ND - Fatal error!");
                    Console.WriteLine("  ORDER(I) < 1.");
                    Console.WriteLine("  For I = " + i + "");
                    Console.WriteLine("  ORDER(I) = " + order[i] + "");
                    return 1;
            }
        }

        int[] iwork = new int[dim_num];
        double[] work = new double[dim_num];

        for (dim = 0; dim < dim_num; dim++)
        {
            iwork[dim] = 1;
        }

        for (;;)
        {
            int k = 1;

            double w1 = 1.0;
            for (i = 0; i < dim_num; i++)
            {
                int m1 = iwork[i];
                work[i] = xtab[i + (m1 - 1) * dim_num];
                w1 *= weight[i + (m1 - 1) * dim_num];
            }

            result += w1 * func(dim_num, work);
            eval_num += 1;

            while (iwork[k - 1] == order[k - 1])
            {
                iwork[k - 1] = 1;
                k += 1;

                if (dim_num < k)
                {
                    return result;
                }
            }

            iwork[k - 1] += 1;
        }
    }
}