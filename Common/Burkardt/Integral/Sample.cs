using System;

namespace Burkardt.IntegralNS
{
    public static class Sample
    {
        public static void sample_nd(Func<int, double[], double> func, int k1, int k2,
                int dim_num, ref double[] est1, ref double[] err1, ref double[] dev1, ref double[] est2,
                ref double[] err2, ref double[] dev2, ref int eval_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SAMPLE_ND estimates a multidimensional integral using sampling.
            //
            //  Discussion:
            //
            //    This routine computes two sequences of integral estimates, EST1 
            //    and EST2, for indices K going from K1 to K2.  These estimates are 
            //    produced by the generation of 'random' abscissas in the region.  
            //    The process can become very expensive if high accuracy is needed.
            //
            //    The total number of function evaluations is
            //    4*(K1^DIM_NUM+(K1+1)^DIM_NUM+...+(K2-1)^DIM_NUM+K2^DIM_NUM), and K2
            //    should be chosen so as to make this quantity reasonable.
            //    In most situations, EST2(K) are much better estimates than EST1(K).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    01 March 2007
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
            //    Input, int K1, the beginning index for the iteration.
            //    1 <= K1 <= K2.
            //
            //    Input, int K2, the final index for the iteration.  K1 <= K2.
            //    Increasing K2 increases the accuracy of the calculation,
            //    but vastly increases the work and running time of the code.
            //
            //    Input, int DIM_NUM, the spatial dimension.  1 <= DIM_NUM <= 10.
            //
            //    Output, double EST1[K2].  Entries K1 through K2 contain
            //    successively better estimates of the integral.
            //
            //    Output, double ERR1[K2].  Entries K1 through K2 contain
            //    the corresponding estimates of the integration errors.
            //
            //    Output, double DEV1[K2].  Entries K1 through K2 contain
            //    estimates of the reliability of the the integration.
            //    If consecutive values DEV1(K) and DEV1(K+1) do not differ
            //    by more than 10 percent, then ERR1(K) can be taken as
            //    a reliable upper bound on the difference between EST1(K)
            //    and the true value of the integral.
            //
            //    Output, double EST2[K2].  Entries K2 through K2 contain
            //    successively better estimates of the integral.
            //
            //    Output, double ERR2[K2].  Entries K2 through K2 contain
            //    the corresponding estimates of the integration errors.
            //
            //    Output, double DEV2[K2].  Entries K2 through K2 contain
            //    estimates of the reliability of the the integration.
            //    If consecutive values DEV2(K) and DEV2(K+2) do not differ
            //    by more than 10 percent, then ERR2(K) can be taken as
            //    a reliable upper bound on the difference between EST2(K)
            //    and the true value of the integral.
            //
            //    Output, int eval_NUM, the number of function evaluations.
            //
        {
            int DIM_MAX = 10;

            double ak;
            double ak1;
            double akn;
            double[] al =
            {
                0.4142135623730950,
                0.7320508075688773,
                0.2360679774997897,
                0.6457513110645906,
                0.3166247903553998,
                0.6055512754639893,
                0.1231056256176605,
                0.3589989435406736,
                0.7958315233127195,
                0.3851648071345040
            };
            double b;
            double[] be;
            double bk;
            double d1;
            double d2;
            double[] dex;
            int dim;
            double g;
            double[] ga;
            int i;
            int j;
            int k;
            int key;
            bool more;
            double[] p1;
            double[] p2;
            double[] p3;
            double[] p4;
            double s1;
            double s2;
            double t;
            double y1;
            double y2;
            double y3;
            double y4;

            eval_num = 0;
            //
            //  Check input.
            //
            if (dim_num < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("SAMPLE_ND - Fatal error!");
                Console.WriteLine("  DIM_NUM must be at least 1,");
                Console.WriteLine("  but DIM_NUM = " + dim_num + "");
                return;
            }

            if (DIM_MAX < dim_num)
            {
                Console.WriteLine("");
                Console.WriteLine("SAMPLE_ND - Fatal error!");
                Console.WriteLine("  DIM_NUM must be no more than DIM_MAX = " + DIM_MAX + "");
                Console.WriteLine("  but DIM_NUM = " + dim_num + "");
                return;
            }

            if (k1 < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("SAMPLE_ND - Fatal error!");
                Console.WriteLine("  K1 must be at least 1, but K1 = " + k1 + "");
                return;
            }

            if (k2 < k1)
            {
                Console.WriteLine("");
                Console.WriteLine("SAMPLE_ND - Fatal error!");
                Console.WriteLine("  K1 may not be greater than K2, but ");
                Console.WriteLine("  K1 = " + k1 + "");
                Console.WriteLine("  K2 = " + k2 + "");
                return;
            }

            be = new double[dim_num];
            dex = new double[dim_num];
            ga = new double[dim_num];
            p1 = new double[dim_num];
            p2 = new double[dim_num];
            p3 = new double[dim_num];
            p4 = new double[dim_num];

            for (dim = 0; dim < dim_num; dim++)
            {
                be[dim] = al[dim];
            }

            for (dim = 0; dim < dim_num; dim++)
            {
                ga[dim] = al[dim];
            }

            for (dim = 0; dim < dim_num; dim++)
            {
                dex[dim] = 0.0;
            }

            for (k = k1; k <= k2; k++)
            {
                ak = (double) (k);
                key = 0;
                ak1 = ak - 1.1;
                s1 = 0.0;
                d1 = 0.0;
                s2 = 0.0;
                d2 = 0.0;
                akn = Math.Pow(ak, dim_num);
                t = Math.Sqrt(Math.Pow(ak, dim_num)) * ak;
                bk = 1.0 / ak;

                for (;;)
                {
                    key = key + 1;

                    if (key != 1)
                    {
                        key = key - 1;
                        more = false;
                        for (j = 0; j < dim_num; j++)
                        {
                            if (dex[j] <= ak1)
                            {
                                dex[j] = dex[j] + 1.0;
                                more = true;
                                break;
                            }

                            dex[j] = 0.0;
                        }

                        if (!more)
                        {
                            break;
                        }
                    }

                    for (i = 0; i < dim_num; i++)
                    {
                        b = be[i] + al[i];
                        if (1.0 < b)
                        {
                            b = b - 1.0;
                        }

                        g = ga[i] + b;
                        if (1.0 < g)
                        {
                            g = g - 1.0;
                        }

                        be[i] = b + al[i];
                        if (1.0 < be[i])
                        {
                            be[i] = be[i] - 1.0;
                        }

                        ga[i] = be[i] + g;
                        if (1.0 < ga[i])
                        {
                            ga[i] = ga[i] - 1.0;
                        }

                        p1[i] = (dex[i] + g) * bk;
                        p2[i] = (dex[i] + 1.0 - g) * bk;
                        p3[i] = (dex[i] + ga[i]) * bk;
                        p4[i] = (dex[i] + 1.0 - ga[i]) * bk;
                    }

                    y1 = func(dim_num, p1);
                    eval_num = eval_num + 1;
                    //
                    //  There may be an error in the next two lines,
                    //  but oddly enough, that is how the original reads
                    //
                    y3 = func(dim_num, p2);
                    eval_num = eval_num + 1;
                    y2 = func(dim_num, p3);
                    eval_num = eval_num + 1;
                    y4 = func(dim_num, p4);
                    eval_num = eval_num + 1;

                    s1 = s1 + y1 + y2;
                    d1 = d1 + (y1 - y2) * (y1 - y2);
                    s2 = s2 + y3 + y4;
                    d2 = d2 + (y1 + y3 - y2 - y4) * (y1 + y3 - y2 - y4);
                }

                est1[k - 1] = 0.5 * s1 / akn;
                err1[k - 1] = 1.5 * Math.Sqrt(d1) / akn;
                dev1[k - 1] = err1[k - 1] * t;
                est2[k - 1] = 0.25 * (s1 + s2) / akn;
                err2[k - 1] = 0.75 * Math.Sqrt(d2) / akn;
                dev2[k - 1] = err2[k - 1] * t * ak;
            }

        }
    }
}