using System;
using Burkardt.IntegralNS;

namespace nIntTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************8080
            //
            //  Purpose:
            //
            //    MAIN is the main program for NINTLIB_TEST.
            //
            //  Discussion:
            //
            //    NINTLIB_TEST tests the NINTLIB library.
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
        {
            int TEST_NUM = 3;

            double a;
            double b;
            int dim_num;
            int[] dim_num_test = {2, 3, 4};
            int test;

            Console.WriteLine("");
            Console.WriteLine("NINTLIB_TEST");
            Console.WriteLine("  Test the NINTLIB library.");

            a = 0.0;
            b = 1.0;

            Console.WriteLine("");
            Console.WriteLine("TESTND");
            Console.WriteLine("  Test routines for estimating the integral of");
            Console.WriteLine("  of F(X) in the hypercube [A,B]**DIM_NUM.");
            Console.WriteLine("");

            for (test = 0; test < TEST_NUM; test++)
            {
                dim_num = dim_num_test[test];

                Console.WriteLine("");
                Console.WriteLine("");
                Console.WriteLine("  DIM_NUM = " + dim_num + "");
                Console.WriteLine("");
                Console.WriteLine("");
                Console.WriteLine("  A(1:DIM_NUM) = " + a + "");
                Console.WriteLine("  B(1:DIM_NUM) = " + b + "");

                Console.WriteLine("");
                Console.WriteLine("");
                Console.WriteLine("  F(X(1:DIM_NUM)) = 1");
                Console.WriteLine("");

                testnd(dim_num, f1dn);

                Console.WriteLine("");
                Console.WriteLine("");
                Console.WriteLine("  F(X(1:DIM_NUM)) = sum ( X(1:DIM_NUM) )");
                Console.WriteLine("");

                testnd(dim_num, fxdn);

                Console.WriteLine("");
                Console.WriteLine("");
                Console.WriteLine("  F(X(1:DIM_NUM)) = sum( X(1:DIM_NUM)^2 )");
                Console.WriteLine("");

                testnd(dim_num, fx2dn);

                Console.WriteLine("");
                Console.WriteLine("");
                Console.WriteLine("  F(X(1:DIM_NUM)) = sum ( X(1:DIM_NUM)^3 )");
                Console.WriteLine("");

                testnd(dim_num, fx3dn);

                Console.WriteLine("");
                Console.WriteLine("");
                Console.WriteLine("  F(X(1:DIM_NUM)) = exp(sum(X(1:DIM_NUM)))");
                Console.WriteLine("");

                testnd(dim_num, fedn);

                Console.WriteLine("");
                Console.WriteLine("");
                Console.WriteLine("  F(X(1:DIM_NUM)) = 1/(1+sum(X(1:DIM_NUM)^2))");
                Console.WriteLine("");

                testnd(dim_num, fbdn);
            }

            //
            //  Terminate.
            //
            Console.WriteLine("");
            Console.WriteLine("NINTLIB_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void testnd(int dim_num, Func<int, double[], double> func)

            //****************************************************************************8080
            //
            //  Purpose:
            //
            //    TESTND tests the integrators on a particular function.
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
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, double FUNC ( int dim_num, double x[] ), evaluates
            //    the function to be integrated.
            //
        {
            test01(dim_num, func);
            test02(dim_num, func);
            test03(dim_num, func);
            test04(dim_num, func);
            if (dim_num == 2)
            {
                test05(dim_num, func);
            }

            test06(dim_num, func);
        }

        static void test01(int dim_num, Func<int, double[], double> func)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 tests BOX_ND.
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
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, double FUNC ( int dim_num, double x[] ), evaluates the function
            //    to be integrated.
            //
        {
            int ORDER = 5;

            int eval_num = 0;
            int i;
            double result;
            double[] wtab =
            {
                0.236926885056189087514264040720,
                0.478628670499366468041291514836,
                0.568888888888888888888888888889,
                0.478628670499366468041291514836,
                0.236926885056189087514264040720
            };
            double[] wtab2 = new double[ORDER];
            double[] xtab =
            {
                -0.906179845938663992797626878299,
                -0.538469310105683091036314420700,
                0.0,
                0.538469310105683091036314420700,
                0.906179845938663992797626878299
            };
            double[] xtab2 = new double[ORDER];
            //
            //  Adjust the quadrature rule from [-1,1] to [0,1]:
            //
            for (i = 0; i < ORDER; i++)
            {
                xtab2[i] = (xtab[i] + 1.0) / 2.0;
            }

            for (i = 0; i < ORDER; i++)
            {
                wtab2[i] = 0.5 * wtab[i];
            }

            result = nInt.box_nd(func, dim_num, ORDER, xtab2, wtab2, ref eval_num);

            Console.WriteLine("  BOX_ND:         "
                              + result.ToString("0.############").PadLeft(20)
                              + "  " + eval_num.ToString().PadLeft(8) + "");

        }

        static void test02(int dim_num, Func<int, double[], double> func)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST02 tests P5_ND.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, double FUNC ( int dim_num, double x[] ), evaluates the function
            //    to be integrated.
            //
        {
            double[] a;
            double[] b;
            int dim;
            int eval_num = 0;
            double result;
            //
            //  Set the integration limits.
            //
            a = new double[dim_num];
            b = new double[dim_num];

            for (dim = 0; dim < dim_num; dim++)
            {
                a[dim] = 0.0;
            }

            for (dim = 0; dim < dim_num; dim++)
            {
                b[dim] = 1.0;
            }

            result = nInt.p5_nd(func, dim_num, a, b, ref eval_num);

            Console.WriteLine("  P5_ND:          "
                              + result.ToString("0.############").PadLeft(20)
                              + "  " + eval_num.ToString().PadLeft(8) + "");
        }

        static void test03(int dim_num, Func<int, double[], double> func)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST03 tests ROMBERG_ND.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, double FUNC ( int dim_num, double x[] ), evaluates the function
            //    to be integrated.
            //
        {
            double[] a;
            double[] b;
            int dim;
            int eval_num = 0;
            int ind = 0;
            int it_max = 3;
            double result;
            int[] sub_num;
            double tol;
            //
            //  Set the integration limits.
            //
            a = new double[dim_num];
            b = new double[dim_num];
            sub_num = new int[dim_num];

            for (dim = 0; dim < dim_num; dim++)
            {
                a[dim] = 0.0;
            }

            for (dim = 0; dim < dim_num; dim++)
            {
                b[dim] = 1.0;
            }

            tol = 0.001;
            for (dim = 0; dim < dim_num; dim++)
            {
                sub_num[dim] = 10;
            }

            result = nInt.romberg_nd(func, a, b, dim_num, sub_num, it_max, tol,
                ref ind, ref eval_num);

            Console.WriteLine("  ROMBERG_ND:     "
                              + result.ToString("0.############").PadLeft(20)
                              + "  " + eval_num.ToString().PadLeft(8) + "");
        }

        static void test04(int dim_num, Func<int, double[], double> func)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST04 tests SAMPLE_ND.
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
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, double FUNC ( int dim_num, double x[] ), evaluates the function
            //    to be integrated.
            //
        {
            int K2 = 4;

            double[] dev1 = new double[K2];
            double[] dev2 = new double[K2];
            double[] err1 = new double[K2];
            double[] est1 = new double[K2];
            double[] est2 = new double[K2];
            double[] err2 = new double[K2];
            int eval_num = 0;
            int k1;

            k1 = 1;

            nInt.sample_nd(func, k1, K2, dim_num, est1, err1, dev1, est2, err2,
                dev2, ref eval_num);

            Console.WriteLine("  SAMPLE_ND:      "
                              + est2[K2 - 1].ToString("0.############").PadLeft(20)
                              + "  " + eval_num.ToString().PadLeft(8) + "");

        }

        static void test05(int dim_num, Func<int, double[], double> func)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST05 demonstrates how to refine multi-dimensional integration results.
            //
            //  Discussion:
            //
            //    This routine is only set up for DIM_NUM = 2 for now.
            //
            //    We are given a routine, NDP5, which will integrate over a
            //    DIM_NUM dimensional hypercube using a fixed method.  In order to
            //    improve the approximation to an integral, we can subdivide
            //    the hypercube and call NDP5 to integrate again over each of
            //    these regions.
            //
            //    The information that we gather can be used to tell us when
            //    to expect that we have achieved a certain degree of accuracy.
            //
            //    With a little more work, we could make this code adaptive.
            //    That is, it would only refine SOME of the subregions, where
            //    the approximation to the integral was still not good enough.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, integer DIM_NUM, the spatial dimension.
            //
            //    Input, double FUNC ( int dim_num, double x[] ), evaluates the function
            //    to be integrated.
            //
        {
            double[] a;
            double[] b;
            int dim;
            int eval_num = 0;
            int eval_total;
            int i;
            int igrid;
            int j;
            int ngrid;
            double result;
            double result_total;
            double[] xlo;
            double[] xhi;

            a = new double[dim_num];
            b = new double[dim_num];
            xlo = new double[dim_num];
            xhi = new double[dim_num];

            for (dim = 0; dim < dim_num; dim++)
            {
                xlo[dim] = 0.0;
            }

            for (dim = 0; dim < dim_num; dim++)
            {
                xhi[dim] = 1.0;
            }

            for (igrid = 1; igrid <= 6; igrid++)
            {
                ngrid = (int) Math.Pow(2, igrid - 1);

                result_total = 0.0;
                eval_total = 0;

                for (i = 1; i <= ngrid; i++)
                {
                    a[0] = ((double) (ngrid - i + 1) * xlo[0]
                            + (double) (i - 1) * xhi[0])
                           / (double) (ngrid);

                    b[0] = ((double) (ngrid - i) * xlo[0]
                            + (double) (i) * xhi[0])
                           / (double) (ngrid);

                    for (j = 1; j <= ngrid; j++)
                    {
                        a[1] = ((double) (ngrid - j + 1) * xlo[1]
                                + (double) (j - 1) * xhi[1])
                               / (double) (ngrid);

                        b[1] = ((double) (ngrid - j) * xlo[1]
                                + (double) (j) * xhi[1])
                               / (double) (ngrid);

                        result = nInt.p5_nd(func, dim_num, a, b, ref eval_num);

                        result_total = result_total + result;
                        eval_total = eval_total + eval_num;
                    }
                }

                Console.WriteLine("  P5_ND+:         "
                                  + result_total.ToString("0.############").PadLeft(20)
                                  + "  " + eval_total.ToString().PadLeft(8) + "");

            }
        }

        static void test06(int dim_num, Func<int, double[], double> func)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST06 tests MONTE_CARLO_ND.
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
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, double FUNC ( int dim_num, double x[] ), evaluates the function
            //    to be integrated.
            //
        {
            double[] a;
            double[] b;
            int dim;
            int eval_num;
            double result;
            int seed;
            int test;
            int test_num = 3;

            seed = 123456789;
            //
            //  Set the integration limits.
            //
            a = new double[dim_num];
            b = new double[dim_num];

            for (dim = 0; dim < dim_num; dim++)
            {
                a[dim] = 0.0;
            }

            for (dim = 0; dim < dim_num; dim++)
            {
                b[dim] = 1.0;
            }

            for (test = 1; test <= test_num; test++)
            {
                eval_num = (int) Math.Pow(8, test) * 10000;

                result = nInt.monte_carlo_nd(func, dim_num, a, b, eval_num, ref seed);

                Console.WriteLine("  MONTE_CARLO_ND: "
                                  + result.ToString("0.############").PadLeft(20)
                                  + "  " + eval_num.ToString().PadLeft(8) + "");
            }
        }

        static double fbdn(int dim_num, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FBDN(X(1:DIM_NUM)) = 1 / ( 1 + sum ( X(1:DIM_NUM)**2 ) )
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
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, double X[DIM_NUM], the argument.
            //
            //    Output, double FBDN, the value of the function at X.
            //
        {
            double arg;
            int dim;
            double value;

            arg = 0.0;
            for (dim = 0; dim < dim_num; dim++)
            {
                arg = arg + x[dim] * x[dim];
            }

            value = 1.0 / (1.0 + arg);

            return value;
        }

        static double fedn(int dim_num, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FEDN(X(1:DIM_NUM)) = EXP ( sum ( X(1:DIM_NUM) ) )
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
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, double X[DIM_NUM], the argument.
            //
            //    Output, double FEDN, the value of the function at X.
            //
        {
            double arg;
            int dim;
            double value;

            arg = 0.0;
            for (dim = 0; dim < dim_num; dim++)
            {
                arg = arg + x[dim];
            }

            value = Math.Exp(arg);

            return value;
        }

        static double f1dn(int dim_num, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    F1DN(X(1:DIM_NUM)) = 1.
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
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, double X[DIM_NUM], the argument.
            //
            //    Output, double F1DN, the value of the function at X.
            //
        {
            double value;

            value = 1.0;

            return value;
        }

        static double fxdn(int dim_num, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FXDN(X(1:DIM_NUM)) = sum ( X(1:DIM_NUM) )
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
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, double X[DIM_NUM], the argument.
            //
            //    Output, double FXDN, the value of the function at X.
            //
        {
            double arg;
            int dim;
            double value;

            arg = 0.0;
            for (dim = 0; dim < dim_num; dim++)
            {
                arg = arg + x[dim];
            }

            value = arg;

            return value;
        }

        static double fx2dn(int dim_num, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FX2DN(X(1:DIM_NUM)) = sum ( X(1:DIM_NUM)**2 )
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
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, double X[DIM_NUM], the argument.
            //
            //    Output, double FX2DN, the value of the function at X.
            //
        {
            double arg;
            int dim;
            double value;

            arg = 0.0;
            for (dim = 0; dim < dim_num; dim++)
            {
                arg = arg + x[dim] * x[dim];
            }

            value = arg;

            return value;
        }

        static double fx3dn(int dim_num, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FX3DN(X(1:DIM_NUM)) = sum ( X(1:DIM_NUM)**3 )
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
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, double X[DIM_NUM], the argument.
            //
            //    Output, double FX3DN, the value of the function at X.
            //
        {
            double arg;
            int dim;
            double value;

            arg = 0.0;
            for (dim = 0; dim < dim_num; dim++)
            {
                arg = arg + x[dim] * x[dim] * x[dim];
            }

            value = arg;

            return value;
        }
    }
}