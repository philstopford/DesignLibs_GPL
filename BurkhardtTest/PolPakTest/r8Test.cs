using System;
using System.Numerics;
using Burkardt.RandomNS;
using Burkardt.Types;

namespace PolPakTest
{
    public static class r8Test
    {
        public static void r8_agm_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_AGM_TEST tests R8_AGM.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    27 April 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a = 0;
            double b = 0;
            double fx = 0;
            double fx2 = 0;
            int n_data;

            Console.WriteLine("");
            Console.WriteLine("R8_AGM_TEST");
            Console.WriteLine("  R8_AGM computes the arithmetic geometric mean.");
            Console.WriteLine("");
            Console.WriteLine("           A           B         "
                              + "   AGM                       AGM               Diff");
            Console.WriteLine("                             "
                              + "      (Tabulated)             R8_AGM(A,B)");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                Burkardt.Values.AGM.agm_values(ref n_data, ref a, ref b, ref fx);

                if (n_data == 0)
                {
                    break;
                }

                fx2 = typeMethods.r8_agm(a, b);

                Console.WriteLine("  " + a.ToString("0.######").PadLeft(10)
                                       + "  " + b.ToString("0.######").PadLeft(10)
                                       + "  " + fx.ToString("0.################").PadLeft(24)
                                       + "  " + fx2.ToString("0.################").PadLeft(24)
                                       + "  " + Math.Abs(fx - fx2).ToString("0.######").PadLeft(10) + "");
            }

        }

        public static void r8_beta_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_BETA_TEST tests R8_BETA.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    01 January 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fxy = 0;
            double fxy2 = 0;
            int n_data = 0;
            double x = 0;
            double y = 0;

            Console.WriteLine("");
            Console.WriteLine("R8_BETA_TEST:");
            Console.WriteLine("  R8_BETA evaluates the Beta function.");
            Console.WriteLine("");
            Console.WriteLine("     X      Y        Exact F       R8_BETA(X,Y)");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                Burkardt.Values.Beta.beta_values(ref n_data, ref x, ref y, ref fxy);

                if (n_data == 0)
                {
                    break;
                }

                fxy2 = typeMethods.r8_beta(x, y);

                Console.WriteLine("  "
                                  + x.ToString().PadLeft(10) + "  "
                                  + y.ToString().PadLeft(10) + "  "
                                  + fxy.ToString().PadLeft(10) + "  "
                                  + fxy2.ToString().PadLeft(10) + "");
            }

        }

        public static void r8_choose_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_CHOOSE_TEST tests R8_CHOOSE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 June 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double cnk = 0;
            int k = 0;
            int n = 0;

            Console.WriteLine("");
            Console.WriteLine("R8_CHOOSE_TEST");
            Console.WriteLine("  R8_CHOOSE evaluates C(N,K).");
            Console.WriteLine("");
            Console.WriteLine("   N     K    CNK");
            Console.WriteLine("");

            for (n = 0; n <= 4; n++)
            {
                for (k = 0; k <= n; k++)
                {
                    cnk = typeMethods.r8_choose(n, k);

                    Console.WriteLine("  "
                                      + n.ToString().PadLeft(6) + "  "
                                      + k.ToString().PadLeft(6) + "  "
                                      + cnk.ToString().PadLeft(6) + "");
                }
            }

        }

        public static void r8_cube_root_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    r8_cube_root_test() tests r8_cube_root().
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 May 2021
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i = 0;
            double x1 = 0;
            double y = 0;
            double x2 = 0;
            Rand48 rand48 = new Rand48();

            Console.WriteLine("");
            Console.WriteLine("r8_cube_root_test()");
            Console.WriteLine("  r8_cube_root() computes the cube root of an R8.");
            Console.WriteLine("");
            Console.WriteLine("       X               Y               Y^3");
            Console.WriteLine("");

            for (i = 1; i <= 10; i++)
            {
                x1 = -10.0 + 20.0 * rand48.Next();
                y = typeMethods.r8_cube_root(x1);
                x2 = Math.Pow(y, 3);
                Console.WriteLine(x1.ToString().PadLeft(14) + "  "
                                                            + y.ToString().PadLeft(14) + "  "
                                                            + x2.ToString().PadLeft(14) + "");
            }

        }

        public static void r8_erf_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_ERF_TEST tests R8_ERF.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 February 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fx = 0;
            double fx2 = 0;
            int n_data = 0;
            double x = 0;

            Console.WriteLine("");
            Console.WriteLine("R8_ERF_TEST:");
            Console.WriteLine("  R8_ERF evaluates the error function.");
            Console.WriteLine("");
            Console.WriteLine("     X      Exact F     R8_ERF(X)");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                Burkardt.Values.ErrorFunc.erf_values(ref n_data, ref x, ref fx);

                if (n_data == 0)
                {
                    break;
                }

                fx2 = typeMethods.r8_erf(x);

                Console.WriteLine("  "
                                  + x.ToString().PadLeft(8) + "  "
                                  + fx.ToString().PadLeft(14) + "  "
                                  + fx2.ToString().PadLeft(14) + "");
            }

        }

        public static void r8_erf_inverse_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_ERF_INVERSE_TEST tests R8_ERF_INVERSE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 August 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fx = 0;
            int n_data = 0;
            double x1 = 0;
            double x2 = 0;

            Console.WriteLine("");
            Console.WriteLine("R8_ERF_INVERSE_TEST");
            Console.WriteLine("  R8_ERF_INVERSE inverts the error function.");
            Console.WriteLine("");
            Console.WriteLine("    FX           X1           X2");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                Burkardt.Values.ErrorFunc.erf_values(ref n_data, ref x1, ref fx);

                if (n_data == 0)
                {
                    break;
                }

                x2 = typeMethods.r8_erf_inverse(fx);

                Console.WriteLine("  "
                                  + fx.ToString().PadLeft(8) + "  "
                                  + x1.ToString().PadLeft(14) + "  "
                                  + x2.ToString().PadLeft(14) + "");
            }

        }

        public static void r8_euler_constant_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_EULER_CONSTANT_TEST tests R8_EULER_CONSTANT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 January 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double g = 0;
            double g_approx = 0;
            int i = 0;
            int n = 0;
            double n_r8 = 0;
            int test = 0;

            g = typeMethods.r8_euler_constant();

            Console.WriteLine("");
            Console.WriteLine("R8_EULER_CONSTANT_TEST:");
            Console.WriteLine("  R8_EULER_CONSTANT returns the Euler-Mascheroni constant");
            Console.WriteLine("  sometimes denoted by 'gamma'.");
            Console.WriteLine("");
            Console.WriteLine("  gamma = limit ( N -> oo ) ( sum ( 1 <= I <= N ) 1 / I ) - log ( N )");
            Console.WriteLine("");
            Console.WriteLine("  Numerically, g = " + g + "");
            Console.WriteLine("");
            Console.WriteLine("         N      Partial Sum    |gamma - partial sum|");
            Console.WriteLine("");

            n = 1;
            for (test = 0; test <= 20; test++)
            {
                n_r8 = (double)(n);
                g_approx = -Math.Log(n_r8);
                for (i = 1; i <= n; i++)
                {
                    g_approx = g_approx + 1.0 / (double)(i);
                }

                Console.WriteLine("  " + n.ToString().PadLeft(8)
                                       + "  " + g_approx.ToString().PadLeft(14)
                                       + "  " + Math.Abs(g_approx - g).ToString().PadLeft(14) + "");
                n = n * 2;
            }

        }

        public static void r8_factorial_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_FACTORIAL_TEST tests R8_FACTORIAL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 June 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fn = 0;
            int n_data;
            int n = 0;

            Console.WriteLine("");
            Console.WriteLine("R8_FACTORIAL_TEST:");
            Console.WriteLine("  R8_FACTORIAL evaluates the factorial function.");
            Console.WriteLine("");
            Console.WriteLine("     N       Exact F       R8_FACTORIAL(N)");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                Burkardt.Values.Factorial.r8_factorial_values(ref n_data, ref n, ref fn);

                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + n.ToString().PadLeft(4) + "  "
                                  + fn.ToString().PadLeft(14) + "  "
                                  + typeMethods.r8_factorial(n).ToString().PadLeft(14) + "");
            }

        }

        public static void r8_factorial_log_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_FACTORIAL_LOG_TEST tests R8_FACTORIAL_LOG.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 June 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fn = 0;
            int n_data = 0;
            int n = 0;

            Console.WriteLine("");
            Console.WriteLine("R8_FACTORIAL_LOG_TEST:");
            Console.WriteLine("  R8_FACTORIAL_LOG evaluates the logarithm of the");
            Console.WriteLine("  factorial function.");
            Console.WriteLine("");
            Console.WriteLine("     N	   Exact F	 R8_FACTORIAL_LOG(N)");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                Burkardt.Values.Factorial.r8_factorial_log_values(ref n_data, ref n, ref fn);

                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + n.ToString().PadLeft(5) + "  "
                                  + fn.ToString().PadLeft(14) + "  "
                                  + typeMethods.r8_factorial_log(n).ToString().PadLeft(14) + "");

            }

        }

        public static void r8_gamma_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_GAMMA_TEST tests R8_GAMMA.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 April 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fx1 = 0;
            double fx2 = 0;
            int n_data;
            double x = 0;

            Console.WriteLine("");
            Console.WriteLine("R8_GAMMA_TEST:");
            Console.WriteLine("   R8_GAMMA evaluates the Gamma function.");
            Console.WriteLine("");
            Console.WriteLine("      X            GAMMA(X)     R8_GAMMA(X)");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                Burkardt.Values.Gamma.gamma_values(ref n_data, ref x, ref fx1);

                if (n_data == 0)
                {
                    break;
                }

                fx2 = typeMethods.r8_gamma(x);

                Console.WriteLine("  " + x.ToString().PadLeft(12)
                                       + "  " + fx1.ToString("0.################").PadLeft(24)
                                       + "  " + fx2.ToString("0.################").PadLeft(24) + "");
            }

        }

        public static void r8_hyper_2f1_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_HYPER_2F1_TEST tests R8_HYPER_2F1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 February 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a = 0;
            double b = 0;
            double c = 0;
            double fx = 0;
            double fx2 = 0;
            int n_data = 0;
            double x = 0;

            Console.WriteLine("");
            Console.WriteLine(" R8_HYPER_2F1_TEST:");
            Console.WriteLine("   R8_HYPER_2F1 evaluates the hypergeometric function 2F1.");
            Console.WriteLine("");
            Console.WriteLine("      A       B       C       X      " +
                              " 2F1                       2F1                     DIFF");
            Console.WriteLine("                                     " + "(tabulated)               (computed)");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                Burkardt.Values.Hypergeometric.hyper_2f1_values(ref n_data, ref a, ref b, ref c, ref x, ref fx);

                if (n_data == 0)
                {
                    break;
                }

                fx2 = typeMethods.r8_hyper_2f1(a, b, c, x);

                Console.WriteLine("  " + a.ToString("0.##").PadLeft(6)
                                       + "  " + b.ToString("0.##").PadLeft(6)
                                       + "  " + c.ToString("0.##").PadLeft(6)
                                       + "  " + x.ToString("0.##").PadLeft(6)
                                       + "  " + fx.ToString("0.################").PadLeft(24)
                                       + "  " + fx2.ToString("0.################").PadLeft(24)
                                       + "  " + Math.Abs(fx - fx2).ToString("0.####").PadLeft(10) + "");
            }
        }

        public static void r8_psi_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_PSI_TEST tests R8_PSI.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 February 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fx = 0;
            double fx2 = 0;
            int n_data = 0;
            double x = 0;

            Console.WriteLine("");
            Console.WriteLine("R8_PSI_TEST:");
            Console.WriteLine("  R8_PSI evaluates the Psi function.");
            Console.WriteLine("");
            Console.WriteLine("         X                  Psi(X)           "
                              + "         Psi(X)          DIFF");
            Console.WriteLine("                         (Tabulated)         "
                              + "       (R8_PSI)");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                Burkardt.Values.Psi.psi_values(ref n_data, ref x, ref fx);

                if (n_data == 0)
                {
                    break;
                }

                fx2 = typeMethods.r8_psi(x);

                Console.WriteLine("  " + x.ToString("0.##").PadLeft(8)
                                       + "  " + fx.ToString("0.################").PadLeft(24)
                                       + "  " + fx2.ToString("0.################").PadLeft(24)
                                       + "  " + Math.Abs(fx - fx2).ToString("0.####").PadLeft(10) + "");

            }

        }

        public static void r8poly_degree_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8POLY_DEGREE_TEST tests R8POLY_DEGREE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 January 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] c1 = { 1.0, 2.0, 3.0, 4.0 };
            double[] c2 = { 1.0, 2.0, 3.0, 0.0 };
            double[] c3 = { 1.0, 2.0, 0.0, 4.0 };
            double[] c4 = { 1.0, 0.0, 0.0, 0.0 };
            double[] c5 = { 0.0, 0.0, 0.0, 0.0 };
            int d;
            int m;

            Console.WriteLine("");
            Console.WriteLine("R8POLY_DEGREE_TEST");
            Console.WriteLine("  R8POLY_DEGREE determines the degree of an R8POLY.");

            m = 3;

            typeMethods.r8poly_print(m, c1, "  The R8POLY:");
            d = typeMethods.r8poly_degree(m, ref c1);
            Console.WriteLine("  Dimensioned degree = " + m + ",  Actual degree = " + d + "");

            typeMethods.r8poly_print(m, c2, "  The R8POLY:");
            d = typeMethods.r8poly_degree(m, ref c2);
            Console.WriteLine("  Dimensioned degree = " + m + ",  Actual degree = " + d + "");

            typeMethods.r8poly_print(m, c3, "  The R8POLY:");
            d = typeMethods.r8poly_degree(m, ref c3);
            Console.WriteLine("  Dimensioned degree = " + m + ",  Actual degree = " + d + "");

            typeMethods.r8poly_print(m, c4, "  The R8POLY:");
            d = typeMethods.r8poly_degree(m, ref c4);
            Console.WriteLine("  Dimensioned degree = " + m + ",  Actual degree = " + d + "");

            typeMethods.r8poly_print(m, c5, "  The R8POLY:");
            d = typeMethods.r8poly_degree(m, ref c5);
            Console.WriteLine("  Dimensioned degree = " + m + ",  Actual degree = " + d + "");

        }

        public static void r8poly_print_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8POLY_PRINT_TEST tests R8POLY_PRINT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 January 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] c = { 2.0, -3.4, 56.0, 0.0, 0.78, 9.0 };
            int m = 5;

            Console.WriteLine("");
            Console.WriteLine("R8POLY_PRINT_TEST");
            Console.WriteLine("  R8POLY_PRINT prints an R8POLY.");

            typeMethods.r8poly_print(m, c, "  The R8POLY:");

        }

        public static void r8poly_value_horner_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8POLY_VALUE_HORNER_TEST tests R8POLY_VALUE_HORNER.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    02 January 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] c = { 24.0, -50.0, +35.0, -10.0, 1.0 };
            int i;
            int m = 4;
            int n = 16;
            double p;
            double[] x;
            double x_hi;
            double x_lo;

            Console.WriteLine("");
            Console.WriteLine("R8POLY_VALUE_HORNER_TEST");
            Console.WriteLine("  R8POLY_VALUE_HORNER evaluates a polynomial at");
            Console.WriteLine("  one point, using Horner's method.");

            typeMethods.r8poly_print(m, c, "  The polynomial coefficients:");

            x_lo = 0.0;
            x_hi = 5.0;
            x = typeMethods.r8vec_linspace_new(n, x_lo, x_hi);

            Console.WriteLine("");
            Console.WriteLine("   I    X    P(X)");
            Console.WriteLine("");

            for (i = 0; i < n; i++)
            {
                p = typeMethods.r8poly_value_horner(m, c, x[i]);
                Console.WriteLine("  " + i.ToString().PadLeft(2)
                                       + "  " + x[i].ToString().PadLeft(8)
                                       + "  " + p.ToString().PadLeft(14) + "");
            }

        }

    }
}