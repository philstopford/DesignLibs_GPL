using System;
using Burkardt;
using Burkardt.IntegralNS;
using Burkardt.MatrixNS;
using Burkardt.PolynomialNS;
using Burkardt.Probability;
using Burkardt.Quadrature;
using Burkardt.Types;
using Burkardt.Uniform;

namespace GegenbauerPolynomialTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GEGENBAUER_POLYNOMIAL_TEST tests the GEGENBAUER_POLYNOMIAL library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 November 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("GEGENBAUER_POLYNOMIAL_TEST:");
            Console.WriteLine("  Test the GEGENBAUER_POLYNOMIAL library.");

            gegenbauer_alpha_check_test();
            gegenbauer_ek_compute_test();
            gegenbauer_integral_test();
            gegenbauer_polynomial_value_test();
            gegenbauer_ss_compute_test();

            imtqlx_test();

            r8_hyper_2f1_test();
            r8_uniform_ab_test();

            Console.WriteLine("");
            Console.WriteLine("GEGENBAUER_POLYNOMIAL_TEST:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void gegenbauer_alpha_check_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GEGENBAUER_ALPHA_CHECK_TEST compares GEGENBAUER_ALPHA_CHECK.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 November 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double alpha;
            bool check;
            int n;
            int seed;

            Console.WriteLine("");
            Console.WriteLine("GEGENBAUER_ALPHA_CHECK_TEST");
            Console.WriteLine("  GEGENBAUER_ALPHA_CHECK checks that ALPHA is legal;");
            Console.WriteLine("");
            Console.WriteLine("       ALPHA   Check?");
            Console.WriteLine("");

            seed = 123456789;

            for (n = 1; n <= 10; n++)
            {
                alpha = UniformRNG.r8_uniform_ab(-5.0, +5.0, ref seed);
                check = GegenbauerPolynomial.gegenbauer_alpha_check(alpha);
                Console.WriteLine("  " + alpha.ToString().PadLeft(10)
                                       + "       " + check.ToString().PadLeft(1) + "");
            }

            return;
        }

        static void gegenbauer_ek_compute_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GEGENBAUER_EK_COMPUTE_TEST tests GEGENBAUER_EK_COMPUTE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 November 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double alpha;
            int i;
            int n;
            int prec;
            double[] w;
            double[] x;

            alpha = 0.5;

            Console.WriteLine("");
            Console.WriteLine("GEGENBAUER_EK_COMPUTE_TEST");
            Console.WriteLine("  GEGENBAUER_EK_COMPUTE computes a Gauss-Gegenbauer rule;");
            Console.WriteLine("");
            Console.WriteLine("  with ALPHA = " + alpha + "");
            Console.WriteLine("  and integration interval [-1,+1]");
            Console.WriteLine("");
            Console.WriteLine("                  W               X");

            for (n = 1; n <= 10; n++)
            {
                Console.WriteLine("");

                w = new double[n];
                x = new double[n];

                GegenbauerPolynomial.gegenbauer_ek_compute(n, alpha, ref x, ref w);

                for (i = 0; i < n; i++)
                {
                    Console.WriteLine("          "
                                      + "  " + w[i].ToString("0.################").PadLeft(14)
                                      + "  " + x[i].ToString("0.################").PadLeft(14) + "");
                }
            }
        }

        static void gegenbauer_integral_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GEGENBAUER_INTEGRAL_TEST tests GEGENBAUER_INTEGRAL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 June 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double alpha;
            int n;
            int prec;
            double value;

            alpha = 0.25;

            Console.WriteLine("");
            Console.WriteLine("GEGENBAUER_INTEGRAL_TEST");
            Console.WriteLine("  GEGENBAUER_INTEGRAL evaluates");
            Console.WriteLine("  Integral ( -1 < x < +1 ) x^n * (1-x*x)^alpha dx");
            Console.WriteLine("");
            Console.WriteLine("         N         Value");
            Console.WriteLine("");

            for (n = 0; n <= 10; n++)
            {
                value = Integral.gegenbauer_integral(n, alpha);
                Console.WriteLine("  " + n.ToString("0.################").PadLeft(8)
                                       + "  " + value.ToString("0.################").PadLeft(24) + "");
            }
        }

        static void gegenbauer_polynomial_value_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GEGENBAUER_POLYNOMIAL_VALUE_TEST tests GEGENBAUER_POLYNOMIAL_VALUE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 November 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double alpha = 0;
            double[] c;
            double fx = 0;
            int m = 0;
            int n = 0;
            int n_data;
            double[] x = new double[1];
            double xscalar = 0;

            Console.WriteLine("");
            Console.WriteLine("GEGENBAUER_POLYNOMIAL_VALUE_TEST:");
            Console.WriteLine("  GEGENBAUER_POLYNOMIAL_VALUE evaluates the Gegenbauer polynomial.");
            Console.WriteLine("");
            Console.WriteLine("       M     ALPHA         X           GPV    GEGENBAUER");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                Burkardt.TestValues.Gegenbauer.gegenbauer_poly_values(ref n_data, ref m, ref alpha, ref xscalar, ref fx);

                if (n_data == 0)
                {
                    break;
                }

                //
                //  Since GEGENBAUER_POLYNOMIAL_VALUE expects a vector input X, we have to
                //  do a little "rewrapping" of the input.
                //
                n = 1;
                x[0] = xscalar;
                c = GegenbauerPolynomial.gegenbauer_polynomial_value(m, n, alpha, x);
                Console.WriteLine("  " + m.ToString("0.################").PadLeft(6)
                                       + "  " + alpha.ToString("0.################").PadLeft(8)
                                       + "  " + x[0].ToString("0.################").PadLeft(8)
                                       + "  " + fx.ToString("0.################").PadLeft(12)
                                       + "  " + c[m + 0 * (m + 1)].ToString("0.################").PadLeft(12) + "");
            }
        }

        static void gegenbauer_ss_compute_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GEGENBAUER_SS_COMPUTE_TEST tests GEGENBAUER_SS_COMPUTE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    10 June 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double alpha;
            int i;
            int n;
            int prec;
            double[] w;
            double[] x;

            alpha = 0.5;

            Console.WriteLine("");
            Console.WriteLine("GEGENBAUER_SS_COMPUTE_TEST");
            Console.WriteLine("  GEGENBAUER_SS_COMPUTE computes a Gauss-Gegenbauer rule;");
            Console.WriteLine("");
            Console.WriteLine("  with ALPHA = " + alpha + "");
            Console.WriteLine("");
            Console.WriteLine("                  W               X");

            for (n = 1; n <= 10; n++)
            {
                Console.WriteLine("");

                w = new double[n];
                x = new double[n];

                GegenbauerQuadrature.gegenbauer_ss_compute(n, alpha, ref x, ref w);

                for (i = 0; i < n; i++)
                {
                    Console.WriteLine("          "
                                      + "  " + w[i].ToString("0.################").PadLeft(14)
                                      + "  " + x[i].ToString("0.################").PadLeft(14) + "");
                }
            }
        }

        static void imtqlx_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    IMTQLX_TEST tests IMTQLX.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    10 June 2015
            //
            //  Author:
            //
            //    John Burkardt.
            //
        {
            double angle;
            double[] d = new double[5];
            double[] e = new double[5];
            int i;
            double[] lam = new double[5];
            double[] lam2 = new double[5];
            int n = 5;
            double[] qtz = new double[5];
            double[] z = new double[5];

            Console.WriteLine("");
            Console.WriteLine("IMTQLX_TEST");
            Console.WriteLine("  IMTQLX takes a symmetric tridiagonal matrix A");
            Console.WriteLine("  and computes its eigenvalues LAM.");
            Console.WriteLine("  It also accepts a vector Z and computes Q'*Z,");
            Console.WriteLine("  where Q is the matrix that diagonalizes A.");

            for (i = 0; i < n; i++)
            {
                d[i] = 2.0;
            }

            for (i = 0; i < n - 1; i++)
            {
                e[i] = -1.0;
            }

            e[n - 1] = 0.0;
            for (i = 0; i < n; i++)
            {
                z[i] = 1.0;
            }

            //
            //  On input, LAM is D, and QTZ is Z.
            //
            for (i = 0; i < n; i++)
            {
                lam[i] = d[i];
            }

            for (i = 0; i < n; i++)
            {
                qtz[i] = z[i];
            }

            IMTQLX.imtqlx(n, ref lam, ref e, ref qtz);

            typeMethods.r8vec_print(n, lam, "  Computed eigenvalues:");

            for (i = 0; i < n; i++)
            {
                angle = (double) (i + 1) * Math.PI / (double) (2 * (n + 1));
                lam2[i] = 4.0 * Math.Pow(Math.Sin(angle), 2);
            }

            typeMethods.r8vec_print(n, lam2, "  Exact eigenvalues:");

            typeMethods.r8vec_print(n, z, "  Vector Z:");
            typeMethods.r8vec_print(n, qtz, "  Vector Q'*Z:");
        }

        static void r8_hyper_2f1_test()

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
            int n_data;
            double x = 0;

            Console.WriteLine("");
            Console.WriteLine(" R8_HYPER_2F1_TEST:");
            Console.WriteLine("   R8_HYPER_2F1 evaluates the hypergeometric function 2F1.");
            Console.WriteLine("");
            Console.WriteLine("      A       B       C       X      " +
                              " 2F1                       2F1                     DIFF");
            Console.WriteLine("                                     " +
                              "(tabulated)               (computed)");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                Burkardt.TestValues.Hypergeometric.hyper_2f1_values(ref n_data, ref a, ref b, ref c, ref x, ref fx);

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

        static void r8_uniform_ab_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_UNIFORM_AB_TEST tests R8_UNIFORM_AB.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 April 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a;
            double b;
            double c;
            int i;
            int seed;

            b = 10.0;
            c = 25.0;
            seed = 17;

            Console.WriteLine("");
            Console.WriteLine("R8_UNIFORM_AB_TEST");
            Console.WriteLine("  R8_UNIFORM_AB produces a random real in a given range.");
            Console.WriteLine("");
            Console.WriteLine("  Using range " + b + " <= A <= " + c + ".");
            Console.WriteLine("");

            Console.WriteLine("");
            Console.WriteLine("      I       A");
            Console.WriteLine("");
            for (i = 0; i < 10; i++)
            {
                a = UniformRNG.r8_uniform_ab(b, c, ref seed);
                Console.WriteLine(i.ToString().PadLeft(6) + " "
                                                          + a.ToString().PadLeft(10) + "");
            }
        }
    }
}