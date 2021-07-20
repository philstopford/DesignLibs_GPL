using System;
using Burkardt.Laguerre;
using Burkardt.PolynomialNS;
using Burkardt.Types;

namespace LaguerrePolynomialTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for LAGUERRE_POLYNOMIAL_TEST.
            //
            //  Discussion:
            //
            //    LAGUERRE_POLYNOMIAL_TEST tests the LAGUERRE_POLYNOMIAL library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 March 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double b;
            int e;
            int p;

            Console.WriteLine("");
            Console.WriteLine("LAGUERRE_POLYNOMIAL_TEST:");
            Console.WriteLine("  Test the LAGUERRE_POLYNOMIAL library.");

            laguerre_polynomial_test01();
            laguerre_polynomial_test02();
            laguerre_polynomial_test03();
            laguerre_polynomial_test04();
            laguerre_polynomial_test05();
            laguerre_polynomial_test06();

            p = 5;
            b = 0.0;
            laguerre_polynomial_test07(p, b);

            p = 5;
            b = 1.0;
            laguerre_polynomial_test07(p, b);

            p = 5;
            e = 0;
            laguerre_polynomial_test08(p, e);

            p = 5;
            e = 1;
            laguerre_polynomial_test08(p, e);

            Console.WriteLine("");
            Console.WriteLine("LAGUERRE_POLYNOMIAL_TEST:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        public static void laguerre_polynomial_test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LAGUERRE_POLYNOMIAL_TEST01 tests L_POLYNOMIAL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 March 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int n_data;
            double e;
            double fx1 = 0;
            double fx2;
            double[] fx2_vec;
            int n = 0;
            double x = 0;
            double[] x_vec = new double[1];

            Console.WriteLine("");
            Console.WriteLine("LAGUERRE_POLYNOMIAL_TEST01:");
            Console.WriteLine("  L_POLYNOMIAL_VALUES stores values of");
            Console.WriteLine("  the Laguerre polynomials.");
            Console.WriteLine("  L_POLYNOMIAL evaluates the polynomial.");
            Console.WriteLine("");
            Console.WriteLine("                        Tabulated                 Computed");
            Console.WriteLine("     N        X           L(N,X)                    L(N,X)                     Error");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                Laguerre.l_polynomial_values(ref n_data, ref n, ref x, ref fx1);

                if (n_data == 0)
                {
                    break;
                }

                x_vec[0] = x;
                fx2_vec = Laguerre.l_polynomial(1, n, x_vec);
                fx2 = fx2_vec[n];

                e = fx1 - fx2;

                Console.WriteLine("  " + n.ToString().PadLeft(4)
                                       + "  " + x.ToString().PadLeft(12)
                                       + "  " + fx1.ToString("0.################").PadLeft(24)
                                       + "  " + fx2.ToString("0.################").PadLeft(24)
                                       + "  " + e.ToString().PadLeft(8) + "");
            }
        }

        public static void laguerre_polynomial_test02()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LAGUERRE_POLYNOMIAL_TEST02 tests L_POLYNOMIAL_COEFFICIENTS.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    11 March 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 10;

            double[] c;
            int i;
            int j;

            Console.WriteLine("");
            Console.WriteLine("LAGUERRE_POLYNOMIAL_TEST02");
            Console.WriteLine("  L_POLYNOMIAL_COEFFICIENTS determines Laguerre polynomial coefficients.");

            c = Laguerre.l_polynomial_coefficients(N);

            for (i = 0; i <= N; i++)
            {
                Console.WriteLine("");
                Console.WriteLine("  L(" + i + ",x) =");
                Console.WriteLine("");
                for (j = i; 0 <= j; j--)
                {
                    if (c[i + j * (N + 1)] == 0.0)
                    {
                    }
                    else if (j == 0)
                    {
                        Console.WriteLine(c[i + j * (N + 1)].ToString().PadLeft(14) + "");
                        ;
                    }
                    else if (j == 1)
                    {
                        Console.WriteLine(c[i + j * (N + 1)].ToString().PadLeft(14) + " * x");
                    }
                    else
                    {
                        Console.WriteLine(c[i + j * (N + 1)].ToString().PadLeft(14) + " * x^" + j + "");
                    }
                }
            }
        }

        public static void laguerre_polynomial_test03()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LAGUERRE_POLYNOMIAL_TEST03 tests L_POLYNOMIAL_ZEROS.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 March 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int degree;
            double[] hz;
            string title;
            double[] z;

            Console.WriteLine("");
            Console.WriteLine("LAGUERRE_POLYNOMIAL_TEST03:");
            Console.WriteLine("  L_POLYNOMIAL_ZEROS computes the zeros of L(n,x)");
            Console.WriteLine("  Check by calling L_POLYNOMIAL there.");

            for (degree = 1; degree <= 5; degree++)
            {
                z = Laguerre.l_polynomial_zeros(degree);
                title = "  Computed zeros for L(" + degree + ",z):";
                typeMethods.r8vec_print(degree, z, title);

                hz = Laguerre.l_polynomial(degree, degree, z);
                title = "  Evaluate L(" + degree + ",z):";
                typeMethods.r8vec_print(degree, hz, title, aIndex: +degree * degree);
            }
        }

        public static void laguerre_polynomial_test04()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LAGUERRE_POLYNOMIAL_TEST04 tests L_QUADRATURE_RULE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 March 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int e;
            double[] f;
            int i;
            int n;
            double q;
            double q_exact;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("LAGUERRE_POLYNOMIAL_TEST04:");
            Console.WriteLine("  L_QUADRATURE_RULE computes the quadrature rule");
            Console.WriteLine("  associated with L(n,x)");

            n = 7;
            x = new double[n];
            w = new double[n];

            QuadratureRule.l_quadrature_rule(n, ref x, ref w);

            typeMethods.r8vec2_print(n, x, w, "      X            W");

            Console.WriteLine("");
            Console.WriteLine("  Use the quadrature rule to estimate:");
            Console.WriteLine("");
            Console.WriteLine("    Q = Integral ( 0 <= X < +00 ) X^E exp(-X) dx");
            Console.WriteLine("");
            Console.WriteLine("   E       Q_Estimate      Q_Exact");
            Console.WriteLine("");

            f = new double[n];

            for (e = 0; e <= 2 * n - 1; e++)
            {
                if (e == 0)
                {
                    for (i = 0; i < n; i++)
                    {
                        f[i] = 1.0;
                    }
                }
                else
                {
                    for (i = 0; i < n; i++)
                    {
                        f[i] = Math.Pow(x[i], e);
                    }
                }

                q = typeMethods.r8vec_dot_product(n, w, f);
                q_exact = Laguerre.l_integral(e);
                Console.WriteLine("  " + e.ToString().PadLeft(2)
                                       + "  " + q.ToString().PadLeft(14)
                                       + "  " + q_exact.ToString().PadLeft(14) + "");
            }
        }

        public static void laguerre_polynomial_test05()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LAGUERRE_POLYNOMIAL_TEST05 tests LM_POLYNOMIAL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 March 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int n_data;
            double e;
            double fx1 = 0;
            double fx2;
            double[] fx2_vec;
            int m = 0;
            int n = 0;
            double x = 0;
            double[] x_vec = new double[1];

            Console.WriteLine("");
            Console.WriteLine("LAGUERRE_POLYNOMIAL_TEST05:");
            Console.WriteLine("  LM_POLYNOMIAL_VALUES stores values of");
            Console.WriteLine("  the Laguerre polynomials.");
            Console.WriteLine("  LM_POLYNOMIAL evaluates the polynomial.");
            Console.WriteLine("");
            Console.WriteLine("                              Tabulated                 Computed");
            Console.WriteLine(
                "     N     M        X         Lm(N,M,X)                  Lm(N,M,X)                     Error");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                Functions.lm_polynomial_values(ref n_data, ref n, ref m, ref x, ref fx1);

                if (n_data == 0)
                {
                    break;
                }

                x_vec[0] = x;
                fx2_vec = Functions.lm_polynomial(1, n, m, x_vec);
                fx2 = fx2_vec[n];

                e = fx1 - fx2;

                Console.WriteLine("  " + n.ToString().PadLeft(4)
                                       + "  " + m.ToString().PadLeft(4)
                                       + "  " + x.ToString().PadLeft(12)
                                       + "  " + fx1.ToString("0.################").PadLeft(24)
                                       + "  " + fx2.ToString("0.################").PadLeft(24)
                                       + "  " + e.ToString().PadLeft(8) + "");
            }
        }

        public static void laguerre_polynomial_test06()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LAGUERRE_POLYNOMIAL_TEST06 tests LM_POLYNOMIAL_COEFFICIENTS.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    11 March 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 5;

            double[] c;
            int i;
            int j;
            int m;

            Console.WriteLine("");
            Console.WriteLine("LAGUERRE_POLYNOMIAL_TEST06");
            Console.WriteLine("  LM_POLYNOMIAL_COEFFICIENTS determines Laguerre polynomial coefficients.");

            for (m = 0; m <= 4; m++)
            {
                c = Functions.lm_polynomial_coefficients(N, m);

                for (i = 0; i <= N; i++)
                {
                    Console.WriteLine("");
                    Console.WriteLine("  Lm(" + i + "," + m + ",x) =");
                    Console.WriteLine("");
                    for (j = i; 0 <= j; j--)
                    {
                        if (c[i + j * (N + 1)] == 0.0)
                        {
                        }
                        else if (j == 0)
                        {
                            Console.WriteLine(c[i + j * (N + 1)].ToString().PadLeft(14) + "");
                            ;
                        }
                        else if (j == 1)
                        {
                            Console.WriteLine(c[i + j * (N + 1)].ToString().PadLeft(14) + " * x");
                        }
                        else
                        {
                            Console.WriteLine(c[i + j * (N + 1)].ToString().PadLeft(14) + " * x^" + j + "");
                        }
                    }
                }
            }
        }

        public static void laguerre_polynomial_test07(int p, double b)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LAGUERRE_POLYNOMIAL_TEST07 tests L_EXPONENTIAL_PRODUCT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 March 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int P, the maximum degree of the polynomial 
            //    factors.
            //
            //    Input, double B, the coefficient of X in the exponential factor.
            //
        {
            double[] table;

            Console.WriteLine("");
            Console.WriteLine("LAGUERREE_POLYNOMIAL_TEST07");
            Console.WriteLine("  Compute an exponential product table for L(n,x):");
            Console.WriteLine("");
            Console.WriteLine("  Tij = integral ( 0 <= x < +oo ) exp(b*x) Ln(i,x) Ln(j,x) exp(-x) dx");
            Console.WriteLine("");
            Console.WriteLine("  Maximum degree P = " + p + "");
            Console.WriteLine("  Exponential argument coefficient B = " + b + "");

            table = Laguerre.l_exponential_product(p, b);

            typeMethods.r8mat_print(p + 1, p + 1, table, "  Exponential product table:");
        }

        public static void laguerre_polynomial_test08(int p, int e)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LAGUERRE_POLYNOMIAL_TEST08 tests L_POWER_PRODUCT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 March 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int P, the maximum degree of the polynomial 
            //    factors.
            //
            //    Input, int E, the exponent of X.
            //
        {
            double[] table;

            Console.WriteLine("");
            Console.WriteLine("LAGUERRE_POLYNOMIAL_TEST08");
            Console.WriteLine("  Compute a power product table for L(n,x).");
            Console.WriteLine("");
            Console.WriteLine("  Tij = integral ( 0 <= x < +oo ) x^e L(i,x) L(j,x) exp(-x) dx");

            Console.WriteLine("");
            Console.WriteLine("  Maximum degree P = " + p + "");
            Console.WriteLine("  Exponent of X, E = " + e + "");

            table = Laguerre.l_power_product(p, e);

            typeMethods.r8mat_print(p + 1, p + 1, table, "  Power product table:");
        }
    }
}