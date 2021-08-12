using System;
using Burkardt;
using Burkardt.IntegralNS;
using Burkardt.PolynomialNS;
using Burkardt.Quadrature;
using Burkardt.Types;

namespace HermitePolynomialTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for HERMITE_POLYNOMIAL_TEST.
            //
            //  Discussion:
            //
            //    HERMITE_POLYNOMIAL_TEST tests the HERMITE_POLYNOMIAL library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 March 2012
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
            Console.WriteLine("HERMITE_POLYNOMIAL_TEST:");
            Console.WriteLine("  Test the HERMITE_POLYNOMIAL library.");

            hermite_polynomial_test01();
            hermite_polynomial_test02();
            hermite_polynomial_test03();
            hermite_polynomial_test04();
            hermite_polynomial_test05();
            hermite_polynomial_test06();
            hermite_polynomial_test07();

            p = 5;
            b = 0.0;
            hermite_polynomial_test08(p, b);

            p = 5;
            b = 1.0;
            hermite_polynomial_test08(p, b);

            p = 5;
            e = 0;
            hermite_polynomial_test09(p, e);

            p = 5;
            e = 1;
            hermite_polynomial_test09(p, e);

            p = 5;
            b = 0.0;
            hermite_polynomial_test10(p, b);

            p = 5;
            b = 1.0;
            hermite_polynomial_test10(p, b);

            p = 5;
            e = 0;
            hermite_polynomial_test11(p, e);

            p = 5;
            e = 1;
            hermite_polynomial_test11(p, e);

            p = 5;
            b = 0.0;
            hermite_polynomial_test12(p, b);

            p = 5;
            b = 1.0;
            hermite_polynomial_test12(p, b);

            p = 5;
            e = 0;
            hermite_polynomial_test13(p, e);

            p = 5;
            e = 1;
            hermite_polynomial_test13(p, e);

            hermite_polynomial_test14();

            hermite_polynomial_test15();
            //
            //  Terminate.
            //
            Console.WriteLine("");
            Console.WriteLine("HERMITE_POLYNOMIAL_TEST:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void hermite_polynomial_test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HERMITE_POLYNOMIAL_TEST01 tests H_POLYNOMIAL_VALUE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    25 February 2012
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
            Console.WriteLine("HERMITE_POLYNOMIAL_TEST01:");
            Console.WriteLine("  H_POLYNOMIAL_VALUES stores values of");
            Console.WriteLine("  the physicist's Hermite polynomials.");
            Console.WriteLine("  H_POLYNOMIAL_VALUE evaluates the polynomial.");
            Console.WriteLine("");
            Console.WriteLine("                        Tabulated                 Computed");
            Console.WriteLine("     N        X           H(N,X)                    H(N,X)                     Error");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                Hermite.h_polynomial_values(ref n_data, ref n, ref x, ref fx1);

                if (n_data == 0)
                {
                    break;
                }

                x_vec[0] = x;
                fx2_vec = Hermite.h_polynomial_value(1, n, x_vec);
                fx2 = fx2_vec[n];

                e = fx1 - fx2;

                Console.WriteLine("  " + n.ToString().PadLeft(4)
                                       + "  " + x.ToString().PadLeft(12)
                                       + "  " + fx1.ToString("0.################").PadLeft(24)
                                       + "  " + fx2.ToString("0.################").PadLeft(24)
                                       + "  " + e.ToString().PadLeft(8) + "");
            }
        }

        static void hermite_polynomial_test02()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HERMITE_POLYNOMIAL_TEST02 tests HE_POLYNOMIAL_VALUE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 February 2012
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
            Console.WriteLine("HERMITE_POLYNOMIAL_TEST02:");
            Console.WriteLine("  HE_POLYNOMIAL_VALUES stores values of");
            Console.WriteLine("  the probabilist's Hermite polynomials.");
            Console.WriteLine("  HE_POLYNOMIAL_VALUE evaluates the polynomial.");
            Console.WriteLine("");
            Console.WriteLine("                        Tabulated                 Computed");
            Console.WriteLine("     N        X          He(N,X)                   He(N,X)                     Error");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                Hermite.he_polynomial_values(ref n_data, ref n, ref x, ref fx1);

                if (n_data == 0)
                {
                    break;
                }

                x_vec[0] = x;
                fx2_vec = Hermite.he_polynomial_value(1, n, x_vec);
                fx2 = fx2_vec[n];

                e = fx1 - fx2;

                Console.WriteLine("  " + n.ToString().PadLeft(4)
                                       + "  " + x.ToString().PadLeft(12)
                                       + "  " + fx1.ToString("0.################").PadLeft(24)
                                       + "  " + fx2.ToString("0.################").PadLeft(24)
                                       + "  " + e.ToString().PadLeft(8) + "");
            }
        }

        static void hermite_polynomial_test03()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HERMITE_POLYNOMIAL_TEST03 tests HF_FUNCTION_VALUE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 February 2012
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
            Console.WriteLine("HERMITE_POLYNOMIAL_TEST03:");
            Console.WriteLine("  HF_FUNCTION_VALUES stores values of");
            Console.WriteLine("  the Hermite function Hf(n,x).");
            Console.WriteLine("  HF_FUNCTION_VALUE evaluates the function.");
            Console.WriteLine("");
            Console.WriteLine("                        Tabulated                 Computed");
            Console.WriteLine("     N        X          Hf(N,X)                   Hf(N,X)                   Error");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                Hermite.hf_function_values(ref n_data, ref n, ref x, ref fx1);

                if (n_data == 0)
                {
                    break;
                }

                x_vec[0] = x;
                fx2_vec = Hermite.hf_function_value(1, n, x_vec);
                fx2 = fx2_vec[n];

                e = fx1 - fx2;

                Console.WriteLine("  " + n.ToString().PadLeft(4)
                                       + "  " + x.ToString().PadLeft(12)
                                       + "  " + fx1.ToString("0.################").PadLeft(24)
                                       + "  " + fx2.ToString("0.################").PadLeft(24)
                                       + "  " + e.ToString("0.######").PadLeft(14) + "");
            }
        }

        static void hermite_polynomial_test04()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HERMITE_POLYNOMIAL_TEST04 tests H_POLYNOMIAL_ZEROS.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    19 October 2014
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
            Console.WriteLine("HERMITE_POLYNOMIAL_TEST04:");
            Console.WriteLine("  H_POLYNOMIAL_ZEROS computes the zeros of H(n,x)");
            Console.WriteLine("  Check by calling H_POLYNOMIAL there.");

            for (degree = 1; degree <= 5; degree++)
            {
                z = Hermite.h_polynomial_zeros(degree);
                title = "  Computed zeros for H(" + (degree) + ",z):";
                typeMethods.r8vec_print(degree, z, title);

                hz = Hermite.h_polynomial_value(degree, degree, z);
                title = "  Evaluate H(" + (degree) + ",z):";
                typeMethods.r8vec_print(degree, hz, title, aIndex: +degree * degree);

            }
        }

        static void hermite_polynomial_test05()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HERMITE_POLYNOMIAL_TEST05 tests HE_POLYNOMIAL_ZEROS.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 February 2012
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
            Console.WriteLine("HERMITE_POLYNOMIAL_TEST05:");
            Console.WriteLine("  HE_POLYNOMIAL_ZEROS computes the zeros of He(n,x)");
            Console.WriteLine("  Check by calling HE_POLYNOMIAL there.");

            for (degree = 1; degree <= 5; degree++)
            {
                z = Hermite.he_polynomial_zeros(degree);
                title = "  Computed zeros for He(" + (degree) + ",z):";
                typeMethods.r8vec_print(degree, z, title);

                hz = Hermite.he_polynomial_value(degree, degree, z);
                title = "  Evaluate He(" + (degree) + ",z):";
                typeMethods.r8vec_print(degree, hz, title, aIndex: +degree * degree);
            }
        }

        static void hermite_polynomial_test06()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HERMITE_POLYNOMIAL_TEST06 tests H_QUADRATURE_RULE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    08 March 2012
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
            Console.WriteLine("HERMITE_POLYNOMIAL_TEST06:");
            Console.WriteLine("  H_QUADRATURE_RULE computes the quadrature rule");
            Console.WriteLine("  associated with H(n,x)");

            n = 7;
            x = new double[n];
            w = new double[n];

            HermiteQuadrature.h_quadrature_rule(n, ref x, ref w);

            typeMethods.r8vec2_print(n, x, w, "      X            W");

            Console.WriteLine("");
            Console.WriteLine("  Use the quadrature rule to estimate:");
            Console.WriteLine("");
            Console.WriteLine("    Q = Integral ( -oo < X < +00 ) X^E exp(-X^2) dx");
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
                q_exact = Integral.h_integral(e);
                Console.WriteLine("  " + e.ToString().PadLeft(2)
                                       + "  " + q.ToString().PadLeft(14)
                                       + "  " + q_exact.ToString().PadLeft(14) + "");
            }
        }

        static void hermite_polynomial_test07()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HERMITE_POLYNOMIAL_TEST07 tests HE_QUADRATURE_RULE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    08 March 2012
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
            Console.WriteLine("HERMITE_POLYNOMIAL_TEST07:");
            Console.WriteLine("  HE_QUADRATURE_RULE computes the quadrature rule");
            Console.WriteLine("  associated with He(n,x)");

            n = 7;
            x = new double[n];
            w = new double[n];

            HermiteQuadrature.he_quadrature_rule(n, ref x, ref w);

            typeMethods.r8vec2_print(n, x, w, "      X            W");

            Console.WriteLine("");
            Console.WriteLine("  Use the quadrature rule to estimate:");
            Console.WriteLine("");
            Console.WriteLine("    Q = Integral ( -oo < X < +00 ) X^E exp(-X^2) dx");
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
                q_exact = Integral.he_integral(e);
                Console.WriteLine("  " + e.ToString().PadLeft(2)
                                       + "  " + q.ToString().PadLeft(14)
                                       + "  " + q_exact.ToString().PadLeft(14) + "");
            }
        }

        static void hermite_polynomial_test08(int p, double b)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HERMITE_POLYNOMIAL_TEST08 tests HN_EXPONENTIAL_PRODUCT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 February 2012
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
            Console.WriteLine("HERMITE_POLYNOMIAL_TEST08");
            Console.WriteLine("  Compute a normalized physicist''s Hermite exponential product table.");
            Console.WriteLine("");
            Console.WriteLine("  Tij = integral ( -oo < X < +oo ) exp(B*X) Hn(I,X) Hn(J,X) exp(-X*X) dx");
            Console.WriteLine("");
            Console.WriteLine("  where Hn(I,X) = normalized physicist''s Hermite polynomial of degree I.");

            Console.WriteLine("");
            Console.WriteLine("  Maximum degree P = " + p + "");
            Console.WriteLine("  Exponential argument coefficient B = " + b + "");

            table = Hermite.hn_exponential_product(p, b);

            typeMethods.r8mat_print(p + 1, p + 1, table, "  Exponential product table:");
        }

        static void hermite_polynomial_test09(int p, int e)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HERMITE_POLYNOMIAL_TEST09 tests HN_POWER_PRODUCT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 February 2012
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
            Console.WriteLine("HERMITE_POLYNOMIAL_TEST09");
            Console.WriteLine("  Compute a normalized physicist''s Hermite power product table.");
            Console.WriteLine("");
            Console.WriteLine("  Tij = integral ( -oo < X < +oo ) X^E Hn(I,X) Hn(J,X) exp(-X*X) dx");
            Console.WriteLine("");
            Console.WriteLine("  where Hn(I,X) = normalized physicist''s Hermite polynomial of degree I.");

            Console.WriteLine("");
            Console.WriteLine("  Maximum degree P = " + p + "");
            Console.WriteLine("  Exponent of X, E = " + e + "");

            table = Hermite.hn_power_product(p, e);

            typeMethods.r8mat_print(p + 1, p + 1, table, "  Power product table:");
        }

        static void hermite_polynomial_test10(int p, double b)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HERMITE_POLYNOMIAL_TEST10 tests HEN_EXPONENTIAL_PRODUCT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 February 2012
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
            Console.WriteLine("HERMITE_POLYNOMIAL_TEST10");
            Console.WriteLine("  Compute a normalized probabilist''s Hermite exponential product table.");
            Console.WriteLine("");
            Console.WriteLine("  Tij = integral ( -oo < X < +oo ) exp(B*X) Hen(I,X) Hen(J,X) exp(-0.5*X*X) dx");
            Console.WriteLine("");
            Console.WriteLine("  where Hen(I,X) = normalized probabilist''s Hermite polynomial of degree I.");

            Console.WriteLine("");
            Console.WriteLine("  Maximum degree P = " + p + "");
            Console.WriteLine("  Exponential argument coefficient B = " + b + "");

            table = Hermite.hen_exponential_product(p, b);

            typeMethods.r8mat_print(p + 1, p + 1, table, "  Exponential product table:");
        }

        static void hermite_polynomial_test11(int p, int e)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HERMITE_POLYNOMIAL_TEST11 tests HEN_POWER_PRODUCT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 February 2012
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
            Console.WriteLine("HERMITE_POLYNOMIAL_TEST11");
            Console.WriteLine("  Compute a normalized probabilist''s Hermite power product table.");
            Console.WriteLine("");
            Console.WriteLine("  Tij = integral ( -oo < X < +oo ) X^E Hen(I,X) Hen(J,X) exp(-X*X) dx");
            Console.WriteLine("");
            Console.WriteLine("  where Hen(I,X) = normalized probabilist''s Hermite polynomial of degree I.");

            Console.WriteLine("");
            Console.WriteLine("  Maximum degree P = " + p + "");
            Console.WriteLine("  Exponent of X, E = " + e + "");

            table = Hermite.hen_power_product(p, e);

            typeMethods.r8mat_print(p + 1, p + 1, table, "  Power product table:");
        }

        static void hermite_polynomial_test12(int p, double b)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HERMITE_POLYNOMIAL_TEST12 tests HF_EXPONENTIAL_PRODUCT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 February 2012
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
            Console.WriteLine("HERMITE_POLYNOMIAL_TEST12");
            Console.WriteLine("  Compute a Hermite function exponential product table.");
            Console.WriteLine("");
            Console.WriteLine("  Tij = integral ( -oo < X < +oo ) exp(B*X) Hf(I,X) Hf(J,X) dx");
            Console.WriteLine("");
            Console.WriteLine("  where Hf(I,X) = Hermite function of \"degree\" I.");

            Console.WriteLine("");
            Console.WriteLine("  Maximum degree P = " + p + "");
            Console.WriteLine("  Exponential argument coefficient B = " + b + "");

            table = Hermite.hf_exponential_product(p, b);

            typeMethods.r8mat_print(p + 1, p + 1, table, "  Exponential product table:");
        }

        static void hermite_polynomial_test13(int p, int e)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HERMITE_POLYNOMIAL_TEST13 tests HF_POWER_PRODUCT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 February 2012
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
            Console.WriteLine("HERMITE_POLYNOMIAL_TEST13");
            Console.WriteLine("  Compute a Hermite function product table.");
            Console.WriteLine("");
            Console.WriteLine("  Tij = integral ( -oo < X < +oo ) X^E Hf(I,X) Hf(J,X) exp(-X*X) dx");
            Console.WriteLine("");
            Console.WriteLine("  where Hf(I,X) = Hermite function of \"degree\" I.");

            Console.WriteLine("");
            Console.WriteLine("  Maximum degree P = " + p + "");
            Console.WriteLine("  Exponent of X, E = " + e + "");

            table = Hermite.hf_power_product(p, e);

            typeMethods.r8mat_print(p + 1, p + 1, table, "  Power product table:");
        }

        static void hermite_polynomial_test14()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HERMITE_POLYNOMIAL_TEST14 tests H_POLYNOMIAL_COEFFICIENTS.
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
            double[] c;
            int i;
            int j;
            int n = 10;

            Console.WriteLine("");
            Console.WriteLine("HERMITE_POLYNOMIAL_TEST14");
            Console.WriteLine("  H_POLYNOMIAL_COEFFICIENTS determines physicist's Hermite polynomial coefficients.");

            c = Hermite.h_polynomial_coefficients(n);

            for (i = 0; i <= n; i++)
            {
                Console.WriteLine("");
                Console.WriteLine("  H(" + i + ",x) =");
                Console.WriteLine("");
                for (j = i; 0 <= j; j--)
                {
                    if (c[i + j * (n + 1)] == 0.0)
                    {
                    }
                    else if (j == 0)
                    {
                        Console.WriteLine(c[i + j * (n + 1)].ToString().PadLeft(14) + "");
                        ;
                    }
                    else if (j == 1)
                    {
                        Console.WriteLine(c[i + j * (n + 1)].ToString().PadLeft(14) + " * x");
                    }
                    else
                    {
                        Console.WriteLine(c[i + j * (n + 1)].ToString().PadLeft(14) + " * x^" + j + "");
                    }
                }
            }
        }

        static void hermite_polynomial_test15()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HERMITE_POLYNOMIAL_TEST15 tests HE_POLYNOMIAL_COEFFICIENTS.
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
            double[] c;
            int i;
            int j;
            int n = 10;

            Console.WriteLine("");
            Console.WriteLine("HERMITE_POLYNOMIAL_TEST15");
            Console.WriteLine("  HE_POLYNOMIAL_COEFFICIENTS determines probabilist's Hermite polynomial coefficients.");

            c = Hermite.he_polynomial_coefficients(n);

            for (i = 0; i <= n; i++)
            {
                Console.WriteLine("");
                Console.WriteLine("  He(" + i + ") =");
                Console.WriteLine("");
                for (j = i; 0 <= j; j--)
                {
                    if (c[i + j * (n + 1)] == 0.0)
                    {
                    }
                    else if (j == 0)
                    {
                        Console.WriteLine(c[i + j * (n + 1)].ToString().PadLeft(14) + "");
                        ;
                    }
                    else if (j == 1)
                    {
                        Console.WriteLine(c[i + j * (n + 1)].ToString().PadLeft(14) + " * x");
                    }
                    else
                    {
                        Console.WriteLine(c[i + j * (n + 1)].ToString().PadLeft(14) + " * x^" + j + "");
                    }
                }
            }
        }
    }
}