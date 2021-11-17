using System;
using Burkardt;
using Burkardt.PolynomialNS;
using Burkardt.Quadrature;
using Burkardt.Types;

namespace LegendrePolynomialTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for LEGENDRE_POLYNOMIAL_TEST.
        //
        //  Discussion:
        //
        //    LEGENDRE_POLYNOMIAL_TEST tests the LEGENDRE_POLYNOMIAL library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 March 2016
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
        Console.WriteLine("LEGENDRE_POLYNOMIAL_TEST:");
        Console.WriteLine("  Test the LEGENDRE_POLYNOMIAL library.");

        p = 5;
        b = 0.0;
        p_exponential_product_test(p, b);

        p = 5;
        b = 1.0;
        p_exponential_product_test(p, b);

        p_integral_test();

        p_polynomial_coefficients_test();
        p_polynomial_prime_test();
        p_polynomial_prime2_test();
        p_polynomial_value_test();
        p_polynomial_zeros_test();

        p = 5;
        e = 0;
        p_power_product_test(p, e);

        p = 5;
        e = 1;
        p_power_product_test(p, e);

        p_quadrature_rule_test();

        pm_polynomial_value_test();

        pmn_polynomial_value_test();

        pmns_polynomial_value_test();

        p = 5;
        pn_pair_product_test(p);
        pn_polynomial_coefficients_test();
        pn_polynomial_value_test();
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("LEGENDRE_POLYNOMIAL_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    public static void p_integral_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P_INTEGRAL_TEST tests P_INTEGRAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 March 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n;
        double value = 0;

        Console.WriteLine("");
        Console.WriteLine("P_INTEGRAL_TEST:");
        Console.WriteLine("  P_INTEGRAL returns the integral of P(n,x) over [-1,+1].");
        Console.WriteLine("");
        Console.WriteLine("     N        Integral");
        Console.WriteLine("");

        for (n = 0; n <= 10; n++)
        {
            value = Legendre.p_integral(n);

            Console.WriteLine("  " + n.ToString().PadLeft(4)
                                   + "  " + value.ToString().PadLeft(14) + "");
        }

    }

    public static void p_polynomial_value_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P_POLYNOMIAL_VALUE_TEST tests P_POLYNOMIAL_VALUE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 March 2012
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
        Console.WriteLine("P_POLYNOMIAL_VALUE_TEST:");
        Console.WriteLine("  P_POLYNOMIAL_VALUE evaluates the Legendre polynomial P(n,x).");
        Console.WriteLine("");
        Console.WriteLine("                        Tabulated                 Computed");
        Console.WriteLine("     N        X           P(N,X)                    P(N,X)                     Error");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            Burkardt.Values.Legendre.p_polynomial_values(ref n_data, ref n, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            x_vec[0] = x;
            fx2_vec = Legendre.p_polynomial_value(1, n, x_vec);
            fx2 = fx2_vec[n];

            e = fx1 - fx2;

            Console.WriteLine("  " + n.ToString().PadLeft(4)
                                   + "  " + x.ToString().PadLeft(12)
                                   + "  " + fx1.ToString("0.################").PadLeft(24)
                                   + "  " + fx2.ToString("0.################").PadLeft(24)
                                   + "  " + e.ToString().PadLeft(8) + "");
        }
    }

    public static void p_polynomial_prime_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P_POLYNOMIAL_PRIME_TEST tests P_POLYNOMIAL_PRIME.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 May 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a;
        double b;
        int i;
        int j;
        int m;
        int n;
        double[] vp;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("P_POLYNOMIAL_PRIME_TEST:");
        Console.WriteLine("  P_POLYNOMIAL_PRIME evaluates the derivative of the");
        Console.WriteLine("  Legendre polynomial P(n,x).");
        Console.WriteLine("");
        Console.WriteLine("                        Computed");
        Console.WriteLine("     N        X           P'(N,X)");
        Console.WriteLine("");

        m = 11;
        a = -1.0;
        b = +1.0;
        x = typeMethods.r8vec_linspace_new(m, a, b);

        n = 5;
        vp = Legendre.p_polynomial_prime(m, n, x);

        for (i = 0; i < m; i++)
        {
            Console.WriteLine("");
            for (j = 0; j <= n; j++)
            {
                Console.WriteLine("  " + j.ToString().PadLeft(4)
                                       + "  " + x[i].ToString().PadLeft(12)
                                       + "  " + vp[i + j * m].ToString("0.################").PadLeft(24) + "");
            }
        }
    }

    public static void p_polynomial_prime2_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P_POLYNOMIAL_PRIME2_TEST tests P_POLYNOMIAL_PRIME2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 May 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a;
        double b;
        int i;
        int j;
        int m;
        int n;
        double[] vpp;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("P_POLYNOMIAL_PRIME2_TEST:");
        Console.WriteLine("  P_POLYNOMIAL_PRIME2 evaluates the second derivative of the");
        Console.WriteLine("  Legendre polynomial P(n,x).");
        Console.WriteLine("");
        Console.WriteLine("                        Computed");
        Console.WriteLine("     N        X           P\"(N,X)");
        Console.WriteLine("");

        m = 11;
        a = -1.0;
        b = +1.0;
        x = typeMethods.r8vec_linspace_new(m, a, b);

        n = 5;
        vpp = Legendre.p_polynomial_prime2(m, n, x);

        for (i = 0; i < m; i++)
        {
            Console.WriteLine("");
            for (j = 0; j <= n; j++)
            {
                Console.WriteLine("  " + j.ToString().PadLeft(4)
                                       + "  " + x[i].ToString().PadLeft(12)
                                       + "  " + vpp[i + j * m].ToString("0.################").PadLeft(24) + "");
            }
        }
    }

    public static void p_polynomial_coefficients_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P_POLYNOMIAL_COEFFICIENTS_TEST tests P_POLYNOMIAL_COEFFICIENTS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 March 2012
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
        Console.WriteLine("P_POLYNOMIAL_COEFFICIENTS_TEST");
        Console.WriteLine("  P_POLYNOMIAL_COEFFICIENTS: coefficients of Legendre polynomial P(n,x).");

        c = Legendre.p_polynomial_coefficients(n);

        for (i = 0; i <= n; i++)
        {
            Console.WriteLine("");
            Console.WriteLine("  P(" + i + ",x) =");
            Console.WriteLine("");
            for (j = i; 0 <= j; j--)
            {
                switch (c[i + j * (n + 1)])
                {
                    case 0.0:
                        break;
                    default:
                        switch (j)
                        {
                            case 0:
                                Console.WriteLine(c[i + j * (n + 1)].ToString().PadLeft(14) + "");
                                ;
                                break;
                            case 1:
                                Console.WriteLine(c[i + j * (n + 1)].ToString().PadLeft(14) + " * x");
                                break;
                            default:
                                Console.WriteLine(c[i + j * (n + 1)].ToString().PadLeft(14) + " * x^" + j + "");
                                break;
                        }

                        break;
                }
            }
        }
    }

    public static void p_polynomial_zeros_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P_POLYNOMIAL_ZEROS_TEST tests P_POLYNOMIAL_ZEROS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 March 2012
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
        Console.WriteLine("P_POLYNOMIAL_ZEROS_TEST:");
        Console.WriteLine("  P_POLYNOMIAL_ZEROS computes the zeros of P(n,x)");
        Console.WriteLine("  Check by calling P_POLYNOMIAL_VALUE there.");

        for (degree = 1; degree <= 5; degree++)
        {
            z = Legendre.p_polynomial_zeros(degree);
            title = "  Computed zeros for P(" + degree + ",z):";
            typeMethods.r8vec_print(degree, z, title);

            hz = Legendre.p_polynomial_value(degree, degree, z);
            title = "  Evaluate P(" + degree + ",z):";
            typeMethods.r8vec_print(degree, hz, title, +degree * degree);
        }
    }

    public static void p_quadrature_rule_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P_QUADRATURE_RULE_TEST tests P_QUADRATURE_RULE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 March 2012
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
        Console.WriteLine("P_QUADRATURE_RULE_TEST:");
        Console.WriteLine("  P_QUADRATURE_RULE computes the quadrature rule");
        Console.WriteLine("  associated with P(n,x)");

        n = 7;
        x = new double[n];
        w = new double[n];

        LegendreQuadrature.p_quadrature_rule(n, ref x, ref w);

        typeMethods.r8vec2_print(n, x, w, "      X            W");

        Console.WriteLine("");
        Console.WriteLine("  Use the quadrature rule to estimate:");
        Console.WriteLine("");
        Console.WriteLine("    Q = Integral ( -1 <= X <= +1.0 ) X^E dx");
        Console.WriteLine("");
        Console.WriteLine("   E       Q_Estimate      Q_Exact");
        Console.WriteLine("");

        f = new double[n];

        for (e = 0; e <= 2 * n - 1; e++)
        {
            switch (e)
            {
                case 0:
                {
                    for (i = 0; i < n; i++)
                    {
                        f[i] = 1.0;
                    }

                    break;
                }
                default:
                {
                    for (i = 0; i < n; i++)
                    {
                        f[i] = Math.Pow(x[i], e);
                    }

                    break;
                }
            }

            q = typeMethods.r8vec_dot_product(n, w, f);
            q_exact = Legendre.p_integral(e);
            Console.WriteLine("  " + e.ToString().PadLeft(2)
                                   + "  " + q.ToString().PadLeft(14)
                                   + "  " + q_exact.ToString().PadLeft(14) + "");
        }
    }

    public static void p_exponential_product_test(int p, double b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P_EXPONENTIAL_PRODUCT_TEST tests P_EXPONENTIAL_PRODUCT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 March 2012
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
        Console.WriteLine("P_EXPONENTIAL_PRODUCT_TEST");
        Console.WriteLine("  P_EXPONENTIAL_PRODUCT_TEST computes a Legendre exponential product table.");
        Console.WriteLine("");
        Console.WriteLine("  Tij = integral ( -1.0 <= X <= +1.0 ) exp(B*X) P(I,X) P(J,X) dx");
        Console.WriteLine("");
        Console.WriteLine("  where P(I,X) = Legendre polynomial of degree I.");

        Console.WriteLine("");
        Console.WriteLine("  Maximum degree P = " + p + "");
        Console.WriteLine("  Exponential argument coefficient B = " + b + "");

        table = Legendre.p_exponential_product(p, b);

        typeMethods.r8mat_print(p + 1, p + 1, table, "  Exponential product table:");
    }

    public static void p_power_product_test(int p, int e)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P_POWER_PRODUCT_TEST tests P_POWER_PRODUCT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 March 2012
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
        Console.WriteLine("P_POWER_PRODUCT_TEST");
        Console.WriteLine("  P_POWER_PRODUCT_TEST computes a Legendre power product table.");
        Console.WriteLine("");
        Console.WriteLine("  Tij = integral ( -1.0 <= X <= +1.0 ) X^E P(I,X) P(J,X) dx");
        Console.WriteLine("");
        Console.WriteLine("  where P(I,X) = Legendre polynomial of degree I.");

        Console.WriteLine("");
        Console.WriteLine("  Maximum degree P = " + p + "");
        Console.WriteLine("  Exponent of X, E = " + e + "");

        table = Legendre.p_power_product(p, e);

        typeMethods.r8mat_print(p + 1, p + 1, table, "  Power product table:");
    }

    public static void pm_polynomial_value_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PM_POLYNOMIAL_VALUE_TEST tests PM_POLYNOMIAL_VALUE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 March 2012
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
        Console.WriteLine("PM_POLYNOMIAL_VALUE_TEST:");
        Console.WriteLine("  PM_POLYNOMIAL_VALUE evaluates the Legendre polynomial Pm(n,m,x).");
        Console.WriteLine("");
        Console.WriteLine("                             Tabulated                 Computed");
        Console.WriteLine(
            "     N     M        X        Pm(N,M,X)                 Pm(N,M,X)                     Error");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            Burkardt.Values.Legendre.pm_polynomial_values(ref n_data, ref n, ref m, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            x_vec[0] = x;
            fx2_vec = Legendre.pm_polynomial_value(1, n, m, x_vec);
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

    public static void pmn_polynomial_value_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PMN_POLYNOMIAL_VALUE_TEST tests PMN_POLYNOMIAL_VALUE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 March 2012
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
        Console.WriteLine("PMN_POLYNOMIAL_VALUE_TEST:");
        Console.WriteLine("  PMN_POLYNOMIAL_VALUE evaluates the Legendre polynomial Pmn(n,m,x).");
        Console.WriteLine("");
        Console.WriteLine("                             Tabulated                 Computed");
        Console.WriteLine(
            "     N     M        X       Pmn(N,M,X)                Pmn(N,M,X)                     Error");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            Burkardt.Values.Legendre.pmn_polynomial_values(ref n_data, ref n, ref m, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            x_vec[0] = x;
            fx2_vec = Legendre.pmn_polynomial_value(1, n, m, x_vec);
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

    public static void pmns_polynomial_value_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PMNS_POLYNOMIAL_VALUE_TEST tests PMNS_POLYNOMIAL_VALUE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 March 2012
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
        Console.WriteLine("PMNS_POLYNOMIAL_VALUE_TEST:");
        Console.WriteLine("  PMNS_POLYNOMIAL_VALUE evaluates the Legendre polynomial Pmns(n,m,x).");
        Console.WriteLine("");
        Console.WriteLine("                             Tabulated                 Computed");
        Console.WriteLine(
            "     N     M        X       Pmns(N,M,X)                Pmns(N,M,X)                     Error");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            Burkardt.Values.Legendre.pmns_polynomial_values(ref n_data, ref n, ref m, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            x_vec[0] = x;
            fx2_vec = Legendre.pmns_polynomial_value(1, n, m, x_vec);
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

    public static void pn_polynomial_coefficients_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PN_POLYNOMIAL_COEFFICIENTS_TEST tests PN_POLYNOMIAL_COEFFICIENTS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 October 2014
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
        Console.WriteLine("PN_POLYNOMIAL_COEFFICIENTS_TEST");
        Console.WriteLine("  PN_POLYNOMIAL_COEFFICIENTS: coefficients of normalized Legendre polynomial P(n,x).");

        c = Legendre.pn_polynomial_coefficients(n);

        for (i = 0; i <= n; i++)
        {
            Console.WriteLine("");
            Console.WriteLine("  P(" + i + ",x) =");
            Console.WriteLine("");
            for (j = i; 0 <= j; j--)
            {
                switch (c[i + j * (n + 1)])
                {
                    case 0.0:
                        break;
                    default:
                        switch (j)
                        {
                            case 0:
                                Console.WriteLine(c[i + j * (n + 1)].ToString().PadLeft(14) + "");
                                ;
                                break;
                            case 1:
                                Console.WriteLine(c[i + j * (n + 1)].ToString().PadLeft(14) + " * x");
                                break;
                            default:
                                Console.WriteLine(c[i + j * (n + 1)].ToString().PadLeft(14) + " * x^" + j + "");
                                break;
                        }

                        break;
                }
            }
        }
    }

    public static void pn_pair_product_test(int p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PN_PAIR_PRODUCT_TEST tests PN_PAIR_PRODUCT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 March 2012
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
    {
        double[] table;

        Console.WriteLine("");
        Console.WriteLine("PN_PAIR_PRODUCT_TEST:");
        Console.WriteLine("  PN_PAIR_PRODUCT_TEST computes a pair product table for Pn(n,x).");
        Console.WriteLine("");
        Console.WriteLine("  Tij = integral ( -1.0 <= X <= +1.0 ) Pn(I,X) Pn(J,X) dx");
        Console.WriteLine("");
        Console.WriteLine("  where Pn(I,X) = normalized Legendre polynomial of degree I.");
        Console.WriteLine("");
        Console.WriteLine("  Maximum degree P = " + p + "");

        table = Legendre.pn_pair_product(p);

        typeMethods.r8mat_print(p + 1, p + 1, table, "  Pair product table:");
    }

    public static void pn_polynomial_value_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PN_POLYNOMIAL_VALUE_TEST tests PN_POLYNOMIAL_VALUE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 March 2012
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
        Console.WriteLine("PN_POLYNOMIAL_VALUE_TEST:");
        Console.WriteLine("  PN_POLYNOMIAL_VALUE evaluates the normalized Legendre polynomial Pn(n,x).");
        Console.WriteLine("");
        Console.WriteLine("                        Tabulated                 Computed");
        Console.WriteLine("     N        X          Pn(N,X)                   Pn(N,X)                     Error");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            Burkardt.Values.Legendre.pn_polynomial_values(ref n_data, ref n, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            x_vec[0] = x;
            fx2_vec = Legendre.pn_polynomial_value(1, n, x_vec);
            fx2 = fx2_vec[n];

            e = fx1 - fx2;

            Console.WriteLine("  " + n.ToString().PadLeft(4)
                                   + "  " + x.ToString().PadLeft(12)
                                   + "  " + fx1.ToString("0.################").PadLeft(24)
                                   + "  " + fx2.ToString("0.################").PadLeft(24)
                                   + "  " + e.ToString().PadLeft(8) + "");
        }
    }
}