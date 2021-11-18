using System;
using Burkardt.ChebyshevPolynomialNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace ChebyshevPolynomialTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for CHEBYSHEV_POLYNOMIAL_TEST.
        //
        //  Discussion:
        //
        //    CHEBYSHEV_POLYNOMIAL_TEST tests the CHEBYSHEV_POLYNOMIAL library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("CHEBYSHEV_POLYNOMIAL_TEST");
        Console.WriteLine("  Test the CHEBYSHEV_POLYNOMIAL library.");

        test01();
        t_mass_matrix_test();
        t_moment_test();
        t_polynomial_test();
        t_polynomial_ab_test();
        t_polynomial_ab_value_test();
        t_polynomial_coefficients_test();
        t_polynomial_plot_test();
        t_polynomial_value_test();
        t_polynomial_zeros_test();
        t_quadrature_rule_test();
        test07();
        test08();
        test09();
        test10();
        tt_product_test();
        tt_product_integral_test();
        ttt_product_integral_test();
        tu_product_test();

        u_mass_matrix_test();
        u_moment_test();
        u_polynomial_test();
        u_polynomial_ab_test();
        u_polynomial_ab_value_test();
        u_polynomial_coefficients_test();
        u_polynomial_plot_test();
        u_polynomial_value_test();
        u_polynomial_zeros_test();
        u_quadrature_rule_test();
        uu_product_test();
        uu_product_integral_test();

        v_mass_matrix_test();
        v_moment_test();
        v_polynomial_test();
        v_polynomial_ab_test();
        v_polynomial_ab_value_test();
        v_polynomial_coefficients_test();
        v_polynomial_plot_test();
        v_polynomial_value_test();
        v_polynomial_zeros_test();
        v_quadrature_rule_test();
        vv_product_integral_test();

        w_mass_matrix_test();
        w_moment_test();
        w_polynomial_test();
        w_polynomial_ab_test();
        w_polynomial_ab_value_test();
        w_polynomial_coefficients_test();
        w_polynomial_plot_test();
        w_polynomial_value_test();
        w_polynomial_zeros_test();
        w_quadrature_rule_test();
        ww_product_integral_test();

        Console.WriteLine("");
        Console.WriteLine("CHEBYSHEV_POLYNOMIAL_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests T_PROJECT_COEFFICIENTS_DATA.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 June 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a;
        double b;
        double[] c;
        double[] d;
        double[] d2;
        int i;
        int m;
        int n;
        int seed;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("CHEBYSHEV_POLYNOMIAL_TEST01:");
        Console.WriteLine("  T_PROJECT_COEFFICIENTS_DATA estimates the Chebyshev polynomial");
        Console.WriteLine("  coefficients for a function given as data (x,fx).");
        Console.WriteLine("");
        Console.WriteLine("  Here, we use fx = f(x) = x^2 for the data.");
        Console.WriteLine("");
        Console.WriteLine("  Since T(0,x) = 1 and T(2,x) = 2*x^2 - 1, the correct expansion is");
        Console.WriteLine("  f(x) = 1/2 T(0,x) + 0 T(1,x) + 1/2 T(2,x) + 0 * all other polys,");
        Console.WriteLine("  if Chebyshev polynomials are based in [-1,+1].");
        //
        //  Data in [0,1];
        //
        a = 0.0;
        b = 1.0;

        Console.WriteLine("");
        Console.WriteLine("  Chebyshev polynomial will be based in [" + a + "," + b + "]");
        //
        //  Compute sample data.
        //
        m = 20;
        seed = 123456789;
        x = UniformRNG.r8vec_uniform_ab_new(m, a, b, ref seed);
        d = new double[m];
        for (i = 0; i < m; i++)
        {
            d[i] = x[i] * x[i];
        }

        typeMethods.r8vec2_print(m, x, d, "  Data ( X, D ):");

        n = 4;
        c = ChebyshevPolynomial.t_project_coefficients_data(a, b, m, n, x, d);

        typeMethods.r8vec_print(n, c, "  Coefficients of Chebyshev expansion of degree 4.");
        //
        //  Compare Chebyshev expansion and original function.
        //
        d2 = ChebyshevPolynomial.t_project_value_ab(m, n, x, c, a, b);

        Console.WriteLine("");
        Console.WriteLine("   I      X(I)     Data(I)      Chebyshev(X(I))");
        Console.WriteLine("");
        for (i = 0; i < m; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + d[i].ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + d2[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    private static void t_mass_matrix_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    T_MASS_MATRIX_TEST tests T_MASS_MATRIX.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a;
        int n;

        Console.WriteLine("");
        Console.WriteLine("T_MASS_MATRIX_TEST:");
        Console.WriteLine("  T_MASS_MATRIX computes the mass matrix for the");
        Console.WriteLine("  Chebyshev polynomials T(i,x).");
        Console.WriteLine("  A(I,J) = integral ( -1 <=x <= +1 ) T(i,x) T(j,x) / sqrt ( 1 - x^2 ) dx");
        Console.WriteLine("  0    if i is not equal to j;");
        Console.WriteLine("  pi   if i = j = 0;");
        Console.WriteLine("  pi/2 if i = j =/= 0.");

        n = 3;
        a = ChebyshevPolynomial.t_mass_matrix(n);

        typeMethods.r8mat_print(n + 1, n + 1, a, "  T mass matrix:");
    }

    private static void t_moment_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    T_MOMENT_TEST tests T_MOMENT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("T_MOMENT_TEST:");
        Console.WriteLine("  T_MOMENT returns the value of");
        Console.WriteLine("  integral ( -1 <=x <= +1 ) x^e / sqrt ( 1 - x^2 ) dx");
        Console.WriteLine("");
        Console.WriteLine("   E       Integral");
        Console.WriteLine("");
        for (int e = 0; e <= 10; e++)
        {
            double value = ChebyshevPolynomial.t_moment(e);
            Console.WriteLine("  " + e.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + value.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void t_polynomial_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    T_POLYNOMIAL_TEST tests T_POLYNOMIAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        int n = 0;
        double x = 0;
        double[] x_vec = new double[1];

        Console.WriteLine("");
        Console.WriteLine("T_POLYNOMIAL_TEST:");
        Console.WriteLine("  T_POLYNOMIAL evaluates the Chebyshev polynomial T(n,x).");
        Console.WriteLine("");
        Console.WriteLine("                   Tabulated      Computed");
        Console.WriteLine("     N      X        T(n,x)        T(n,x)");
        Console.WriteLine("");

        int n_data = 0;

        for (;;)
        {
            ChebyshevPolynomial.t_polynomial_values(ref n_data, ref n, ref x, ref fx);

            if (n_data == 0)
            {
                break;
            }

            switch (n)
            {
                case < 0:
                    continue;
            }

            x_vec[0] = x;
            double[] fx2 = ChebyshevPolynomial.t_polynomial(1, n, x_vec);

            Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + fx.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + fx2[n].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

        }
    }

    private static void t_polynomial_ab_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    T_POLYNOMIAL_AB_TEST tests T_POLYNOMIAL_AB.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a;
        double b;
        int m;
        int n;
        double[] v;
        double[] x;

        m = 11;
        n = 5;

        Console.WriteLine("");
        Console.WriteLine("T_POLYNOMIAL_AB_TEST:");
        Console.WriteLine("  T_POLYNOMIAL_AB evaluates Chebyshev polynomials TAB(n,x)");
        Console.WriteLine("  shifted from [-1,+1] to the domain [A,B].");
        Console.WriteLine("");
        Console.WriteLine("  Here, we will use the new domain [0,1]");
        Console.WriteLine("  and the desired maximum polynomial degree will be N = 5.");

        a = 0.0;
        b = 1.0;
        x = typeMethods.r8vec_linspace_new(m, a, b);

        v = ChebyshevPolynomial.t_polynomial_ab(a, b, m, n, ref x);

        typeMethods.r8mat_print(m, n + 1, v, "  Tables of T values:");

    }

    private static void t_polynomial_ab_value_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    T_POLYNOMIAL_AB_VALUE_TEST tests T_POLYNOMIAL_AB_VALUE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        int n = 0;
        double x01 = 0;

        Console.WriteLine("");
        Console.WriteLine("T_POLYNOMIAL_AB_VALUE_TEST:");
        Console.WriteLine("  T_POLYNOMIAL_AB_VALUE evaluates Chebyshev polynomials TAB(n,x)");
        Console.WriteLine("  shifted from [-1,+1] to the domain [A,B].");
        Console.WriteLine("");
        Console.WriteLine("  Here, we will use the new domain [0,1].");
        Console.WriteLine("");
        Console.WriteLine("                   Tabulated      Computed");
        Console.WriteLine("     N      X01    T01(n,x)       T01(n,x)");
        Console.WriteLine("");

        double a = 0.0;
        double b = 1.0;

        int n_data = 0;

        for (;;)
        {
            ChebyshevPolynomial.t_polynomial_01_values(ref n_data, ref n, ref x01, ref fx);

            if (n_data == 0)
            {
                break;
            }

            double fx2 = ChebyshevPolynomial.t_polynomial_ab_value(a, b, n, x01);

            Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + x01.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + fx.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + fx2.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }
    //****************************************************************************80

    private static void t_polynomial_coefficients_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    T_POLYNOMIAL_COEFFICIENTS_TEST tests T_POLYNOMIAL_COEFFICIENTS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;

        Console.WriteLine("");
        Console.WriteLine("T_POLYNOMIAL_COEFFICIENTS_TEST");
        Console.WriteLine("  T_POLYNOMIAL_COEFFICIENTS determines the polynomial coefficients ");
        Console.WriteLine("  of T(n,x).");

        int n = 5;

        double[] c = ChebyshevPolynomial.t_polynomial_coefficients(n);

        for (i = 0; i <= n; i++)
        {
            double[] c2 = new double[i + 1];
            for (j = 0; j <= i; j++)
            {
                c2[j] = c[i + j * (n + 1)];
            }

            typeMethods.r8poly_print(i, c2, "");
        }
    }

    private static void t_polynomial_plot_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    T_POLYNOMIAL_PLOT_TEST tests T_POLYNOMIAL_PLOT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n_num = 6;
        int[] n_val = new int[6];

        Console.WriteLine("");
        Console.WriteLine("T_POLYNOMIAL_PLOT_TEST");
        Console.WriteLine("  T_POLYNOMIAL_PLOT plots selected");
        Console.WriteLine("  Chebyshev polynomials T(n,x).");

        for (int i = 0; i <= 5; i++)
        {
            n_val[i] = i;
        }

        string output_filename = "t_polynomial_plot.png";

        ChebyshevPolynomial.t_polynomial_plot(n_num, n_val, output_filename);

    }

    private static void t_polynomial_value_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    T_POLYNOMIAL_VALUE_TEST tests T_POLYNOMIAL_VALUE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        int n = 0;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("T_POLYNOMIAL_VALUE_TEST:");
        Console.WriteLine("  T_POLYNOMIAL_VALUE evaluates the Chebyshev polynomial T(n,x).");
        Console.WriteLine("");
        Console.WriteLine("                   Tabulated      Computed");
        Console.WriteLine("     N      X        T(n,x)        T(n,x)");
        Console.WriteLine("");

        int n_data = 0;

        for (;;)
        {
            ChebyshevPolynomial.t_polynomial_values(ref n_data, ref n, ref x, ref fx);

            if (n_data == 0)
            {
                break;
            }

            double fx2 = ChebyshevPolynomial.t_polynomial_value(n, x);

            Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + fx.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + fx2.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void t_polynomial_zeros_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    T_POLYNOMIAL_ZEROS_TEST tests T_POLYNOMIAL_ZEROS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] fx;
        int i;
        int n;
        const int N_MAX = 5;
        double[] z;

        Console.WriteLine("");
        Console.WriteLine("T_POLYNOMIAL_ZEROS_TEST:");
        Console.WriteLine("  T_POLYNOMIAL_ZEROS returns zeroes of T(n,x).");
        Console.WriteLine("");
        Console.WriteLine("       N      X        T(n,x)");
        Console.WriteLine("");

        for (n = 1; n <= n_max; n++)
        {
            z = ChebyshevPolynomial.t_polynomial_zeros(n);
            fx = ChebyshevPolynomial.t_polynomial(n, n, z);
            for (i = 0; i < n; i++)
            {
                Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + z[i].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + fx[i + n * n].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }

            Console.WriteLine("");
        }
    }

    private static void t_quadrature_rule_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    T_QUADRATURE_RULE_TEST tests T_QUADRATURE_RULE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("T_QUADRATURE_RULE_TEST:");
        Console.WriteLine("  T_QUADRATURE_RULE computes the quadrature rule");
        Console.WriteLine("  associated with T(n,x);");

        int n = 7;
        double[] x = new double[n];
        double[] w = new double[n];

        ChebyshevPolynomial.t_quadrature_rule(n, ref x, ref w);

        typeMethods.r8vec2_print(n, x, w, "    N      X            W");

        Console.WriteLine("");
        Console.WriteLine("  Use the quadrature rule to estimate:");
        Console.WriteLine("");
        Console.WriteLine("    Q = Integral ( -1 <= X <= +1 ) X^E / sqrt ( 1-x^2) dx");
        Console.WriteLine("");
        Console.WriteLine("   E       Q_Estimate      Q_Exact");
        Console.WriteLine("");

        double[] f = new double[n];

        for (int e = 0; e <= 2 * n - 1; e++)
        {
            switch (e)
            {
                case 0:
                {
                    for (int i = 0; i < n; i++)
                    {
                        f[i] = 1.0;
                    }

                    break;
                }
                default:
                {
                    for (int i = 0; i < n; i++)
                    {
                        f[i] = Math.Pow(x[i], e);
                    }

                    break;
                }
            }

            double q = typeMethods.r8vec_dot_product(n, w, f);
            double q_exact = ChebyshevPolynomial.t_moment(e);
            Console.WriteLine("  " + e.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + q.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + q_exact.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void test07()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST07 tests T_PROJECT_COEFFICIENTS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a;
        double b;
        double[] c;
        int n;

        Console.WriteLine("");
        Console.WriteLine("TEST07:");
        Console.WriteLine("  T_PROJECT_COEFFICIENTS computes the Chebyshev coefficients");
        Console.WriteLine("  of a function defined over [-1,+1].");
        Console.WriteLine("  T_PROJECT_COEFFICIENTS_AB works in [A,B].");

        n = 3;
        c = new double[n + 1];
        c = ChebyshevPolynomial.t_project_coefficients(n, Math.Exp);
        typeMethods.r8vec_print(n + 1, c, "  Chebyshev coefficients for exp(x) in [-1,+1]");

        n = 5;
        c = new double[n + 1];
        c = ChebyshevPolynomial.t_project_coefficients(n, Math.Exp);
        typeMethods.r8vec_print(n + 1, c, "  Chebyshev coefficients for exp(x) in [-1,+1]");

        n = 5;
        c = new double[n + 1];
        c = ChebyshevPolynomial.t_project_coefficients(n, Math.Sin);
        typeMethods.r8vec_print(n + 1, c, "  Chebyshev coefficients for sin(x) in [-1,+1]");
        //
        //  Repeat calculation with T_PROJECT_COEFFICIENTS_AB.
        //
        n = 5;
        c = new double[n + 1];
        a = -1.0;
        b = +1.0;
        c = ChebyshevPolynomial.t_project_coefficients_ab(n, Math.Sin, a, b);
        typeMethods.r8vec_print(n + 1, c, "  Chebyshev coefficients for sin(x) in [-1,+1]");
        //
        //  Now try a different interval.
        //
        n = 5;
        c = new double[n + 1];
        a = 0.0;
        b = 1.0;
        c = ChebyshevPolynomial.t_project_coefficients_ab(n, Math.Sqrt, a, b);
        typeMethods.r8vec_print(n + 1, c, "  Chebyshev coefficients for sqrt(x) in [0,+1]");
    }

    private static void test08()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST08 tests T_PROJECT_COEFFICIENTS_DATA.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a;
        double b;
        double[] c;
        double[] d;
        int i;
        int m;
        int n;
        int seed;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST08:");
        Console.WriteLine("  T_PROJECT_COEFFICIENTS_DATA computes the Chebyshev");
        Console.WriteLine("  coefficients of a function defined by data.");
        Console.WriteLine("");
        Console.WriteLine("  We are looking for an approximation that is good in [-1,+1].");
        Console.WriteLine("");
        Console.WriteLine("  Begin by using equally spaced points in [-1,+1].");

        a = -1.0;
        b = +1.0;
        m = 10;
        x = typeMethods.r8vec_linspace_new(m, a, b);
        d = new double[m];
        for (i = 0; i < m; i++)
        {
            d[i] = Math.Exp(x[i]);
        }

        n = 3;
        c = ChebyshevPolynomial.t_project_coefficients_data(a, b, m, n, x, d);
        typeMethods.r8vec_print(n + 1, c, "  Chebyshev coefficients for exp(x) on [-1,+1]");

        a = -1.0;
        b = +1.0;
        m = 10;
        x = typeMethods.r8vec_linspace_new(m, a, b);
        d = new double[m];
        for (i = 0; i < m; i++)
        {
            d[i] = Math.Exp(x[i]);
        }

        n = 5;
        c = ChebyshevPolynomial.t_project_coefficients_data(a, b, m, n, x, d);
        typeMethods.r8vec_print(n + 1, c, "  Chebyshev coefficients for exp(x) on [-1,+1]");

        a = -1.0;
        b = +1.0;
        m = 10;
        x = typeMethods.r8vec_linspace_new(m, a, b);
        d = new double[m];
        for (i = 0; i < m; i++)
        {
            d[i] = Math.Sin(x[i]);
        }

        n = 5;
        c = ChebyshevPolynomial.t_project_coefficients_data(a, b, m, n, x, d);
        typeMethods.r8vec_print(n + 1, c, "  Chebyshev coefficients for sin(x) on [-1,+1]");

        Console.WriteLine("");
        Console.WriteLine("  Now sample equally spaced points in [0,+1].");
        Console.WriteLine("  The approximation still applies to the interval [-1,+1].");

        a = 0.0;
        b = +1.0;
        m = 10;
        x = typeMethods.r8vec_linspace_new(m, a, b);
        d = new double[m];
        for (i = 0; i < m; i++)
        {
            d[i] = Math.Sin(x[i]);
        }

        n = 5;
        c = ChebyshevPolynomial.t_project_coefficients_data(a, b, m, n, x, d);
        typeMethods.r8vec_print(n + 1, c, "  Chebyshev coefficients for sin(x) on [0,+1]");

        a = 0.0;
        b = +1.0;
        m = 10;
        x = typeMethods.r8vec_linspace_new(m, a, b);
        d = new double[m];
        for (i = 0; i < m; i++)
        {
            d[i] = Math.Sqrt(x[i]);
        }

        n = 5;
        c = ChebyshevPolynomial.t_project_coefficients_data(a, b, m, n, x, d);
        typeMethods.r8vec_print(n + 1, c, "  Chebyshev coefficients for sqrt(x) on [0,+1]");

        Console.WriteLine("");
        Console.WriteLine("  Now random points in [-1,+1].");

        a = -1.0;
        b = +1.0;
        m = 10;
        seed = 123456789;
        x = UniformRNG.r8vec_uniform_ab_new(m, a, b, ref seed);
        d = new double[m];
        for (i = 0; i < m; i++)
        {
            d[i] = Math.Sin(x[i]);
        }

        n = 5;
        c = ChebyshevPolynomial.t_project_coefficients_data(a, b, m, n, x, d);
        typeMethods.r8vec_print(n + 1, c, "  Chebyshev coefficients for sin(x) on [-1,+1]");
    }

    private static void test09()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST09 compares a function and projection over [-1,+1].
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a;
        double b;
        double[] c;
        int i;
        int m;
        int n;
        double r;
        double[] v;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST09:");
        Console.WriteLine("  T_PROJECT_COEFFICIENTS computes the Chebyshev interpolant C(F)(n,x)");
        Console.WriteLine("  of a function F(x) defined over [-1,+1].");
        Console.WriteLine("  T_PROJECT_VALUE evaluates that projection.");

        Console.WriteLine("");
        Console.WriteLine("  Compute projections of order N to exp(x) over [-1,+1],");
        Console.WriteLine("");
        Console.WriteLine("   N   Max||F(x)-C(F)(n,x)||");
        Console.WriteLine("");

        a = -1.0;
        b = +1.0;

        for (n = 0; n <= 10; n++)
        {
            c = ChebyshevPolynomial.t_project_coefficients(n, Math.Exp);
            m = 101;
            x = typeMethods.r8vec_linspace_new(m, a, b);
            v = ChebyshevPolynomial.t_project_value(m, n, x, c);
            r = 0.0;
            for (i = 0; i < m; i++)
            {
                r = Math.Max(r, Math.Abs(v[i] - Math.Exp(x[i])));
            }

            Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + r.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void test10()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST10 compares a function and projection over [A,B].
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a;
        double b;
        double[] c;
        int i;
        int m;
        int n;
        double r;
        double[] v;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST10:");
        Console.WriteLine("  T_PROJECT_COEFFICIENTS_AB computes the Chebyshev interpolant C(F)(n,x)");
        Console.WriteLine("  of a function F(x) defined over [A,B].");
        Console.WriteLine("  T_PROJECT_VALUE_AB evaluates that projection.");

        a = 0.0;
        b = 1.5;

        Console.WriteLine("");
        Console.WriteLine("  Compute projections of order N to exp(x) over [" + a + "," + b + "]");
        Console.WriteLine("");
        Console.WriteLine("   N   Max||F(x)-C(F)(n,x)||");
        Console.WriteLine("");

        for (n = 0; n <= 10; n++)
        {
            c = ChebyshevPolynomial.t_project_coefficients_ab(n, Math.Exp, a, b);
            m = 101;
            x = typeMethods.r8vec_linspace_new(m, a, b);
            v = ChebyshevPolynomial.t_project_value_ab(m, n, x, c, a, b);
            r = 0.0;
            for (i = 0; i < m; i++)
            {
                r = Math.Max(r, Math.Abs(v[i] - Math.Exp(x[i])));
            }

            Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + r.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void tt_product_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TT_PRODUCT_TEST tests TT_PRODUCT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;
        double r8_hi;
        double r8_lo;
        int seed;
        int test;
        double ti;
        double titj;
        double tj;
        double x;

        Console.WriteLine("");
        Console.WriteLine("TT_PRODUCT_TEST:");
        Console.WriteLine("  TT_PRODUCT(I,J;X) = T(I,X) * T(J,X)");

        r8_lo = -1.0;
        r8_hi = +1.0;
        seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("   I   J      X               TI              TJ              TI*TJ       TT_PRODUCT");
        Console.WriteLine("");
        for (test = 1; test <= 10; test++)
        {
            x = UniformRNG.r8_uniform_ab(r8_lo, r8_hi, ref seed);
            i = UniformRNG.i4_uniform_ab(0, 6, ref seed);
            ti = ChebyshevPolynomial.t_polynomial_value(i, x);
            j = UniformRNG.i4_uniform_ab(-1, 4, ref seed);
            tj = ChebyshevPolynomial.t_polynomial_value(j, x);
            titj = ChebyshevPolynomial.tt_product(i, j, x);
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + j.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + ti.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + tj.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + (ti * tj).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + titj.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void tt_product_integral_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TT_PRODUCT_INTEGRAL_TEST tests TT_PRODUCT_INTEGRAL.
        //
        //  Discussion:
        //
        //    This process should match the T_MASS_MATRIX computation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a;
        int i;
        int j;
        int n;

        Console.WriteLine("");
        Console.WriteLine("TT_PRODUCT_INTEGRAL_TEST:");
        Console.WriteLine("  TT_PRODUCT_INTEGRAL computes the product integral");
        Console.WriteLine("  of a pair of Chebyshev T polynomials T(i,x) and T(j,x).");
        Console.WriteLine("  A(I,J) = integral ( -1 <=x <= +1 ) T(i,x) T(j,x) / sqrt ( 1 - x^2 ) dx");
        Console.WriteLine("  0    if i is not equal to j;");
        Console.WriteLine("  pi   if i = j = 0;");
        Console.WriteLine("  pi/2 if i = j =/= 0.");

        n = 4;
        a = new double[(n + 1) * (n + 1)];
        for (i = 0; i <= n; i++)
        {
            for (j = 0; j <= n; j++)
            {
                a[i + j * (n + 1)] = ChebyshevPolynomial.tt_product_integral(i, j);
            }
        }

        typeMethods.r8mat_print(n + 1, n + 1, a, "  T(i,x)*T(j,x) integral matrix:");
    }

    private static void ttt_product_integral_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TTT_PRODUCT_INTEGRAL_TEST tests TTT_PRODUCT_INTEGRAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1;
        double fx2;
        int i;
        int j;
        int k;
        int l;
        int n;
        int seed;
        int test;
        int test_num = 20;
        double ti;
        double tj;
        double tk;
        double[] w;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TTT_PRODUCT_INTEGRAL_TEST:");
        Console.WriteLine("  TTT_PRODUCT_INTEGRAL computes the triple integral");
        Console.WriteLine("  Tijk = integral ( -1 <= x <= 1 ) T(i,x) T(j,x) T(k,x) / sqrt ( 1-x^2) dx");
        Console.WriteLine("");
        Console.WriteLine("   I   J   K     Tijk           Tijk");
        Console.WriteLine("                 computed       exact");
        Console.WriteLine("");

        n = 15;
        x = new double[n];
        w = new double[n];

        ChebyshevPolynomial.t_quadrature_rule(n, ref x, ref w);

        seed = 123456789;

        for (test = 1; test <= test_num; test++)
        {
            i = UniformRNG.i4_uniform_ab(2, 6, ref seed);
            j = UniformRNG.i4_uniform_ab(1, 3, ref seed);
            k = UniformRNG.i4_uniform_ab(0, 4, ref seed);
            fx1 = ChebyshevPolynomial.ttt_product_integral(i, j, k);
            fx2 = 0.0;
            for (l = 0; l < n; l++)
            {
                ti = ChebyshevPolynomial.t_polynomial_value(i, x[l]);
                tj = ChebyshevPolynomial.t_polynomial_value(j, x[l]);
                tk = ChebyshevPolynomial.t_polynomial_value(k, x[l]);
                fx2 += w[l] * ti * tj * tk;
            }

            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + j.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + k.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + fx1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + fx2.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void tu_product_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TU_PRODUCT_TEST tests TU_PRODUCT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;
        double r8_hi;
        double r8_lo;
        int seed;
        int test;
        double ti;
        double tiuj;
        double uj;
        double x;

        Console.WriteLine("");
        Console.WriteLine("TU_PRODUCT_TEST:");
        Console.WriteLine("  TU_PRODUCT(I,J;X) = T(I,X) * U(J,X)");

        r8_lo = -1.0;
        r8_hi = +1.0;
        seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("   I   J      X               TI              UJ              TI*UJ       TU_PRODUCT");
        Console.WriteLine("");
        for (test = 1; test <= 10; test++)
        {
            x = UniformRNG.r8_uniform_ab(r8_lo, r8_hi, ref seed);
            i = UniformRNG.i4_uniform_ab(0, 6, ref seed);
            ti = ChebyshevPolynomial.t_polynomial_value(i, x);
            j = UniformRNG.i4_uniform_ab(-1, 4, ref seed);
            uj = ChebyshevPolynomial.u_polynomial_value(j, x);
            tiuj = ChebyshevPolynomial.tu_product(i, j, x);
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + j.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + ti.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + uj.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + (ti * uj).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + tiuj.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

    }

    private static void u_mass_matrix_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    U_MASS_MATRIX_TEST tests U_MASS_MATRIX.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a;
        int n;

        Console.WriteLine("");
        Console.WriteLine("U_MASS_MATRIX_TEST:");
        Console.WriteLine("  U_MASS_MATRIX computes the mass matrix for the");
        Console.WriteLine("  Chebyshev U polynomials U(i,x).");
        Console.WriteLine("  A(I,J) = integral ( -1 <=x <= +1 ) U(i,x) U(j,x) * sqrt ( 1 - x^2 ) dx");
        Console.WriteLine("  0    if i is not equal to j;");
        Console.WriteLine("  pi/2 if i = j.");

        n = 3;
        a = ChebyshevPolynomial.u_mass_matrix(n);

        typeMethods.r8mat_print(n + 1, n + 1, a, "  U mass matrix:");
    }

    private static void u_moment_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    U_MOMENT_TEST tests U_MOMENT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int e;
        double value = 0;

        Console.WriteLine("");
        Console.WriteLine("U_MOMENT_TEST:");
        Console.WriteLine("  U_MOMENT returns the value of");
        Console.WriteLine("  integral ( -1 <=x <= +1 ) x^e * sqrt ( 1 - x^2 ) dx");
        Console.WriteLine("");
        Console.WriteLine("   E       Integral");
        Console.WriteLine("");
        for (e = 0; e <= 10; e++)
        {
            value = ChebyshevPolynomial.u_moment(e);
            Console.WriteLine("  " + e.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + value.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void u_polynomial_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    U_POLYNOMIAL_TEST tests U_POLYNOMIAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        double[] fx2;
        int n = 0;
        double x = 0;
        double[] x_vec = new double[1];

        Console.WriteLine("");
        Console.WriteLine("U_POLYNOMIAL_TEST:");
        Console.WriteLine("  U_POLYNOMIAL evaluates the Chebyshev polynomial U(n,x).");
        Console.WriteLine("");
        Console.WriteLine("                   Tabulated      Computed");
        Console.WriteLine("     N      X        U(n,x)        U(n,x)");
        Console.WriteLine("");

        int n_data = 0;

        for (;;)
        {
            ChebyshevPolynomial.u_polynomial_values(ref n_data, ref n, ref x, ref fx);

            if (n_data == 0)
            {
                break;
            }

            switch (n)
            {
                case < 0:
                    continue;
            }

            x_vec[0] = x;
            fx2 = ChebyshevPolynomial.u_polynomial(1, n, x_vec);

            Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + fx.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + fx2[n].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void u_polynomial_ab_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    U_POLYNOMIAL_AB_TEST tests U_POLYNOMIAL_AB.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a;
        double b;
        int m;
        int n;
        double[] v;
        double[] x;

        m = 11;
        n = 5;

        Console.WriteLine("");
        Console.WriteLine("U_POLYNOMIAL_AB_TEST:");
        Console.WriteLine("  U_POLYNOMIAL_AB evaluates Chebyshev polynomials UAB(n,x)");
        Console.WriteLine("  shifted from [-1,+1] to the domain [A,B].");
        Console.WriteLine("");
        Console.WriteLine("  Here, we will use the new domain [0,1]");
        Console.WriteLine("  and the desired maximum polynomial degree will be N = 5.");

        a = 0.0;
        b = 1.0;
        x = typeMethods.r8vec_linspace_new(m, a, b);

        v = ChebyshevPolynomial.u_polynomial_ab(a, b, m, n, x);

        typeMethods.r8mat_print(m, n + 1, v, "  Tables of U values:");

    }

    private static void u_polynomial_ab_value_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    U_POLYNOMIAL_AB_VALUE_TEST tests U_POLYNOMIAL_AB_VALUE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        int n = 0;
        double x01 = 0;

        Console.WriteLine("");
        Console.WriteLine("U_POLYNOMIAL_AB_VALUE_TEST:");
        Console.WriteLine("  U_POLYNOMIAL_AB_VALUE evaluates Chebyshev polynomials UAB(n,x)");
        Console.WriteLine("  shifted from [-1,+1] to the domain [A,B].");
        Console.WriteLine("");
        Console.WriteLine("  Here, we will use the new domain [0,1].");
        Console.WriteLine("");
        Console.WriteLine("                   Tabulated      Computed");
        Console.WriteLine("     N      X01    U01(n,x)       U01(n,x)");
        Console.WriteLine("");

        double a = 0.0;
        double b = 1.0;

        int n_data = 0;

        for (;;)
        {
            ChebyshevPolynomial.u_polynomial_01_values(ref n_data, ref n, ref x01, ref fx);

            if (n_data == 0)
            {
                break;
            }

            double fx2 = ChebyshevPolynomial.u_polynomial_ab_value(a, b, n, x01);

            Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + x01.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + fx.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + fx2.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void u_polynomial_coefficients_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    U_POLYNOMIAL_COEFFICIENTS_TEST tests U_POLYNOMIAL_COEFFICIENTS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] c;
        double[] c2;
        int i;
        int j;
        int n;

        Console.WriteLine("");
        Console.WriteLine("U_POLYNOMIAL_COEFFICIENTS_TEST");
        Console.WriteLine("  U_POLYNOMIAL_COEFFICIENTS determines the polynomial coefficients ");
        Console.WriteLine("  of U(n,x).");

        n = 5;

        c = ChebyshevPolynomial.u_polynomial_coefficients(n);

        for (i = 0; i <= n; i++)
        {
            c2 = new double[i + 1];
            for (j = 0; j <= i; j++)
            {
                c2[j] = c[i + j * (n + 1)];
            }

            typeMethods.r8poly_print(i, c2, "");
        }
    }

    private static void u_polynomial_plot_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    U_POLYNOMIAL_PLOT_TEST tests U_POLYNOMIAL_PLOT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int n_num = 6;
        int[] n_val = new int[6];
        string output_filename;

        Console.WriteLine("");
        Console.WriteLine("U_POLYNOMIAL_PLOT_TEST");
        Console.WriteLine("  U_POLYNOMIAL_PLOT plots selected");
        Console.WriteLine("  Chebyshev polynomials U(n,x).");

        for (i = 0; i <= 5; i++)
        {
            n_val[i] = i;
        }

        output_filename = "u_polynomial_plot.png";

        ChebyshevPolynomial.u_polynomial_plot(n_num, n_val, output_filename);
    }

    private static void u_polynomial_value_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    U_POLYNOMIAL_VALUE_TEST tests U_POLYNOMIAL_VALUE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        int n = 0;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("U_POLYNOMIAL_VALUE_TEST:");
        Console.WriteLine("  U_POLYNOMIAL_VALUE evaluates the Chebyshev polynomial U(n,x).");
        Console.WriteLine("");
        Console.WriteLine("                   Tabulated      Computed");
        Console.WriteLine("     N      X        U(n,x)        U(n,x)");
        Console.WriteLine("");

        int n_data = 0;

        for (;;)
        {
            ChebyshevPolynomial.u_polynomial_values(ref n_data, ref n, ref x, ref fx);

            if (n_data == 0)
            {
                break;
            }

            double fx2 = ChebyshevPolynomial.u_polynomial_value(n, x);

            Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + fx.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + fx2.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void u_polynomial_zeros_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    U_POLYNOMIAL_ZEROS_TEST tests U_POLYNOMIAL_ZEROS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] fx;
        int i;
        int n;
        const int N_MAX = 5;
        double[] z;

        Console.WriteLine("");
        Console.WriteLine("U_POLYNOMIAL_ZEROS_TEST:");
        Console.WriteLine("  U_POLYNOMIAL_ZEROS returns zeroes of U(n,x).");
        Console.WriteLine("");
        Console.WriteLine("       N      X        U(n,x)");
        Console.WriteLine("");

        for (n = 1; n <= n_max; n++)
        {
            z = ChebyshevPolynomial.u_polynomial_zeros(n);
            fx = ChebyshevPolynomial.u_polynomial(n, n, z);
            for (i = 0; i < n; i++)
            {
                Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + z[i].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + fx[i + n * n].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }

            Console.WriteLine("");
        }
    }

    private static void u_quadrature_rule_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    U_QUADRATURE_RULE_TEST tests U_QUADRATURE_RULE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 April 2012
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
        Console.WriteLine("U_QUADRATURE_RULE_TEST:");
        Console.WriteLine("  U_QUADRATURE_RULE computes the quadrature rule");
        Console.WriteLine("  associated with U(n,x);");

        n = 7;
        x = new double[n];
        w = new double[n];

        ChebyshevPolynomial.u_quadrature_rule(n, ref x, ref w);

        typeMethods.r8vec2_print(n, x, w, "    N      X            W");

        Console.WriteLine("");
        Console.WriteLine("  Use the quadrature rule to estimate:");
        Console.WriteLine("");
        Console.WriteLine("    Q = Integral ( -1 <= X <= +1 ) X^E * sqrt ( 1-x^2) dx");
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
            q_exact = ChebyshevPolynomial.u_moment(e);
            Console.WriteLine("  " + e.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + q.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + q_exact.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void uu_product_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UU_PRODUCT_TEST tests UU_PRODUCT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;
        double r8_hi;
        double r8_lo;
        int seed;
        int test;
        double ui;
        double uiuj;
        double uj;
        double x;

        Console.WriteLine("");
        Console.WriteLine("UU_PRODUCT_TEST:");
        Console.WriteLine("  UU_PRODUCT(I,J;X) = U(I,X) * U(J,X)");

        r8_lo = -1.0;
        r8_hi = +1.0;
        seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("   I   J      X               UI              UJ              UI*UJ       UU_PRODUCT");
        Console.WriteLine("");
        for (test = 1; test <= 10; test++)
        {
            x = UniformRNG.r8_uniform_ab(r8_lo, r8_hi, ref seed);
            i = UniformRNG.i4_uniform_ab(0, 6, ref seed);
            ui = ChebyshevPolynomial.u_polynomial_value(i, x);
            j = UniformRNG.i4_uniform_ab(-1, 4, ref seed);
            uj = ChebyshevPolynomial.u_polynomial_value(j, x);
            uiuj = ChebyshevPolynomial.uu_product(i, j, x);
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + j.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + ui.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + uj.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + (ui * uj).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + uiuj.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void uu_product_integral_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UU_PRODUCT_INTEGRAL_TEST tests UU_PRODUCT_INTEGRAL.
        //
        //  Discussion:
        //
        //    This process should match the U_MASS_MATRIX computation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a;
        int i;
        int j;
        int n;

        Console.WriteLine("");
        Console.WriteLine("UU_PRODUCT_INTEGRAL_TEST:");
        Console.WriteLine("  UU_PRODUCT_INTEGRAL computes the product integral");
        Console.WriteLine("  of a pair of Chebyshev U polynomials U(i,x) and U(j,x).");
        Console.WriteLine("  A(I,J) = integral ( -1 <=x <= +1 ) U(i,x) U(j,x) * sqrt ( 1 - x^2 ) dx");
        Console.WriteLine("  0    if i is not equal to j;");
        Console.WriteLine("  pi/2 if i = j.");

        n = 4;
        a = new double[(n + 1) * (n + 1)];
        for (i = 0; i <= n; i++)
        {
            for (j = 0; j <= n; j++)
            {
                a[i + j * (n + 1)] = ChebyshevPolynomial.uu_product_integral(i, j);
            }
        }

        typeMethods.r8mat_print(n + 1, n + 1, a, "  U(i,x)*U(j,x) integral matrix:");

    }

    private static void v_mass_matrix_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    V_MASS_MATRIX_TEST tests V_MASS_MATRIX.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a;
        int n;

        Console.WriteLine("");
        Console.WriteLine("V_MASS_MATRIX_TEST:");
        Console.WriteLine("  V_MASS_MATRIX computes the mass matrix for the");
        Console.WriteLine("  Chebyshev polynomials V(i,x).");
        Console.WriteLine("  A(I,J) = integral ( -1 <=x <= +1 ) V(i,x) V(j,x) sqrt(1+x)/sqrt(1-x) dx");
        Console.WriteLine("  0  if i is not equal to j;");
        Console.WriteLine("  pi if i = j.");

        n = 3;
        a = ChebyshevPolynomial.v_mass_matrix(n);

        typeMethods.r8mat_print(n + 1, n + 1, a, "  V mass matrix:");

    }

    private static void v_moment_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    V_MOMENT_TEST tests V_MOMENT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int e;
        double value = 0;

        Console.WriteLine("");
        Console.WriteLine("V_MOMENT_TEST:");
        Console.WriteLine("  V_MOMENT returns the value of");
        Console.WriteLine("  integral ( -1 <=x <= +1 ) x^e * sqrt(1+x)/sqrt(1-x) dx");
        Console.WriteLine("");
        Console.WriteLine("   E       Integral");
        Console.WriteLine("");
        for (e = 0; e <= 10; e++)
        {
            value = ChebyshevPolynomial.v_moment(e);
            Console.WriteLine("  " + e.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + value.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void v_polynomial_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    V_POLYNOMIAL_TEST tests V_POLYNOMIAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        int n = 0;
        double x = 0;
        double[] x_vec = new double[1];

        Console.WriteLine("");
        Console.WriteLine("V_POLYNOMIAL_TEST:");
        Console.WriteLine("  V_POLYNOMIAL evaluates the Chebyshev polynomial V(n,x).");
        Console.WriteLine("");
        Console.WriteLine("                   Tabulated      Computed");
        Console.WriteLine("     N      X        V(n,x)        V(n,x)");
        Console.WriteLine("");

        int n_data = 0;

        for (;;)
        {
            ChebyshevPolynomial.v_polynomial_values(ref n_data, ref n, ref x, ref fx);

            if (n_data == 0)
            {
                break;
            }

            switch (n)
            {
                case < 0:
                    continue;
            }

            x_vec[0] = x;
            double[] fx2 = ChebyshevPolynomial.v_polynomial(1, n, x_vec);

            Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + fx.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + fx2[n].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

        }
    }

    private static void v_polynomial_ab_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    V_POLYNOMIAL_AB_TEST tests V_POLYNOMIAL_AB.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a;
        double b;
        int m;
        int n;
        double[] v;
        double[] x;

        m = 11;
        n = 5;

        Console.WriteLine("");
        Console.WriteLine("V_POLYNOMIAL_AB_TEST:");
        Console.WriteLine("  V_POLYNOMIAL_AB evaluates Chebyshev polynomials VAB(n,x)");
        Console.WriteLine("  shifted from [-1,+1] to the domain [A,B].");
        Console.WriteLine("");
        Console.WriteLine("  Here, we will use the new domain [0,1]");
        Console.WriteLine("  and the desired maximum polynomial degree will be N = 5.");

        a = 0.0;
        b = 1.0;
        x = typeMethods.r8vec_linspace_new(m, a, b);

        v = ChebyshevPolynomial.v_polynomial_ab(a, b, m, n, x);

        typeMethods.r8mat_print(m, n + 1, v, "  Tables of V values:");
    }

    private static void v_polynomial_ab_value_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    V_POLYNOMIAL_AB_VALUE_TEST tests V_POLYNOMIAL_AB_VALUE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        int n = 0;
        double x01 = 0;

        Console.WriteLine("");
        Console.WriteLine("V_POLYNOMIAL_AB_VALUE_TEST:");
        Console.WriteLine("  V_POLYNOMIAL_AB_VALUE evaluates Chebyshev polynomials VAB(n,x)");
        Console.WriteLine("  shifted from [-1,+1] to the domain [A,B].");
        Console.WriteLine("");
        Console.WriteLine("  Here, we will use the new domain [0,1].");
        Console.WriteLine("");
        Console.WriteLine("                   Tabulated      Computed");
        Console.WriteLine("     N      X01    V01(n,x)       V01(n,x)");
        Console.WriteLine("");

        double a = 0.0;
        double b = 1.0;

        int n_data = 0;

        for (;;)
        {
            ChebyshevPolynomial.v_polynomial_01_values(ref n_data, ref n, ref x01, ref fx);

            if (n_data == 0)
            {
                break;
            }

            double fx2 = ChebyshevPolynomial.v_polynomial_ab_value(a, b, n, x01);

            Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + x01.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + fx.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + fx2.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void v_polynomial_coefficients_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    V_POLYNOMIAL_COEFFICIENTS_TEST tests V_POLYNOMIAL_COEFFICIENTS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("V_POLYNOMIAL_COEFFICIENTS_TEST");
        Console.WriteLine("  V_POLYNOMIAL_COEFFICIENTS determines the polynomial coefficients ");
        Console.WriteLine("  of V(n,x).");

        int n = 5;

        double[] c = ChebyshevPolynomial.v_polynomial_coefficients(n);

        for (int i = 0; i <= n; i++)
        {
            double[] c2 = new double[i + 1];
            for (int j = 0; j <= i; j++)
            {
                c2[j] = c[i + j * (n + 1)];
            }

            typeMethods.r8poly_print(i, c2, "");
        }
    }

    private static void v_polynomial_plot_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    V_POLYNOMIAL_PLOT_TEST tests V_POLYNOMIAL_PLOT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n_num = 6;
        int[] n_val = new int[6];
        string output_filename;

        Console.WriteLine("");
        Console.WriteLine("V_POLYNOMIAL_PLOT_TEST");
        Console.WriteLine("  V_POLYNOMIAL_PLOT plots selected");
        Console.WriteLine("  Chebyshev polynomials V(n,x).");

        for (int i = 0; i <= 5; i++)
        {
            n_val[i] = i;
        }

        output_filename = "v_polynomial_plot.png";

        ChebyshevPolynomial.v_polynomial_plot(n_num, n_val, output_filename);
    }

    private static void v_polynomial_value_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    V_POLYNOMIAL_VALUE_TEST tests V_POLYNOMIAL_VALUE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        int n = 0;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("V_POLYNOMIAL_VALUE_TEST:");
        Console.WriteLine("  V_POLYNOMIAL_VALUE evaluates the Chebyshev polynomial V(n,x).");
        Console.WriteLine("");
        Console.WriteLine("                   Tabulated      Computed");
        Console.WriteLine("     N      X        V(n,x)        V(n,x)");
        Console.WriteLine("");

        int n_data = 0;

        for (;;)
        {
            ChebyshevPolynomial.v_polynomial_values(ref n_data, ref n, ref x, ref fx);

            if (n_data == 0)
            {
                break;
            }

            double fx2 = ChebyshevPolynomial.v_polynomial_value(n, x);

            Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + fx.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + fx2.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void v_polynomial_zeros_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    V_POLYNOMIAL_ZEROS_TEST tests V_POLYNOMIAL_ZEROS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N_MAX = 5;

        Console.WriteLine("");
        Console.WriteLine("V_POLYNOMIAL_ZEROS_TEST:");
        Console.WriteLine("  V_POLYNOMIAL_ZEROS returns zeroes of V(n,x).");
        Console.WriteLine("");
        Console.WriteLine("       N      X        V(n,x)");
        Console.WriteLine("");

        for (int n = 1; n <= n_max; n++)
        {
            double[] z = ChebyshevPolynomial.v_polynomial_zeros(n);
            double[] fx = ChebyshevPolynomial.v_polynomial(n, n, z);
            for (int i = 0; i < n; i++)
            {
                Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + z[i].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + fx[i + n * n].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }

            Console.WriteLine("");
        }
    }

    private static void v_quadrature_rule_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    V_QUADRATURE_RULE_TEST tests V_QUADRATURE_RULE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("V_QUADRATURE_RULE_TEST:");
        Console.WriteLine("  V_QUADRATURE_RULE computes the quadrature rule");
        Console.WriteLine("  associated with V(n,x);");

        int n = 7;
        double[] x = new double[n];
        double[] w = new double[n];

        ChebyshevPolynomial.v_quadrature_rule(n, ref x, ref w);

        typeMethods.r8vec2_print(n, x, w, "    N      X            W");

        Console.WriteLine("");
        Console.WriteLine("  Use the quadrature rule to estimate:");
        Console.WriteLine("");
        Console.WriteLine("    Q = Integral ( -1 <= X <= +1 ) X^E * sqrt(1+x)/sqrt(1-x) dx");
        Console.WriteLine("");
        Console.WriteLine("   E       Q_Estimate      Q_Exact");
        Console.WriteLine("");

        double[] f = new double[n];

        for (int e = 0; e <= 2 * n - 1; e++)
        {
            switch (e)
            {
                case 0:
                {
                    for (int i = 0; i < n; i++)
                    {
                        f[i] = 1.0;
                    }

                    break;
                }
                default:
                {
                    for (int i = 0; i < n; i++)
                    {
                        f[i] = Math.Pow(x[i], e);
                    }

                    break;
                }
            }

            double q = typeMethods.r8vec_dot_product(n, w, f);
            double q_exact = ChebyshevPolynomial.v_moment(e);
            Console.WriteLine("  " + e.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + q.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + q_exact.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void vv_product_integral_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VV_PRODUCT_INTEGRAL_TEST tests VV_PRODUCT_INTEGRAL.
        //
        //  Discussion:
        //
        //    This process should match the V_MASS_MATRIX computation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a;
        int i;
        int j;
        int n;

        Console.WriteLine("");
        Console.WriteLine("VV_PRODUCT_INTEGRAL_TEST:");
        Console.WriteLine("  VV_PRODUCT_INTEGRAL computes the product integral");
        Console.WriteLine("  of a pair of Chebyshev V polynomials V(i,x) and V(j,x).");
        Console.WriteLine("  A(I,J) = integral ( -1 <=x <= +1 ) V(i,x) V(j,x) * sqrt ( 1 + x ) / sqrt ( 1 - x ) dx");
        Console.WriteLine("  0  if i is not equal to j;");
        Console.WriteLine("  pi if i = j.");

        n = 4;
        a = new double[(n + 1) * (n + 1)];
        for (i = 0; i <= n; i++)
        {
            for (j = 0; j <= n; j++)
            {
                a[i + j * (n + 1)] = ChebyshevPolynomial.vv_product_integral(i, j);
            }
        }

        typeMethods.r8mat_print(n + 1, n + 1, a, "  V(i,x)*V(j,x) integral matrix:");
    }

    private static void w_mass_matrix_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    W_MASS_MATRIX_TEST tests W_MASS_MATRIX.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a;
        int n;

        Console.WriteLine("");
        Console.WriteLine("W_MASS_MATRIX_TEST:");
        Console.WriteLine("  W_MASS_MATRIX computes the mass matrix for the");
        Console.WriteLine("  Chebyshev polynomials W(i,x).");
        Console.WriteLine("  A(I,J) = integral ( -1 <=x <= +1 ) W(i,x) W(j,x) sqrt(1-x)/sqrt(1+x) dx");
        Console.WriteLine("  0  if i is not equal to j;");
        Console.WriteLine("  pi if i = j.");

        n = 3;
        a = ChebyshevPolynomial.w_mass_matrix(n);

        typeMethods.r8mat_print(n + 1, n + 1, a, "  W mass matrix:");

    }

    private static void w_moment_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    W_MOMENT_TEST tests W_MOMENT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int e;
        double value = 0;

        Console.WriteLine("");
        Console.WriteLine("W_MOMENT_TEST:");
        Console.WriteLine("  W_MOMENT returns the value of");
        Console.WriteLine("  integral ( -1 <=x <= +1 ) x^e * sqrt(1-x)/sqrt(1+x) dx");
        Console.WriteLine("");
        Console.WriteLine("   E       Integral");
        Console.WriteLine("");
        for (e = 0; e <= 10; e++)
        {
            value = ChebyshevPolynomial.w_moment(e);
            Console.WriteLine("  " + e.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + value.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void w_polynomial_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    W_POLYNOMIAL_TEST tests W_POLYNOMIAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        int n = 0;
        double x = 0;
        double[] x_vec = new double[1];

        Console.WriteLine("");
        Console.WriteLine("W_POLYNOMIAL_TEST:");
        Console.WriteLine("  W_POLYNOMIAL evaluates the Chebyshev polynomial W(n,x).");
        Console.WriteLine("");
        Console.WriteLine("                   Tabulated      Computed");
        Console.WriteLine("     N      X        W(n,x)        W(n,x)");
        Console.WriteLine("");

        int n_data = 0;

        for (;;)
        {
            ChebyshevPolynomial.w_polynomial_values(ref n_data, ref n, ref x, ref fx);

            if (n_data == 0)
            {
                break;
            }

            switch (n)
            {
                case < 0:
                    continue;
            }

            x_vec[0] = x;
            double[] fx2 = ChebyshevPolynomial.w_polynomial(1, n, x_vec);

            Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + fx.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + fx2[n].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void w_polynomial_ab_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    W_POLYNOMIAL_AB_TEST tests W_POLYNOMIAL_AB.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a;
        double b;
        int m;
        int n;
        double[] v;
        double[] x;

        m = 11;
        n = 5;

        Console.WriteLine("");
        Console.WriteLine("W_POLYNOMIAL_AB_TEST:");
        Console.WriteLine("  W_POLYNOMIAL_AB evaluates Chebyshev polynomials WAB(n,x)");
        Console.WriteLine("  shifted from [-1,+1] to the domain [A,B].");
        Console.WriteLine("");
        Console.WriteLine("  Here, we will use the new domain [0,1]");
        Console.WriteLine("  and the desired maximum polynomial degree will be N = 5.");

        a = 0.0;
        b = 1.0;
        x = typeMethods.r8vec_linspace_new(m, a, b);

        v = ChebyshevPolynomial.w_polynomial_ab(a, b, m, n, x);

        typeMethods.r8mat_print(m, n + 1, v, "  Tables of W values:");
    }

    private static void w_polynomial_ab_value_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    W_POLYNOMIAL_AB_VALUE_TEST tests W_POLYNOMIAL_AB_VALUE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        int n = 0;
        double x01 = 0;

        Console.WriteLine("");
        Console.WriteLine("W_POLYNOMIAL_AB_VALUE_TEST:");
        Console.WriteLine("  W_POLYNOMIAL_AB_VALUE evaluates Chebyshev polynomials WAB(n,x)");
        Console.WriteLine("  shifted from [-1,+1] to the domain [A,B].");
        Console.WriteLine("");
        Console.WriteLine("  Here, we will use the new domain [0,1].");
        Console.WriteLine("");
        Console.WriteLine("                   Tabulated      Computed");
        Console.WriteLine("     N      X01    W01(n,x)       W01(n,x)");
        Console.WriteLine("");

        double a = 0.0;
        double b = 1.0;

        int n_data = 0;

        for (;;)
        {
            ChebyshevPolynomial.w_polynomial_01_values(ref n_data, ref n, ref x01, ref fx);

            if (n_data == 0)
            {
                break;
            }

            double fx2 = ChebyshevPolynomial.w_polynomial_ab_value(a, b, n, x01);

            Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + x01.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + fx.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + fx2.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void w_polynomial_coefficients_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    W_POLYNOMIAL_COEFFICIENTS_TEST tests W_POLYNOMIAL_COEFFICIENTS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("W_POLYNOMIAL_COEFFICIENTS_TEST");
        Console.WriteLine("  W_POLYNOMIAL_COEFFICIENTS determines the polynomial coefficients ");
        Console.WriteLine("  of W(n,x).");

        int n = 5;

        double[] c = ChebyshevPolynomial.w_polynomial_coefficients(n);

        for (int i = 0; i <= n; i++)
        {
            double[] c2 = new double[i + 1];
            for (int j = 0; j <= i; j++)
            {
                c2[j] = c[i + j * (n + 1)];
            }

            typeMethods.r8poly_print(i, c2, "");
        }
    }

    private static void w_polynomial_plot_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    W_POLYNOMIAL_PLOT_TEST tests W_POLYNOMIAL_PLOT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int n_num = 6;
        int[] n_val = new int[6];
        string output_filename;

        Console.WriteLine("");
        Console.WriteLine("W_POLYNOMIAL_PLOT_TEST");
        Console.WriteLine("  W_POLYNOMIAL_PLOT plots selected");
        Console.WriteLine("  Chebyshev polynomials W(n,x).");

        for (i = 0; i <= 5; i++)
        {
            n_val[i] = i;
        }

        output_filename = "w_polynomial_plot.png";

        ChebyshevPolynomial.w_polynomial_plot(n_num, n_val, output_filename);
    }

    private static void w_polynomial_value_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    W_POLYNOMIAL_VALUE_TEST tests W_POLYNOMIAL_VALUE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        int n = 0;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("W_POLYNOMIAL_VALUE_TEST:");
        Console.WriteLine("  W_POLYNOMIAL_VALUE evaluates the Chebyshev polynomial W(n,x).");
        Console.WriteLine("");
        Console.WriteLine("                   Tabulated      Computed");
        Console.WriteLine("     N      X        W(n,x)        W(n,x)");
        Console.WriteLine("");

        int n_data = 0;

        for (;;)
        {
            ChebyshevPolynomial.w_polynomial_values(ref n_data, ref n, ref x, ref fx);

            if (n_data == 0)
            {
                break;
            }

            double fx2 = ChebyshevPolynomial.w_polynomial_value(n, x);

            Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + fx.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + fx2.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void w_polynomial_zeros_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    W_POLYNOMIAL_ZEROS_TEST tests W_POLYNOMIAL_ZEROS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int n;
        const int N_MAX = 5;

        Console.WriteLine("");
        Console.WriteLine("W_POLYNOMIAL_ZEROS_TEST:");
        Console.WriteLine("  W_POLYNOMIAL_ZEROS returns zeroes of W(n,x).");
        Console.WriteLine("");
        Console.WriteLine("       N      X        W(n,x)");
        Console.WriteLine("");

        for (n = 1; n <= n_max; n++)
        {
            double[] z = ChebyshevPolynomial.w_polynomial_zeros(n);
            double[] fx = ChebyshevPolynomial.w_polynomial(n, n, z);
            for (i = 0; i < n; i++)
            {
                Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + z[i].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + fx[i + n * n].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }

            Console.WriteLine("");
        }
    }

    private static void w_quadrature_rule_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    W_QUADRATURE_RULE_TEST tests W_QUADRATURE_RULE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("W_QUADRATURE_RULE_TEST:");
        Console.WriteLine("  W_QUADRATURE_RULE computes the quadrature rule");
        Console.WriteLine("  associated with W(n,x);");

        int n = 7;
        double[] x = new double[n];
        double[] w = new double[n];

        ChebyshevPolynomial.w_quadrature_rule(n, ref x, ref w);

        typeMethods.r8vec2_print(n, x, w, "    N      X            W");

        Console.WriteLine("");
        Console.WriteLine("  Use the quadrature rule to estimate:");
        Console.WriteLine("");
        Console.WriteLine("    Q = Integral ( -1 <= X <= +1 ) X^E * sqrt(1-x)/sqrt(1+x) dx");
        Console.WriteLine("");
        Console.WriteLine("   E       Q_Estimate      Q_Exact");
        Console.WriteLine("");

        double[] f = new double[n];

        for (int e = 0; e <= 2 * n - 1; e++)
        {
            switch (e)
            {
                case 0:
                {
                    for (int i = 0; i < n; i++)
                    {
                        f[i] = 1.0;
                    }

                    break;
                }
                default:
                {
                    for (int i = 0; i < n; i++)
                    {
                        f[i] = Math.Pow(x[i], e);
                    }

                    break;
                }
            }

            double q = typeMethods.r8vec_dot_product(n, w, f);
            double q_exact = ChebyshevPolynomial.w_moment(e);
            Console.WriteLine("  " + e.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + q.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + q_exact.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void ww_product_integral_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WW_PRODUCT_INTEGRAL_TEST tests WW_PRODUCT_INTEGRAL.
        //
        //  Discussion:
        //
        //    This process should match the W_MASS_MATRIX computation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("WW_PRODUCT_INTEGRAL_TEST:");
        Console.WriteLine("  WW_PRODUCT_INTEGRAL computes the product integral");
        Console.WriteLine("  of a pair of Chebyshev W polynomials W(i,x) and W(j,x).");
        Console.WriteLine("  A(I,J) = integral ( -1 <=x <= +1 ) W(i,x) W(j,x) * sqrt ( 1 - x ) / sqrt ( 1 + x ) dx");
        Console.WriteLine("  0  if i is not equal to j;");
        Console.WriteLine("  pi if i = j.");

        int n = 4;
        double[] a = new double[(n + 1) * (n + 1)];
        for (int i = 0; i <= n; i++)
        {
            for (int j = 0; j <= n; j++)
            {
                a[i + j * (n + 1)] = ChebyshevPolynomial.ww_product_integral(i, j);
            }
        }

        typeMethods.r8mat_print(n + 1, n + 1, a, "  W(i,x)*W(j,x) integral matrix:");
    }
}