using System;
using System.Globalization;
using Burkardt.IntegralNS;
using Burkardt.Types;

namespace JacobiPolynomialTest;

using Polynomial = Burkardt.PolynomialNS.Jacobi;
using Quadrature = Burkardt.Quadrature.JacobiQuadrature;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for JACOBI_POLYNOMIAL_TEST.
        //
        //  Discussion:
        //
        //    JACOBI_POLYNOMIAL_TEST tests the JACOBI_POLYNOMIAL library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("JACOBI_POLYNOMIAL_TEST");
        Console.WriteLine("  Test the JACOBI_POLYNOMIAL library.");

        test01();
        test02();
        test03();
        test04();
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("JACOBI_POLYNOMIAL_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests J_POLYNOMIAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a = 0;
        double b = 0;
        double fx1 = 0;
        int n = 0;
        double x = 0;
        double[] x_vec = new double[1];

        Console.WriteLine("");
        Console.WriteLine("TEST01:");
        Console.WriteLine("  J_POLYNOMIAL_VALUES stores values of");
        Console.WriteLine("  the Jacobi polynomials.");
        Console.WriteLine("  J_POLYNOMIAL evaluates the polynomial.");
        Console.WriteLine("");
        Console.WriteLine("                                    Tabulated                 Computed" +
                          "     N     A     B        X           J(N,A,B,X)                    J(N,A,B,X)                     Error");
        Console.WriteLine("");

        int n_data = 0;

        for (;;)
        {
            Burkardt.Values.Jacobi.jacobi_poly_values(ref n_data, ref n, ref a, ref b, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            int m = 1;
            x_vec[0] = x;
            double[] fx2_vec = Polynomial.j_polynomial(m, n, a, b, x_vec);
            double fx2 = fx2_vec[0 + n * 1];
            double e = fx1 - fx2;

            Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + a.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + b.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + fx1.ToString(CultureInfo.InvariantCulture).PadLeft(24)
                                   + "  " + fx2.ToString(CultureInfo.InvariantCulture).PadLeft(24)
                                   + "  " + e.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");

        }
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests J_POLYNOMIAL_ZEROS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a_test =  {
                0.5, 1.0, 2.0
            }
            ;
        double[] b_test =  {
                0.5, 1.5, 0.5
            }
            ;
        int test;
        const int test_num = 3;

        Console.WriteLine("");
        Console.WriteLine("TEST02:");
        Console.WriteLine("  J_POLYNOMIAL_ZEROS computes the zeros of J(n,a,b,x);");
        Console.WriteLine("  Check by calling J_POLYNOMIAL there.");

        for (test = 0; test < test_num; test++)
        {
            double a = a_test[test];
            double b = b_test[test];

            int degree;
            for (degree = 1; degree <= 5; degree++)
            {
                double[] z = Polynomial.j_polynomial_zeros(degree, a, b);
                string title = "Zeros for J(" + degree.ToString(CultureInfo.InvariantCulture) + ","
                               + a.ToString(CultureInfo.InvariantCulture) + "," + b.ToString(CultureInfo.InvariantCulture) + ")";
                typeMethods.r8vec_print(degree, z, title);

                double[] hz = Polynomial.j_polynomial(degree, degree, a, b, z);
                title = "Evaluate J(" + degree.ToString(CultureInfo.InvariantCulture) + ","
                        + a.ToString(CultureInfo.InvariantCulture) + "," + b.ToString(CultureInfo.InvariantCulture) + ")";
                typeMethods.r8vec_print(degree, hz, title, aIndex:  + degree * degree);
            }
        }
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests J_QUADRATURE_RULE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("TEST03:");
        Console.WriteLine("  J_QUADRATURE_RULE computes the quadrature rule");
        Console.WriteLine("  associated with J(n,a,b,x);");

        const int n = 7;
        const double a = 1.0;
        const double b = 2.5;

        double[] x = new double[n];
        double[] w = new double[n];

        Quadrature.j_quadrature_rule(n, a, b, ref x, ref w);

        typeMethods.r8vec2_print(n, x, w, "      X            W");

        Console.WriteLine("");
        Console.WriteLine("  Use the quadrature rule to estimate:");
        Console.WriteLine("");
        Console.WriteLine("    Q = Integral (-1<x<+1) J(i,a,b,x) J(j,a,b,x) (1-x)^a (1+x)^b dx");
        Console.WriteLine("");
        Console.WriteLine("   I   J      Q_Estimate         Q_Exact");
        Console.WriteLine("");

        for (i = 0; i <= 5; i++)
        {
            double[] ji = Polynomial.j_polynomial(n, i, a, b, x);
            int j;
            for (j = i; j <= 5; j++)
            {
                double[] jj = Polynomial.j_polynomial(n, j, a, b, x);
                double q = 0.0;
                int k;
                for (k = 0; k < n; k++)
                {
                    q += w[k] * ji[k + i * n] * jj[k + j * n];
                }

                double q_exact = Integral.j_double_product_integral(i, j, a, b);
                Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                       + "  " + j.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                       + "  " + q.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + q_exact.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }
        }
    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests J_DOUBLE_PRODUCT_INTEGRAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("TEST04:");
        Console.WriteLine("  J_DOUBLE_PRODUCT_INTEGRAL returns the weighted integral of");
        Console.WriteLine("  J(i,a,b,x) * J(j,a,b,x);");

        const double a = 1.0;
        const double b = 2.5;

        Console.WriteLine("");
        Console.WriteLine("    Q = Integral (-1<x<+1) J(i,a,b,x) J(j,a,b,x) (1-x)^a (1+x)^b dx");
        Console.WriteLine("");
        Console.WriteLine("   I   J      Q");
        Console.WriteLine("");

        for (i = 0; i <= 5; i++)
        {
            int j;
            for (j = i; j <= 5; j++)
            {
                double q = Integral.j_double_product_integral(i, j, a, b);
                Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                       + "  " + j.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                       + "  " + q.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }
        }
    }
}