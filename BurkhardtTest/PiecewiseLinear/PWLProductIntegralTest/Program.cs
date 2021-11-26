using System;
using System.Globalization;
using Burkardt.PiecewiseLinear;

namespace PWLProductIntegralTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for pwl_PRODUCT_INTEGRAL_TEST.
        //
        //  Discussion:
        //
        //    pwl_PRODUCT_INTEGRAL_TEST tests pwl_PRODUCT_INTEGRAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 April 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("pwl_PRODUCT_INTEGRAL_TEST");
        Console.WriteLine("  Test the pwl_PRODUCT_INTEGRAL_INTEGRAL library.");

        test01();
        test02();
        test03();
        test04();

        Console.WriteLine("");
        Console.WriteLine("pwl_PRODUCT_INTEGRAL_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests pwl_PRODUCT_INTEGRAL.
        //
        //  Discussion:
        //
        //    For the first test, we use the same single "piece" for both F and G.
        //    Hence, we are actually integrating X^2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 April 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int F_NUM = 2;
        const int G_NUM = 2;

        double[] f_v = { 0.0, 5.0 };
        double[] f_x = { 0.0, 5.0 };
        double[] g_v = { 0.0, 5.0 };
        double[] g_x = { 0.0, 5.0 };
        int i;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Test pwl_PRODUCT_INTEGRAL on a very simple problem.");
        Console.WriteLine("  F and G are both defined over a single common");
        Console.WriteLine("  interval, so that F(X) = G(X) = X.");
        Console.WriteLine("");
        Console.WriteLine("           A           B      Integral        Exact");
        Console.WriteLine("");

        const double a = 1.0;
        for (i = 1; i <= 5; i++)
        {
            double b = i;
            double integral = ProductIntegral.pwl_product_integral(a, b, F_NUM, f_x, f_v, G_NUM,
                g_x, g_v);
            double exact = (b * b * b - a * a * a) / 3.0;
            Console.WriteLine("  " + a.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + b.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + integral.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + exact.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests pwl_PRODUCT_INTEGRAL.
        //
        //  Discussion:
        //
        //    For this test, we use multiple "pieces" for both F and G,
        //    but we define the values so that we are still actually integrating X^2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 April 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int F_NUM = 3;
        const int G_NUM = 4;

        double[] f_v = { 0.0, 2.0, 5.0 };
        double[] f_x = { 0.0, 2.0, 5.0 };
        double[] g_v = { 0.0, 1.5, 3.0, 5.0 };
        double[] g_x = { 0.0, 1.5, 3.0, 5.0 };
        int i;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  Test pwl_PRODUCT_INTEGRAL on a simple problem.");
        Console.WriteLine("  F and G are both defined over separate, multiple");
        Console.WriteLine("  intervals, but still true that F(X) = G(X) = X.");
        Console.WriteLine("");
        Console.WriteLine("           A           B      Integral        Exact");
        Console.WriteLine("");

        const double a = 1.0;
        for (i = 1; i <= 5; i++)
        {
            double b = i;
            double integral = ProductIntegral.pwl_product_integral(a, b, F_NUM, f_x, f_v, G_NUM,
                g_x, g_v);
            double exact = (b * b * b - a * a * a) / 3.0;
            Console.WriteLine("  " + a.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + b.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + integral.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + exact.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests pwl_PRODUCT_INTEGRAL.
        //
        //  Discussion:
        //
        //    For this test, F(X) and G(X) are piecewise linear interpolants to
        //    SIN(X) and 2 * COS(X), so we know the exact value of the integral
        //    of the product of the original functions, but this is only an estimate 
        //    of the exact value of the integral of the product of the interpolants.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 April 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int F_NUM = 11;
        const int G_NUM = 31;

        double[] f_v = new double[F_NUM];
        double[] f_x = new double[F_NUM];
        double[] g_v = new double[G_NUM];
        double[] g_x = new double[G_NUM];
        int i;

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  Test pwl_PRODUCT_INTEGRAL on a simple problem.");
        Console.WriteLine("  F and G are defined over separate, multiple");
        Console.WriteLine("  intervals.");
        Console.WriteLine("");
        Console.WriteLine("  F(X) interpolates SIN(X),");
        Console.WriteLine("  G(X) interpolates 2*COS(X).");
        Console.WriteLine("");
        Console.WriteLine("  We compare:");
        Console.WriteLine("");
        Console.WriteLine("  INTEGRAL, our value for the integral,");
        Console.WriteLine("  QUAD, a quadrature estimate for the integral, and");
        Console.WriteLine("  CLOSE, the value of the integral of 2*COS(X)*SIN(X)");
        Console.WriteLine("");
        Console.WriteLine("           A           B      Integral        Quad            Close");
        Console.WriteLine("");

        for (i = 0; i < F_NUM; i++)
        {
            f_x[i] = ((F_NUM - i - 1) * 0.0
                      + i * Math.PI)
                     / (F_NUM - 1);
            f_v[i] = Math.Sin(f_x[i]);
        }

        for (i = 0; i < G_NUM; i++)
        {
            g_x[i] = ((G_NUM - i - 1) * 0.0
                      + i * Math.PI)
                     / (G_NUM - 1);
            g_v[i] = 2.0 * Math.Cos(g_x[i]);
        }

        const double a = 0.0;
        for (i = 0; i <= 6; i++)
        {
            double b = i * Math.PI / 6.0;
            double integral = ProductIntegral.pwl_product_integral(a, b, F_NUM, f_x, f_v,
                G_NUM, g_x, g_v);
            double exact = -(Math.Cos(2.0 * b) - Math.Cos(2.0 * a)) / 2.0;
            const int quad_num = 2000;
            double quad = ProductIntegral.pwl_product_quad(a, b, F_NUM, f_x, f_v, G_NUM,
                g_x, g_v, quad_num);
            Console.WriteLine("  " + a.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + b.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + integral.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + exact.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests pwl_PRODUCT_INTEGRAL.
        //
        //  Discussion:
        //
        //    For this test, we compute the integrals of a hat function with itself,
        //    and a hat function with its neighbor.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 April 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int F_NUM = 3;
        const int G_NUM = 3;

        double[] f_v = { 0.0, 1.0, 0.0 };
        double[] f_x = { 0.0, 1.0, 2.0 };
        double[] g_v = { 1.0, 0.0, 0.0 };
        double[] g_x = { 0.0, 1.0, 2.0 };

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  Test pwl_PRODUCT_INTEGRAL.");
        Console.WriteLine("  The nodes are at 0, 1, and 2.");
        Console.WriteLine("  F(X) = ( 0, 1, 0 ).");
        Console.WriteLine("  G(X) = ( 1, 0, 0 ).");
        Console.WriteLine("");

        const double a = 0.0;
        const double b = 2.0;

        double integral = ProductIntegral.pwl_product_integral(a, b, F_NUM, f_x, f_v, F_NUM,
            f_x, f_v);

        Console.WriteLine("  Integral F(X) * F(X) dx = " + integral + "");

        integral = ProductIntegral.pwl_product_integral(a, b, F_NUM, f_x, f_v, G_NUM,
            g_x, g_v);

        Console.WriteLine("  Integral F(X) * G(X) dx = " + integral + "");

        integral = ProductIntegral.pwl_product_integral(a, b, G_NUM, g_x, g_v, G_NUM,
            g_x, g_v);

        Console.WriteLine("  Integral G(X) * G(X) dx = " + integral + "");
    }
}