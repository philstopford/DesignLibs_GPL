using System;
using Burkardt.IntegralNS;
using Burkardt.Quadrature;

namespace SquareExactnessTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SQUARE_EXACTNESS_TEST.
        //
        //  Discussion:
        //
        //    SQUARE_EXACTNESS_TEST tests the SQUARE_EXACTNESS library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {

        Console.WriteLine("");
        Console.WriteLine("SQUARE_EXACTNESS_TEST");
        Console.WriteLine("  Test the SQUARE_EXACTNESS library.");

        test01();
        test02();

        Console.WriteLine("");
        Console.WriteLine("SQUARE_EXACTNESS_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests product Gauss-Legendre rules for the Legendre 2D integral.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a = new double[2];
        double[] b = new double[2];
        int l;
        int n;
        int n_1d;
        int p_max;
        int t;
        double[] w;
        double[] x;
        double[] y;

        a[0] = -1.0;
        a[1] = -1.0;
        b[0] = +1.0;
        b[1] = +1.0;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Product Gauss-Legendre rules for the 2D Legendre integral.");
        Console.WriteLine("  Density function rho(x) = 1.");
        Console.WriteLine("  Region: -1 <= x <= +1.");
        Console.WriteLine("  Region: -1 <= y <= +1.");
        Console.WriteLine("  Level: L");
        Console.WriteLine("  Exactness: 2*L+1");
        Console.WriteLine("  Order: N = (L+1)*(L+1)");

        for (l = 0; l <= 5; l++)
        {
            n_1d = l + 1;
            n = n_1d * n_1d;
            t = 2 * l + 1;

            w = new double[n];
            x = new double[n];
            y = new double[n];

            LegendreQuadrature.legendre_2d_set(a, b, n_1d, n_1d, ref x, ref y, ref w);

            p_max = t + 1;
            Integral.legendre_2d_exactness(a, b, n, ref x, ref y, ref w, p_max);
        }

    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests Padua rules for the Legendre 2D integral.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a = new double[2];
        double[] b = new double[2];
        int l;
        int n;
        int p_max;
        double[] w;
        double[] x;
        double[] y;

        a[0] = -1.0;
        a[1] = -1.0;
        b[0] = +1.0;
        b[1] = +1.0;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  Padua rule for the 2D Legendre integral.");
        Console.WriteLine("  Density function rho(x) = 1.");
        Console.WriteLine("  Region: -1 <= x <= +1.");
        Console.WriteLine("  Region: -1 <= y <= +1.");
        Console.WriteLine("  Level: L");
        Console.WriteLine("  Exactness: L+1 when L is 0,");
        Console.WriteLine("             L   otherwise.");
        Console.WriteLine("  Order: N = (L+1)*(L+2)/2");

        for (l = 0; l <= 5; l++)
        {
            n = (l + 1) * (l + 2) / 2;

            w = new double[n];
            x = new double[n];
            y = new double[n];

            Burkardt.Values.Padua.padua_point_set(l, ref x, ref y);
            Burkardt.Values.Padua.padua_weight_set(l, ref w);

            p_max = l switch
            {
                0 => l + 2,
                _ => l + 1
            };

            Integral.legendre_2d_exactness(a, b, n, ref x, ref y, ref w, p_max);

        }

    }
}