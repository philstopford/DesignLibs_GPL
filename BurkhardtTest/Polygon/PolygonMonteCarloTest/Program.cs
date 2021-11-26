using System;
using System.Globalization;
using Burkardt.MonomialNS;
using Burkardt.Polygon;
using Burkardt.Types;

namespace PolygonMonteCarloTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for POLYGON_MONTE_CARLO_TEST.
        //
        //  Discussion:
        //
        //    POLYGON_MONTE_CARLO_TEST tests the POLYGON_MONTE_CARLO library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int nv1 = 4;
        double[] v1 =
        {
            -1.0, -1.0,
            1.0, -1.0,
            1.0, 1.0,
            -1.0, 1.0
        };

        Console.WriteLine("");
        Console.WriteLine("POLYGON_MONTE_CARLO_TEST");
        Console.WriteLine("  Test the POLYGON_MONTE_CARLO library.");

        test01(nv1, v1);

        Console.WriteLine("");
        Console.WriteLine("POLYGON_MONTE_CARLO_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01(int nv, double[] v)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 estimates integrals over a polygon in 2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] e = new int[2];
        int[] e_test =
        {
            0, 0,
            2, 0,
            0, 2,
            4, 0,
            2, 2,
            0, 4,
            6, 0
        };
        int j;
        double result;
        string cout;
        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Use POLYGON_SAMPLE to estimate integrals");
        Console.WriteLine("  over the interior of a polygon in 2D.");

        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("         N" +
                          "        1" +
                          "              X^2 " +
                          "             Y^2" +
                          "             X^4" +
                          "           X^2Y^2" +
                          "             Y^4" +
                          "           X^6");

        int n = 1;

        while (n <= 65536)
        {
            double[] x = MonteCarlo.polygon_sample(nv, v, n, ref seed);

            cout = "  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(8);

            for (j = 0; j < 7; j++)
            {
                e[0] = e_test[0 + j * 2];
                e[1] = e_test[1 + j * 2];

                double[] value = Monomial.monomial_value(2, n, e, x);

                result = MonteCarlo.polygon_area(nv, v) * typeMethods.r8vec_sum(n, value) / n;
                cout += "  " + result.ToString(CultureInfo.InvariantCulture).PadLeft(14);
            }

            Console.WriteLine(cout);

            n = 2 * n;

        }

        Console.WriteLine("");
        cout = "     Exact";
        for (j = 0; j < 7; j++)
        {
            e[0] = e_test[0 + j * 2];
            e[1] = e_test[1 + j * 2];

            result = MonteCarlo.polygon_monomial_integral(nv, v, e);
            cout += "  " + result.ToString(CultureInfo.InvariantCulture).PadLeft(14);
        }

        Console.WriteLine(cout);

    }
}