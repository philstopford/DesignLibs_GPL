using System;
using System.Globalization;
using Burkardt.MonomialNS;
using Burkardt.PyramidNS;
using Burkardt.Types;

namespace PyramidMonteCarloTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for PYRAMID_MONTE_CARLO_TEST.
        //
        //  Discussion:
        //
        //    PYRAMID_MONTE_CARLO_TEST tests the PYRAMID_MONTE_CARLO library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 April 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("PYRAMID_MONTE_CARLO_TEST");
        Console.WriteLine("  Test the PYRAMID_MONTE_CARLO library.");

        test01();

        Console.WriteLine("");
        Console.WriteLine("PYRAMID_MONTE_CARLO_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 compares exact and estimated monomial integrals.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 April 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int M = 3;
        const int TEST_NUM = 10;

        int[] e = new int[M];
        int[] e_test =
        {
            0, 0, 0,
            0, 0, 1,
            2, 0, 0,
            0, 2, 0,
            0, 0, 2,
            2, 0, 1,
            0, 2, 1,
            0, 0, 3,
            2, 2, 0,
            2, 0, 2
        };
        int i;
        int j;
        double result;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Use PYRAMID01_SAMPLE to estimate integrals");
        Console.WriteLine("  over the interior of the unit pyramid in 3D.");

        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("         N"
                          + "        1"
                          + "               Z"
                          + "             X^2"
                          + "             Y^2"
                          + "             Z^2"
                          + "            X^2Z"
                          + "            Y^2Z"
                          + "             Z^3"
                          + "          X^2Y^2"
                          + "          X^2Z^2");
        Console.WriteLine("");

        int n = 1;

        while (n <= 65536)
        {
            string cout = "  " + n.ToString().PadLeft(8);

            double[] x = MonteCarlo.pyramid01_sample(n, ref seed);

            for (j = 0; j < TEST_NUM; j++)
            {
                for (i = 0; i < M; i++)
                {
                    e[i] = e_test[i + j * M];
                }

                double[] value = Monomial.monomial_value(M, n, e, x);

                result = MonteCarlo.pyramid01_volume() * typeMethods.r8vec_sum(n, value) / n;
                cout += "  " + result.ToString(CultureInfo.InvariantCulture).PadLeft(14);
            }

            Console.WriteLine(cout);

            n = 2 * n;
        }

        Console.WriteLine("");
        string cout2 = "     Exact";

        for (j = 0; j < 10; j++)
        {
            for (i = 0; i < M; i++)
            {
                e[i] = e_test[i + j * M];
            }

            result = MonteCarlo.pyramid01_integral(e);
            cout2 += "  " + result.ToString(CultureInfo.InvariantCulture).PadLeft(14);
        }

        Console.WriteLine(cout2);
    }
}