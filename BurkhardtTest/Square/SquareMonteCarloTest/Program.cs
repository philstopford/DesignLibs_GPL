using System;
using System.Globalization;
using Burkardt.MonomialNS;
using Burkardt.Square;
using Burkardt.Types;

namespace SquareMonteCarloTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SQUARE_MONTE_CARLO_TEST.
        //
        //  Discussion:
        //
        //    SQUARE_MONTE_CARLO_TEST tests the SQUARE_MONTE_CARLO library.
        //    
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("SQUARE_MONTE_CARLO_TEST");
            
        Console.WriteLine("  Test the SQUARE_MONTE_CARLO library.");

        test01();
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("SQUARE_MONTE_CARLO_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");

    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 estimates integrals over the unit square in 2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 January 2014
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
        int i;
        int j;
        string cout;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Use SQUARE01_SAMPLE to estimate integrals");
        Console.WriteLine("  over the interior of the unit square in 2D.");

        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("         N        1              X^2             Y^2" + 
                          "             X^4           X^2Y^2          Y^4          X^6");
        Console.WriteLine("");

        int n = 1;

        while (n <= 65536)
        {
            double[] x = MonteCarlo.square01_sample(n, ref seed);
            cout = "  " + n.ToString().PadLeft(8);
            for (j = 0; j < 7; j++)
            {
                for (i = 0; i < 2; i++)
                {
                    e[i] = e_test[i + j * 2];
                }

                double[] value = Monomial.monomial_value(2, n, e, x);

                double result = MonteCarlo.square01_area() * typeMethods.r8vec_sum(n, value) / n;
                cout += "  " + result.ToString(CultureInfo.InvariantCulture).PadLeft(14);

            }

            Console.WriteLine(cout);

            n = 2 * n;
        }

        Console.WriteLine("");
        cout = "     Exact";
        for (j = 0; j < 7; j++)
        {
            for (i = 0; i < 2; i++)
            {
                e[i] = e_test[i + j * 2];
            }

            double   exact = MonteCarlo.square01_monomial_integral(e);
            cout += "  " + exact.ToString(CultureInfo.InvariantCulture).PadLeft(14);
        }

        Console.WriteLine(cout);

    }
}