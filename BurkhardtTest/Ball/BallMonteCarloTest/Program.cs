using System;
using Burkardt.Ball;
using Burkardt.MonomialNS;
using Burkardt.Types;

namespace BallMonteCarloTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for BALL_MONTE_CARLO_TEST.
        //
        //  Discussion:
        //
        //    BALL_MONTE_CARLO_TEST tests the BALL_MONTE_CARLO library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("BALL_MONTE_CARLO_TEST");
        Console.WriteLine("  Test the BALL_MONTE_CARLO library.");

        test01();

        Console.WriteLine("");
        Console.WriteLine("BALL_MONTE_CARLO_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 uses BALL01_SAMPLE with an increasing number of points.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] e = new int[3];
        int[] e_test =
        {
            0, 0, 0,
            2, 0, 0,
            0, 2, 0,
            0, 0, 2,
            4, 0, 0,
            2, 2, 0,
            0, 0, 4
        };
        int i;
        int j;
        int n;
        double[] result = new double[7];
        int seed;
        double[] value;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Estimate integrals over the interior of the unit ball");
        Console.WriteLine("  using the Monte Carlo method.");

        seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("         N        1              X^2             Y^2 " +
                          "             Z^2             X^4           X^2Y^2           Z^4");
        Console.WriteLine("");

        n = 1;

        string cout = "";

        while (n <= 65536)
        {
            x = MonteCarlo.ball01_sample(n, ref seed);

            cout = "  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(8);
            for (j = 0; j < 7; j++)
            {
                for (i = 0; i < 3; i++)
                {
                    e[i] = e_test[i + j * 3];
                }

                value = Monomial.monomial_value(3, n, e, x);

                result[j] = MonteCarlo.ball01_volume() * typeMethods.r8vec_sum(n, value) / n;
                cout += "  " + result[j].ToString(CultureInfo.InvariantCulture).PadLeft(14);

            }

            Console.WriteLine(cout);

            n = 2 * n;
        }

        Console.WriteLine("");
        cout = "     Exact";
        for (j = 0; j < 7; j++)
        {
            for (i = 0; i < 3; i++)
            {
                e[i] = e_test[i + j * 3];
            }

            result[j] = MonteCarlo.ball01_monomial_integral(e);
            cout += "  " + result[j].ToString(CultureInfo.InvariantCulture).PadLeft(14);
        }

        Console.WriteLine(cout);

    }
}