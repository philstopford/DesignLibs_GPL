using System;
using Burkardt.MonomialNS;
using Burkardt.SimplexNS;
using Burkardt.Types;

namespace SimplexMonteCarloTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SIMPLEX_MONTE_CARLO_TEST.
        //
        //  Discussion:
        //
        //    SIMPLEX_MONTE_CARLO_TEST tests the SIMPLEX_MONTE_CARLO library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("SIMPLEX_MONTE_CARLO_TEST");
        Console.WriteLine("  Test the SIMPLEX_MONTE_CARLO library.");

        test01();
        test02();
        test03();

        Console.WriteLine("");
        Console.WriteLine("SIMPLEX_MONTE_CARLO_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 uses SIMPLEX_UNIT_SAMPLE to estimate integrals in 3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] e = new int [3];
        int[] e_test =
        {
            0, 0, 0,
            1, 0, 0,
            0, 1, 0,
            0, 0, 1,
            2, 0, 0,
            1, 1, 0,
            1, 0, 1,
            0, 2, 0,
            0, 1, 1,
            0, 0, 2
        };
        int i;
        int j;
        const int m = 3;
        int n;
        double result;
        int seed;
        double[] value;
        double[] x;
        string cout = "";

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Use SIMPLEX_UNIT_SAMPLE for a Monte Carlo estimate of an");
        Console.WriteLine("  integral over the interior of the unit simplex in 3D.");

        seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("         N        1               X               Y " +
                          "              Z               X^2              XY             XZ" +
                          "              Y^2             YZ               Z^2");
        Console.WriteLine("");

        n = 1;

        while (n <= 65536)
        {
            x = MonteCarlo.simplex_unit_sample(m, n, ref seed);

            cout = "  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(8);
            for (j = 0; j < 10; j++)
            {
                for (i = 0; i < m; i++)
                {
                    e[i] = e_test[i + j * m];
                }

                value = Monomial.monomial_value(m, n, e, x);

                result = MonteCarlo.simplex_unit_volume(m) * typeMethods.r8vec_sum(n, value)
                         / n;

                cout += "  " + result.ToString(CultureInfo.InvariantCulture).PadLeft(14);

            }

            Console.WriteLine(cout);

            n = 2 * n;

        }

        Console.WriteLine("");
        cout = "     Exact";

        for (j = 0; j < 10; j++)
        {
            for (i = 0; i < m; i++)
            {
                e[i] = e_test[i + j * m];
            }

            result = MonteCarlo.simplex_unit_monomial_integral(m, e);
            cout += "  " + result.ToString(CultureInfo.InvariantCulture).PadLeft(14);
        }

        Console.WriteLine(cout);
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 uses SIMPLEX_UNIT_SAMPLE to estimate integrals in 6D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        // 
        //    16 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] e = new int[6];
        int[] e_test =
        {
            0, 0, 0, 0, 0, 0,
            1, 0, 0, 0, 0, 0,
            0, 2, 0, 0, 0, 0,
            0, 2, 2, 0, 0, 0,
            0, 0, 0, 4, 0, 0,
            2, 0, 0, 0, 2, 2,
            0, 0, 0, 0, 0, 6
        };
        int i;
        int j;
        const int m = 6;
        int n;
        double result;
        int seed;
        double[] value;
        double[] x;
        string cout;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  Use SIMPLEX_UNIT_SAMPLE for a Monte Carlo estimate of an");
        Console.WriteLine("  integral over the interior of the unit simplex in 6D.");

        seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("         N" +
                          "        1      " +
                          "        U      " +
                          "         V^2    " +
                          "         V^2W^2 " +
                          "         X^4    " +
                          "         Y^2Z^2 " +
                          "         Z^6");
        Console.WriteLine("");

        n = 1;

        while (n <= 65536)
        {
            x = MonteCarlo.simplex_unit_sample(m, n, ref seed);

            cout = "  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(8);
            for (j = 0; j < 7; j++)
            {
                for (i = 0; i < m; i++)
                {
                    e[i] = e_test[i + j * m];
                }

                value = Monomial.monomial_value(m, n, e, x);

                result = MonteCarlo.simplex_unit_volume(m) * typeMethods.r8vec_sum(n, value)
                         / n;

                cout += "  " + result.ToString(CultureInfo.InvariantCulture).PadLeft(14);

            }

            Console.WriteLine(cout);

            n = 2 * n;
        }

        Console.WriteLine("");
        cout = "     Exact";

        for (j = 0; j < 7; j++)
        {
            for (i = 0; i < m; i++)
            {
                e[i] = e_test[i + j * m];
            }

            result = MonteCarlo.simplex_unit_monomial_integral(m, e);
            cout += "  " + result.ToString(CultureInfo.InvariantCulture).PadLeft(14);
        }

        Console.WriteLine(cout);
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 uses SIMPLEX_GENERAL_SAMPLE to estimate integrals in 3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 March 2017
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
            1, 0, 0,
            0, 1, 0,
            0, 0, 1,
            2, 0, 0,
            1, 1, 0,
            1, 0, 1,
            0, 2, 0,
            0, 1, 1,
            0, 0, 2
        };
        int i;
        int j;
        const int m = 3;
        int n;
        double result;
        int seed;
        double[] t =
        {
            1.0, 0.0, 0.0,
            2.0, 0.0, 0.0,
            1.0, 2.0, 0.0,
            1.0, 0.0, 3.0
        };
        double[] value;
        double[] x;
        string cout = "";

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  SIMPLEX_GENERAL_SAMPLE computes a Monte Carlo estimate of an");
        Console.WriteLine("  integral over the interior of a general simplex in 3D.");

        Console.WriteLine("");
        Console.WriteLine("  Simplex vertices:");
        Console.WriteLine("");
        for (j = 0; j < 4; j++)
        {
            cout = "";
            for (i = 0; i < 3; i++)
            {
                cout += "  " + t[i + j * 3].ToString(CultureInfo.InvariantCulture).PadLeft(14);
            }

            Console.WriteLine(cout);
        }

        seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("         N        1               X               Y " +
                          "              Z               X^2              XY             XZ" +
                          "              Y^2             YZ               Z^2");
        Console.WriteLine("");

        n = 1;

        while (n <= 65536)
        {
            x = MonteCarlo.simplex_general_sample(m, n, t, ref seed);

            cout = "  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(8);
            for (j = 0; j < 10; j++)
            {
                for (i = 0; i < m; i++)
                {
                    e[i] = e_test[i + j * m];
                }

                value = Monomial.monomial_value(m, n, e, x);

                result = MonteCarlo.simplex_general_volume(m, t) * typeMethods.r8vec_sum(n, value)
                         / n;

                cout += "  " + result.ToString(CultureInfo.InvariantCulture).PadLeft(14);
            }

            Console.WriteLine(cout);

            n = 2 * n;
        }

    }
}