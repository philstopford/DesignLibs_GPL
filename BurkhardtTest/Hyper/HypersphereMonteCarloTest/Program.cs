using System;
using Burkardt.MonomialNS;
using Burkardt.Types;

namespace HypersphereMonteCarloTest;

using MonteCarlo = Burkardt.HyperGeometry.Hypersphere.MonteCarlo;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for HYPERSPHERE_MONTE_CARLO_TEST.
        //
        //  Discussion:
        //
        //    HYPERSPHERE_MONTE_CARLO_TEST tests the HYPERSPHERE_MONTE_CARLO library.
        //    
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("HYPERSPHERE_MONTE_CARLO_TEST");
        Console.WriteLine("  Test the HYPERSPHERE_MONTE_CARLO library.");

        test01();
        test02();

        Console.WriteLine("");
        Console.WriteLine("HYPERSPHERE_MONTE_CARLO_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 uses HYPERSPHERE01_SAMPLE to estimate integrals in 3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 January 2014
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
        int m = 3;
        int n;
        double result;
        int seed;
        double[] value;
        double[] x;
        string cout = "";

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Use HYPERSPHERE01_SAMPLE to estimate integrals");
        Console.WriteLine("  on the surface of the unit hypersphere in M dimensions.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension M = " + m + "");

        seed = 123456789;
        typeMethods.r8vecNormalData data = new();

        Console.WriteLine("");
        Console.WriteLine("         N        1              X^2             Y^2" +
                          "             Z^2             X^4            X^2Y^2          Z^4");
        Console.WriteLine("");

        n = 1;

        while (n <= 65536)
        {
            x = MonteCarlo.hypersphere01_sample(m, n, ref data, ref seed);
            cout = "  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(8);
            for (j = 0; j < 7; j++)
            {
                for (i = 0; i < m; i++)
                {
                    e[i] = e_test[i + j * m];
                }

                value = Monomial.monomial_value(m, n, e, x);

                result = MonteCarlo.hypersphere01_area(m) * typeMethods.r8vec_sum(n, value) / n;
                cout += "  " + result.ToString("0.##########").PadLeft(14);

            }

            Console.WriteLine(cout);

            n = 2 * n;
        }

        Console.WriteLine("");
        cout = "  " + "   Exact";
        for (j = 0; j < 7; j++)
        {
            for (i = 0; i < m; i++)
            {
                e[i] = e_test[i + j * m];
            }

            result = MonteCarlo.hypersphere01_monomial_integral(m, e);
            cout += "  " + result.ToString("0.##########").PadLeft(14);
        }

        Console.WriteLine(cout);

    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 uses HYPERSPHERE01_SAMPLE to estimate integrals in 6D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 January 2014
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
        int m = 6;
        int n;
        double result;
        int seed;
        double[] value;
        double[] x;
        string cout = "";

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  Use HYPERSPHERE01_SAMPLE to estimate integrals");
        Console.WriteLine("  on the surface of the unit hypersphere in M dimensions.");
        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension M = " + m + "");

        seed = 123456789;
        typeMethods.r8vecNormalData data = new();

        Console.WriteLine("");
        Console.WriteLine("         N" +
                          "        1      " +
                          "        U      " +
                          "         V^2   " +
                          "         V^2W^2" +
                          "         X^4   " +
                          "         Y^2Z^2" +
                          "         Z^6");
        Console.WriteLine("");

        n = 1;

        while (n <= 65536)
        {
            x = MonteCarlo.hypersphere01_sample(m, n, ref data, ref seed);
            cout = "  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(8);
            for (j = 0; j < 7; j++)
            {
                for (i = 0; i < m; i++)
                {
                    e[i] = e_test[i + j * m];
                }

                value = Monomial.monomial_value(m, n, e, x);

                result = MonteCarlo.hypersphere01_area(m) * typeMethods.r8vec_sum(n, value) / n;
                cout += "  " + result.ToString("0.##########").PadLeft(14);
            }

            Console.WriteLine(cout);

            n = 2 * n;
        }

        Console.WriteLine("");
        cout = "  " + "   Exact";
        for (j = 0; j < 7; j++)
        {
            for (i = 0; i < m; i++)
            {
                e[i] = e_test[i + j * m];
            }

            result = MonteCarlo.hypersphere01_monomial_integral(m, e);
            cout += "  " + result.ToString("0.##########").PadLeft(14);
        }

        Console.WriteLine(cout);
    }
}