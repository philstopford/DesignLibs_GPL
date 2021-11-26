using System;
using System.Globalization;
using Burkardt.MonomialNS;
using Burkardt.Types;

namespace EllipsoidMonteCarloTest;

using MonteCarlo = Burkardt.Ellipsoid.MonteCarlo;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for ELLIPSOID_MONTE_CARLO_TEST.
        //
        //  Discussion:
        //
        //    ELLIPSOID_MONTE_CARLO_TEST tests the ELLIPSOID_MONTE_CARLO library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("ELLIPSOID_MONTE_CARLO_TEST");
        Console.WriteLine("  Test the ELLIPSOID_MONTE_CARLO library.");

        test01();
        test02();
        test03();
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("ELLIPSOID_MONTE_CARLO_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 uses ELLIPSOID_SAMPLE on a 2D ellipse centered at (0,0).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int M = 2;

        double[] a =  {
                9.0, 1.0,
                1.0, 4.0
            }
            ;
        int[] e = new int[M];
        int[] e_test =  {
                0, 0,
                1, 0,
                0, 1,
                2, 0,
                1, 1,
                0, 2,
                3, 0
            }
            ;
        const double r = 2.0;
        double[] v =  {
                0.0, 0.0
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Use ELLIPSOID_SAMPLE to estimate integrals");
        Console.WriteLine("  in a 2D ellipse x' * A * x <= r^2.");

        Console.WriteLine("");
        typeMethods.r8_print(r, "  Ellipsoid radius R:");
        typeMethods.r8vec_print(M, v, "  Ellipsoid center V:");
        typeMethods.r8mat_print(M, M, a, "  Ellipsoid matrix A:");

        double volume = MonteCarlo.ellipsoid_volume(M, a, v, r);
        Console.WriteLine("");
        typeMethods.r8_print(volume, "  Ellipsoid volume:");

        int seed = 123456789;
        typeMethods.r8vecNormalData data = new();

        Console.WriteLine("");
        Console.WriteLine("         N        1              X               Y  " +
                          "             X^2               XY             Y^2             X^3");
        Console.WriteLine("");

        int n = 1;

        while (n <= 65536)
        {
            double[] x = MonteCarlo.ellipsoid_sample(M, n, a, v, r, ref data, ref seed);

            string cout = n.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  ";
            int j;
            for (j = 0; j < 7; j++)
            {
                int i;
                for (i = 0; i < M; i++)
                {
                    e[i] = e_test[i + j * M];
                }

                double[] value = Monomial.monomial_value(M, n, e, x);

                double result = volume * typeMethods.r8vec_sum(n, value) / n;
                cout += result.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "  ";
            }

            Console.WriteLine(cout);
                
            n = 2 * n;
        }
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 uses ELLIPSOID_SAMPLE on a 2D ellipse centered at (2,3).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int M = 2;

        double[] a =  {
                9.0, 1.0,
                1.0, 4.0
            }
            ;
        int[] e = new int[M];
        int[] e_test =  {
                0, 0,
                1, 0,
                0, 1,
                2, 0,
                1, 1,
                0, 2,
                3, 0
            }
            ;
        const double r = 0.5;
        double[] v =  {
                2.0, 3.0
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  Use ELLIPSOID_SAMPLE to estimate integrals");
        Console.WriteLine("  in a 2D ellipse (x-v)' * A * (x-v) <= r^2.");

        Console.WriteLine("");
        typeMethods.r8_print(r, "  Ellipsoid radius R:");
        typeMethods.r8vec_print(M, v, "  Ellipsoid center V:");
        typeMethods.r8mat_print(M, M, a, "  Ellipsoid matrix A:");

        double volume = MonteCarlo.ellipsoid_volume(M, a, v, r);
        Console.WriteLine("");
        typeMethods.r8_print(volume, "  Ellipsoid volume:");

        int seed = 123456789;
        typeMethods.r8vecNormalData data = new();

        Console.WriteLine("");
        Console.WriteLine("         N        1              X               Y  " + 
                          "             X^2               XY             Y^2             X^3");
        Console.WriteLine("");

        int n = 1;

        while (n <= 65536)
        {
            double[] x = MonteCarlo.ellipsoid_sample(M, n, a, v, r, ref data, ref seed);

            string cout = n.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  ";
            int j;
            for (j = 0; j < 7; j++)
            {
                int i;
                for (i = 0; i < M; i++)
                {
                    e[i] = e_test[i + j * M];
                }

                double[] value = Monomial.monomial_value(M, n, e, x);

                double result = volume * typeMethods.r8vec_sum(n, value) / n;
                cout += result.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "  ";
            }

            Console.WriteLine(cout);

            n = 2 * n;
        }
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 uses ELLIPSOID_SAMPLE on a 3D ellipse centered at (1,2,3).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int M = 3;

        double[] a =  {
                9.0, 6.0, 3.0,
                6.0, 5.0, 4.0,
                3.0, 4.0, 9.0
            }
            ;
        int[] e = new int[M];
        int[] e_test =  {
                0, 0, 0,
                1, 0, 0,
                0, 1, 0,
                0, 0, 1,
                2, 0, 0,
                0, 2, 2,
                0, 0, 3
            }
            ;
        const double r = 0.5;
        double[] v =  {
                1.0, 2.0, 3.0
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  Use ELLIPSOID_SAMPLE to estimate integrals");
        Console.WriteLine("  in a 3D ellipse (x-v)' * A * (x-v) <= r^2.");

        Console.WriteLine("");
        typeMethods.r8_print(r, "  Ellipsoid radius R:");
        typeMethods.r8vec_print(M, v, "  Ellipsoid center V:");
        typeMethods.r8mat_print(M, M, a, "  Ellipsoid matrix A:");

        double volume = MonteCarlo.ellipsoid_volume(M, a, v, r);
        Console.WriteLine("");
        typeMethods.r8_print(volume, "  Ellipsoid volume:");

        int seed = 123456789;
        typeMethods.r8vecNormalData data = new();

        Console.WriteLine("");
        Console.WriteLine("         N        1              X               Y  " +
                          "              Z                X^2            YZ              Z^3");
        Console.WriteLine("");

        int n = 1;

        while (n <= 65536)
        {
            double[] x = MonteCarlo.ellipsoid_sample(M, n, a, v, r, ref data, ref seed);

            string cout = n.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  ";
            int j;
            for (j = 0; j < 7; j++)
            {
                int i;
                for (i = 0; i < M; i++)
                {
                    e[i] = e_test[i + j * M];
                }

                double[] value = Monomial.monomial_value(M, n, e, x);

                double result = volume * typeMethods.r8vec_sum(n, value) / n;
                cout += result.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "  ";
            }

            Console.WriteLine(cout);

            n = 2 * n;
        }
    }
}