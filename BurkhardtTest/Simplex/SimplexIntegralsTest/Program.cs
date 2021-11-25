using System;
using System.Globalization;
using Burkardt.MonomialNS;
using Burkardt.SimplexNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace SimplexIntegralsTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SIMPLEX_INTEGRALS_TEST.
        //
        //  Discussion:
        //
        //    SIMPLEX_INTEGRALS_TEST tests the SIMPLEX_INTEGRALS library.
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
        Console.WriteLine("SIMPLEX_INTEGRALS_TEST");
        Console.WriteLine("  Test the SIMPLEX_INTEGRALS library.");

        test01();
        test02();

        Console.WriteLine("");
        Console.WriteLine("SIMPLEX_INTEGRALS_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 compares exact and estimated integrals in 3D.
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
        const int m = 3;
        const int n = 4192;
        int test;
        const int test_num = 20;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Estimate monomial integrals using Monte Carlo");
        Console.WriteLine("  over the interior of the unit simplex in M dimensions.");
        //
        //  Get sample points.
        //
        int seed = 123456789;
        double[] x = Integrals.simplex01_sample(m, n, ref seed);

        Console.WriteLine("");
        Console.WriteLine("  Number of sample points used is " + n + "");
        //
        //  Randomly choose exponents.
        //
        Console.WriteLine("");
        Console.WriteLine("  We randomly choose the exponents.");
        Console.WriteLine("");
        Console.WriteLine("  Ex  Ey  Ez     MC-Estimate      Exact           Error");
        Console.WriteLine("");

        for (test = 1; test <= test_num; test++)
        {
            int[] e = UniformRNG.i4vec_uniform_ab_new(m, 0, 4, ref seed);

            double[] value = Monomial.monomial_value(m, n, e, x);

            double result = Integrals.simplex01_volume(m) * typeMethods.r8vec_sum(n, value)
                            / n;

            double exact = Integrals.simplex01_monomial_integral(m, e);
            double error = Math.Abs(result - exact);

            Console.WriteLine("  " + e[0].ToString().PadLeft(2)
                                   + "  " + e[1].ToString().PadLeft(2)
                                   + "  " + e[2].ToString().PadLeft(2)
                                   + "  " + result.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + exact.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + error.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 compares exact and estimated integrals in 6D.
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
        const int m = 6;
        const int n = 4192;
        int test;
        const int test_num = 20;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  Estimate monomial integrals using Monte Carlo");
        Console.WriteLine("  over the interior of the unit simplex in M dimensions.");
        //
        //  Get sample points.
        //
        int seed = 123456789;
        double[] x = Integrals.simplex01_sample(m, n, ref seed);

        Console.WriteLine("");
        Console.WriteLine("  Number of sample points used is " + n + "");
        //
        //  Randomly choose exponents.
        //
        Console.WriteLine("");
        Console.WriteLine("  We randomly choose the exponents.");
        Console.WriteLine("");
        Console.WriteLine("  E1  E2  E3  E4  E5  E6     MC-Estimate      Exact           Error");
        Console.WriteLine("");

        for (test = 1; test <= test_num; test++)
        {
            int[] e = UniformRNG.i4vec_uniform_ab_new(m, 0, 4, ref seed);

            double[] value = Monomial.monomial_value(m, n, e, x);

            double result = Integrals.simplex01_volume(m) * typeMethods.r8vec_sum(n, value)
                            / n;

            double exact = Integrals.simplex01_monomial_integral(m, e);
            double error = Math.Abs(result - exact);

            Console.WriteLine("  " + e[0].ToString().PadLeft(2)
                                   + "  " + e[1].ToString().PadLeft(2)
                                   + "  " + e[2].ToString().PadLeft(2)
                                   + "  " + e[3].ToString().PadLeft(2)
                                   + "  " + e[4].ToString().PadLeft(2)
                                   + "  " + e[5].ToString().PadLeft(2)
                                   + "  " + result.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + exact.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + error.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

    }
}