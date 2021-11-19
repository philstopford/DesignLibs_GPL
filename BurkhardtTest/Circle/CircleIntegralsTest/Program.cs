using System;
using System.Globalization;
using Burkardt.CircleNS;
using Burkardt.MonomialNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace CircleIntegralsTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for CIRCLE_INTEGRALS_TEST.
        //
        //  Discussion:
        //
        //    CIRCLE_INTEGRALS_TEST tests the CIRCLE_INTEGRALS library.
        //    
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("CIRCLE_INTEGRALS_TEST");

        Console.WriteLine("  Test the CIRCLE_INTEGRALS library.");

        test01();

        Console.WriteLine("");
        Console.WriteLine("CIRCLE_INTEGRALS_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 uses CIRCLE01_SAMPLE with an increasing number of points.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int m = 2;
        const int n = 4192;
        int test;
        const int test_num = 20;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Use CIRCLE01_SAMPLE to compare exact and");
        Console.WriteLine("  estimated integrals along the circumference");
        Console.WriteLine("  of the unit circle in 2D.");
        //
        //  Get sample points.
        //
        int seed = 123456789;
        double[] x = Integrals.circle01_sample(n, ref seed);

        Console.WriteLine("");
        Console.WriteLine("  Number of sample points used is " + n + "");
        //
        //  Randomly choose X, Y exponents.
        //
        Console.WriteLine("");
        Console.WriteLine("  If any exponent is odd, the integral is zero.");
        Console.WriteLine("  We restrict this test to randomly chosen even exponents.");
        Console.WriteLine("");
        Console.WriteLine("  Ex  Ey     MC-Estimate           Exact      Error");
        Console.WriteLine("");

        for (test = 1; test <= test_num; test++)
        {
            int[] e = UniformRNG.i4vec_uniform_ab_new(m, 0, 5, ref seed);

            int i;
            for (i = 0; i < m; i++)
            {
                e[i] *= 2;
            }

            double[] value = Monomial.monomial_value(m, n, e, x);

            double result = Integrals.circle01_length() * typeMethods.r8vec_sum(n, value)
                            / n;
            double exact = Integrals.circle01_monomial_integral(e);
            double error = Math.Abs(result - exact);

            Console.WriteLine("  " + e[0].ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + e[1].ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + result.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + exact.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + error.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }
    }
}