using System;
using Burkardt.MonomialNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace HypersphereIntegralsTest;

using Integrals = Burkardt.HyperGeometry.Hypersphere.Integrals;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for HYPERSPHERE_INTEGRALS_TEST.
        //
        //  Discussion:
        //
        //    HYPERSPHERE_INTEGRALS_TEST tests the HYPERSPHERE_INTEGRALS library.
        //    
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("HYPERSPHERE_INTEGRALS_TEST");
        Console.WriteLine("  Test the HYPERSPHERE_INTEGRALS library.");

        test01();
        test02();

        Console.WriteLine("");
        Console.WriteLine("HYPERSPHERE_INTEGRALS_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 uses HYPERSPHERE01_SAMPLE to estimate monomial integrands in 3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] e;
        double error;
        double exact;
        int i;
        int m = 3;
        int n = 4192;
        double result;
        int seed;
        int test;
        int test_num = 20;
        double[] value;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Estimate monomial integrals using Monte Carlo");
        Console.WriteLine("  over the surface of the unit hypersphere in 3D.");
        //
        //  Get sample points.
        //
        seed = 123456789;
        typeMethods.r8vecNormalData data = new();
        x = Integrals.hypersphere01_sample(m, n, ref data, ref seed);

        Console.WriteLine("");
        Console.WriteLine("  Number of sample points used is " + n + "");
        //
        //  Randomly choose X,Y,Z exponents between 0 and 8.
        //
        Console.WriteLine("");
        Console.WriteLine("  If any exponent is odd, the integral is zero.");
        Console.WriteLine("  We will restrict this test to randomly chosen even exponents.");
        Console.WriteLine("");
        Console.WriteLine("  Ex  Ey  Ez     MC-Estimate           Exact      Error");
        Console.WriteLine("");

        for (test = 1; test <= test_num; test++)
        {
            e = UniformRNG.i4vec_uniform_ab_new(m, 0, 4, ref seed);

            for (i = 0; i < m; i++)
            {
                e[i] *= 2;
            }

            value = Monomial.monomial_value(m, n, e, x);

            result = Integrals.hypersphere01_area(m) * typeMethods.r8vec_sum(n, value)
                     / n;
            exact = Integrals.hypersphere01_monomial_integral(m, e);
            error = Math.Abs(result - exact);

            Console.WriteLine("  " + e[0]
                                   + "  " + e[1]
                                   + "  " + e[2]
                                   + "  " + result
                                   + "  " + exact
                                   + "  " + error + "");

        }
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 uses HYPERSPHERE01_SAMPLE to estimate monomial integrands in 6D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] e;
        double error;
        double exact;
        int i;
        int m = 6;
        int n = 4192;
        double result;
        int seed;
        int test;
        int test_num = 20;
        double[] value;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  Estimate monomial integrals using Monte Carlo");
        Console.WriteLine("  over the surface of the unit hypersphere in 6D.");
        //
        //  Get sample points.
        //
        seed = 123456789;
        typeMethods.r8vecNormalData data = new();
        x = Integrals.hypersphere01_sample(m, n, ref data, ref seed);

        Console.WriteLine("");
        Console.WriteLine("  Number of sample points used is " + n + "");
        //
        //  Randomly choose X,Y,Z exponents between 0 and 6.
        //
        Console.WriteLine("");
        Console.WriteLine("  If any exponent is odd, the integral is zero.");
        Console.WriteLine("  We will restrict this test to randomly chosen even exponents.");
        Console.WriteLine("");
        Console.WriteLine("  E1  E2  E3  E4  E5  E6     MC-Estimate           Exact      Error");
        Console.WriteLine("");

        for (test = 1; test <= test_num; test++)
        {
            e = UniformRNG.i4vec_uniform_ab_new(m, 0, 3, ref seed);

            for (i = 0; i < m; i++)
            {
                e[i] *= 2;
            }

            value = Monomial.monomial_value(m, n, e, x);

            result = Integrals.hypersphere01_area(m) * typeMethods.r8vec_sum(n, value)
                     / n;
            exact = Integrals.hypersphere01_monomial_integral(m, e);
            error = Math.Abs(result - exact);

            Console.WriteLine("  " + e[0].ToString().PadLeft(2)
                                   + "  " + e[1].ToString().PadLeft(2)
                                   + "  " + e[2].ToString().PadLeft(2)
                                   + "  " + e[3].ToString().PadLeft(2)
                                   + "  " + e[4].ToString().PadLeft(2)
                                   + "  " + e[5].ToString().PadLeft(2)
                                   + "  " + result.ToString().PadLeft(14)
                                   + "  " + exact.ToString().PadLeft(14)
                                   + "  " + error.ToString().PadLeft(10) + "");
        }
    }
}