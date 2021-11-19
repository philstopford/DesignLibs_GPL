using System;
using Burkardt.MonomialNS;
using Burkardt.SphereNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace SphereIntegralsTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SPHERE_INTEGRALS_TEST.
        //
        //  Discussion:
        //
        //    SPHERE_INTEGRALS_TEST tests the SPHERE_INTEGRALS library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("SPHERE_INTEGRALS_TEST:");
        Console.WriteLine("  Test the SPHERE_INTEGRALS library.");

        test01();

        Console.WriteLine("");
        Console.WriteLine("SPHERE_INTEGRALS_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 uses SPHERE01_SAMPLE to estimate monomial integrands.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 January 2014
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
        int n;
        double result;
        int seed;
        int test;
        const int test_num = 20;
        double[] value;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Estimate monomial integrands using Monte Carlo");
        Console.WriteLine("  over the surface of the unit sphere in 3D.");
        //
        //  Get sample points.
        //
        n = 8192;
        seed = 123456789;
        x = Integrals.sphere01_sample(n, ref seed);
        Console.WriteLine("");
        Console.WriteLine("  Number of sample points used is " + n + "");
        //
        //  Randomly choose X,Y,Z exponents between (0,0,0) and (9,9,9).
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

            result = Integrals.sphere01_area() * typeMethods.r8vec_sum(n, value) / n;
            exact = Integrals.sphere01_monomial_integral(e);
            error = Math.Abs(result - exact);

            Console.WriteLine("  " + e[0].ToString().PadLeft(2)
                                   + "  " + e[1].ToString().PadLeft(2)
                                   + "  " + e[2].ToString().PadLeft(2)
                                   + "  " + result.ToString().PadLeft(14)
                                   + "  " + exact.ToString().PadLeft(14)
                                   + "  " + error.ToString().PadLeft(14) + "");
        }
    }
}