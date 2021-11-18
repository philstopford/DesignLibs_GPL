using System;
using Burkardt.Cube;
using Burkardt.MonomialNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace CubeIntegralsTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for CUBE_INTEGRALS_TEST.
        //
        //  Discussion:
        //
        //    CUBE_INTEGRALS_TEST tests the CUBE_INTEGRALS library.
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
        Console.WriteLine("CUBE_INTEGRALS_TEST");
            
        Console.WriteLine("  Test the CUBE_INTEGRALS library.");

        test01();

        Console.WriteLine("");
        Console.WriteLine("CUBE_INTEGRALS_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 estimates integrals over the unit cube in 3D.
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
        int[] e;
        double error;
        double exact;
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
        Console.WriteLine("  Compare exact and estimated integrals");
        Console.WriteLine("  over the interior of the unit cube in 3D.");
        //
        //  Get sample points.
        //
        seed = 123456789;
        x = Integrals.cube01_sample(n, ref seed);
        Console.WriteLine("");
        Console.WriteLine("  Number of sample points is " + n + "");
        //
        //  Randomly choose exponents.
        //
        Console.WriteLine("");
        Console.WriteLine("  Ex  Ey  Ez     MC-Estimate           Exact      Error");
        Console.WriteLine("");

        for (test = 1; test <= test_num; test++)
        {
            e = UniformRNG.i4vec_uniform_ab_new(m, 0, 7, ref seed);

            value = Monomial.monomial_value(m, n, e, x);

            result = Integrals.cube01_volume() * typeMethods.r8vec_sum(n, value) / n;
            exact = Integrals.cube01_monomial_integral(e);
            error = Math.Abs(result - exact);

            Console.WriteLine("  " + e[0].ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + e[1].ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + e[2].ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + result.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + exact.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + error.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }
}