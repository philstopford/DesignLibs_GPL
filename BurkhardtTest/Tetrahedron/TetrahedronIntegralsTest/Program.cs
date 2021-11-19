using System;
using Burkardt.MonomialNS;
using Burkardt.TetrahedronNS;
using Burkardt.Types;

namespace TetrahedronIntegralsTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TETRAHEDRON_INTEGRALS_TEST.
        //
        //  Discussion:
        //
        //    TETRAHEDRON_INTEGRALS_TEST tests the TETRAHEDRON_INTEGRALS library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("TETRAHEDRON_INTEGRALS_TEST");
        Console.WriteLine("  Test the TETRAHEDRON_INTEGRALS library.");

        test01();

        Console.WriteLine("");
        Console.WriteLine("TETRAHEDRON_INTEGRALS_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 uses TETRAHEDRON_SAMPLE_01 to compare exact and estimated integrals.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] e = new int[3];
        double error;
        double exact;
        int i;
        int j;
        int k;
        int m = 3;
        int n = 4192;
        double result;
        int seed;
        double[] value;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Estimate monomial integrals using Monte Carlo");
        Console.WriteLine("  over the interior of the unit tetrahedron in 3D.");
        //
        //  Get sample points.
        //
        seed = 123456789;
        x = Integrals.tetrahedron01_sample(n, ref seed);

        Console.WriteLine("");
        Console.WriteLine("  Number of sample points used is " + n + "");
        //
        // Run through the exponents.
        //
        Console.WriteLine("");
        Console.WriteLine("  Ex  Ey  Ez     MC-Estimate      Exact           Error");
        Console.WriteLine("");

        for (i = 0; i <= 3; i++)
        {
            e[0] = i;
            for (j = 0; j <= 3; j++)
            {
                e[1] = j;
                for (k = 0; k <= 3; k++)
                {
                    e[2] = k;

                    value = Monomial.monomial_value(m, n, e, x);

                    result = Integrals.tetrahedron01_volume() * typeMethods.r8vec_sum(n, value)
                             / n;
                    exact = Integrals.tetrahedron01_monomial_integral(e);
                    error = Math.Abs(result - exact);

                    Console.WriteLine("  " + e[0].ToString().PadLeft(2)
                                           + "  " + e[1].ToString().PadLeft(2)
                                           + "  " + e[2].ToString().PadLeft(2)
                                           + "  " + result.ToString().PadLeft(14)
                                           + "  " + exact.ToString().PadLeft(14)
                                           + "  " + error.ToString().PadLeft(10) + "");

                }
            }
        }
    }
}