using System;
using System.Globalization;
using Burkardt.MonomialNS;
using Burkardt.Types;
using Integrals = Burkardt.Wedge.Integrals;

namespace WedgeIntegralsTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for WEDGE_INTEGRALS_TEST.
        //
        //  Discussion:
        //
        //    WEDGE_INTEGRALS_TEST tests the WEDGE_INTEGRALS library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("WEDGE_INTEGRALS_TEST");
            
        Console.WriteLine("  Test the WEDGE_INTEGRALS library.");

        test01();

        Console.WriteLine("");
        Console.WriteLine("WEDGE_INTEGRALS_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 compares exact and estimated monomial integrals.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int e_max = 6;
        int e3;
        int[] expon = new int[3];
        const int m = 3;
        const int n = 500000;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Compare exact and estimated integrals ");
        Console.WriteLine("  over the unit wedge in 3D.");
        //
        //  Get sample points.
        //
        int seed = 123456789;
        double[] x = Integrals.wedge01_sample(n, ref seed);

        Console.WriteLine("");
        Console.WriteLine("  Number of sample points used is " + n + "");
        Console.WriteLine("");
        Console.WriteLine("   E1  E2  E3     MC-Estimate      Exact           Error");
        Console.WriteLine("");
        //
        //  Check all monomials up to total degree E_MAX.
        //
        for (e3 = 0; e3 <= e_max; e3++)
        {
            expon[2] = e3;
            int e2;
            for (e2 = 0; e2 <= e_max - e3; e2++)
            {
                expon[1] = e2;
                int e1;
                for (e1 = 0; e1 <= e_max - e3 - e2; e1++)
                {
                    expon[0] = e1;

                    double[] value = Monomial.monomial_value(m, n, expon, x);

                    double q = Integrals.wedge01_volume() * typeMethods.r8vec_sum(n, value) / n;
                    double exact = Integrals.wedge01_integral(expon);
                    double error = Math.Abs(q - exact);

                    Console.WriteLine(expon[0].ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                                                                     + expon[1].ToString(CultureInfo.InvariantCulture).PadLeft(2) + "  "
                                                                     + expon[2].ToString(CultureInfo.InvariantCulture).PadLeft(2) + "  "
                                                                     + q.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "  "
                                                                     + exact.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "  "
                                                                     + error.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                }
            }
        }
    }
}