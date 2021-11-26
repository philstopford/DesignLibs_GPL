using System;
using System.Globalization;
using Burkardt.PolynomialNS;

namespace LegendreShiftedPolynomialTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for LEGENDRE_SHIFTED_POLYNOMIAL_TEST.
        //
        //  Discussion:
        //
        //    LEGENDRE_SHIFTED_POLYNOMIAL_TEST tests the LEGENDRE_SHIFTED_POLYNOMIAL library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 March 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("LEGENDRE_SHIFTED_POLYNOMIAL_TEST:");
        Console.WriteLine("  Test the LEGENDRE_SHIFTED_POLYNOMIAL library.");

        p01_polynomial_value_test();
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("LEGENDRE_SHIFTED_POLYNOMIAL_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void p01_polynomial_value_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    P01_POLYNOMIAL_VALUE_TEST tests P01_POLYNOMIAL_VALUE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 March 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx1 = 0;
        int n = 0;
        double x = 0;
        double[] x_vec = new double[1];

        Console.WriteLine("");
        Console.WriteLine("P01_POLYNOMIAL_VALUE_TEST:");
        Console.WriteLine("  P01_POLYNOMIAL_VALUE evaluates the shifted Legendre polynomial P01(n,x).");
        Console.WriteLine("");
        Console.WriteLine("                        Tabulated                 Computed");
        Console.WriteLine("     N        X          P01(N,X)                 P01(N,X)                     Error");
        Console.WriteLine("");

        int n_data = 0;

        for (;;)
        {
            LegendreShifted.p01_polynomial_values(ref n_data, ref n, ref x, ref fx1);

            if (n_data == 0)
            {
                break;
            }

            x_vec[0] = x;
            double[] fx2_vec = LegendreShifted.p01_polynomial_value(1, n, x_vec);
            double fx2 = fx2_vec[n];

            double e = fx1 - fx2;

            Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + fx1.ToString(CultureInfo.InvariantCulture).PadLeft(24)
                                   + "  " + fx2.ToString(CultureInfo.InvariantCulture).PadLeft(24)
                                   + "  " + e.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
        }
    }
}