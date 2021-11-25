using System;
using System.Globalization;
using Burkardt.Values;

namespace TestValuesTest;

public static class GeometricTest
{
    public static void geometric_cdf_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEOMETRIC_CDF_VALUES_TEST tests GEOMETRIC_CDF_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double cdf = 0;
        double p = 0;
        int x = 0;
        Console.WriteLine("");
        Console.WriteLine("GEOMETRIC_CDF_VALUES_TEST:");
        Console.WriteLine("  GEOMETRIC_CDF_VALUES stores values of");
        Console.WriteLine("  the Geometric Probability Cumulative Density Function.");
        Console.WriteLine("");
        Console.WriteLine("      X      P       CDF");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Geometric.geometric_cdf_values(ref n_data, ref x, ref p, ref cdf);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + x.ToString().PadLeft(6) + "  "
                              + p.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + cdf.ToString("0.################").PadLeft(24) + "");
        }
    }
}