using System;
using Burkardt.Values;

namespace TestValuesTest;

public static class BivariateTest
{
    public static void bivariate_normal_cdf_values_test ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BIVARIATE_NORMAL_CDF_VALUES_TEST tests BIVARIATE_NORMAL_CDF_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fxy = 0;
        int n_data;
        double r = 0;
        double x = 0;
        double y = 0;
        Console.WriteLine("");
        Console.WriteLine("BIVARIATE_NORMAL_CDF_VALUES_TEST:");
        Console.WriteLine("  BIVARIATE_NORMAL_CDF_VALUES stores values of");
        Console.WriteLine("  the bivariate normal CDF.");
        Console.WriteLine("");
        Console.WriteLine("      X            Y            R            F(R)(X,Y)");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Bivariate.bivariate_normal_cdf_values(ref n_data, ref x, ref y, ref r, ref fxy);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + x.ToString().PadLeft(12) + "  "
                              + y.ToString().PadLeft(12) + "  "
                              + r.ToString().PadLeft(12) + "  "
                              + fxy.ToString("0.################").PadLeft(24) + fxy + "");
        }
    }
}