using System;
using Burkardt.Values;

namespace TestValuesTest;

public class WeibullTest
{
    public static void weibull_cdf_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WEIBULL_CDF_VALUES_TEST tests WEIBULL_CDF_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 June 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double alpha = 0;
        ;
        double beta = 0;
        double fx = 0;
        int n_data;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("WEIBULL_CDF_VALUES_TEST:");
        Console.WriteLine("  WEIBULL_CDF_VALUES returns values of ");
        Console.WriteLine("  the Weibull Cumulative Density Function.");
        Console.WriteLine("");
        Console.WriteLine("     Alpha   Beta        X   CDF(X)");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Weibull.weibull_cdf_values(ref n_data, ref alpha, ref beta, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + alpha.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + beta.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }
}