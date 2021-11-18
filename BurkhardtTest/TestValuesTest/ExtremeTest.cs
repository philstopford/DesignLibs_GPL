using System;
using Burkardt.Values;

namespace TestValuesTest;

public static class ExtremeTest
{
    public static void extreme_values_cdf_values_test ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXTREME_VALUES_CDF_VALUES_TEST tests EXTREME_VALUES_CDF_VALUES.
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
        double alpha = 0;;
        double beta = 0;
        double fx = 0;
        int n_data;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("EXTREME_VALUES_CDF_VALUES_TEST:");
        Console.WriteLine("  EXTREME_VALUES_CDF_VALUES stores values of ");
        Console.WriteLine("  the extreme values CDF.");
        Console.WriteLine("");
        Console.WriteLine("        Alpha    Beta        X                     FX");
        Console.WriteLine("");
        n_data = 0;
        for ( ; ; )
        {
            Extreme.extreme_values_cdf_values ( ref n_data, ref alpha, ref beta, ref x, ref fx );
            if ( n_data == 0 )
            {
                break;
            }
            Console.WriteLine("  "
                              + alpha.ToString(CultureInfo.InvariantCulture).PadLeft(12)                        + alpha  + "  "
                              + beta.ToString(CultureInfo.InvariantCulture).PadLeft(12)                        + beta   + "  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12)                        + x      + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }
}