using System;
using Burkardt.Values;

namespace TestValuesTest;

public class RayleighTest
{
    public static void rayleigh_cdf_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RAYLEIGH_CDF_VALUES_TEST tests RAYLEIGH_CDF_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 January 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        int n_data;
        double sigma = 0;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("RAYLEIGH_CDF_VALUES_TEST:");
        Console.WriteLine("  RAYLEIGH_CDF_VALUES stores values of");
        Console.WriteLine("  the Rayleigh CDF.");
        Console.WriteLine("");
        Console.WriteLine("      SIGMA        X            CDF(X)");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Rayleigh.rayleigh_cdf_values(ref n_data, ref sigma, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + sigma.ToString().PadLeft(12) + "  "
                              + x.ToString().PadLeft(12) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }
}