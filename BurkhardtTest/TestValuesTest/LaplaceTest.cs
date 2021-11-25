using System;
using System.Globalization;
using Burkardt.Values;

namespace TestValuesTest;

public static class LaplaceTest
{
    public static void laplace_cdf_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAPLACE_CDF_VALUES_TEST tests LAPLACE_CDF_VALUES.
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
        double beta = 0;
        double fx = 0;
        double mu = 0;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("LAPLACE_CDF_VALUES_TEST:");
        Console.WriteLine("  LAPLACE_CDF_VALUES returns values of ");
        Console.WriteLine("  the Laplace Cumulative Density Function.");
        Console.WriteLine("");
        Console.WriteLine("     Mu      Beta         X   CDF(X)");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Laplace.laplace_cdf_values(ref n_data, ref mu, ref beta, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + mu.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + beta.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }
}