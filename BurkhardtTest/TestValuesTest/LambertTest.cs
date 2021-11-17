using System;
using Burkardt.Values;

namespace TestValuesTest;

public static class LambertTest
{
    public static void lambert_w_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAMBERT_W_VALUES_TEST tests LAMBERT_W_VALUES.
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
        double fx = 0;
        int n_data;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("LAMBERT_W_VALUES_TEST:");
        Console.WriteLine("  LAMBERT_W_VALUES stores values of ");
        Console.WriteLine("  the Lambert W function.");
        Console.WriteLine("");
        Console.WriteLine("                X                     W(X)");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Lambert.lambert_w_values(ref n_data, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + x.ToString("0.################").PadLeft(24) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }
}