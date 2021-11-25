using System;
using Burkardt.Values;

namespace TestValuesTest;

public static class GoodwinTest
{
    public static void goodwin_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GOODWIN_VALUES_TEST tests GOODWIN_VALUES.
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
        double fx = 0;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("GOODWIN_VALUES_TEST:");
        Console.WriteLine("  GOODWIN_VALUES stores values of ");
        Console.WriteLine("  the Goodwin function.");
        Console.WriteLine("");
        Console.WriteLine("                X                     FX");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Goodwin.goodwin_values(ref n_data, ref x, ref fx);
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