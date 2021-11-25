using System;
using Burkardt.Values;

namespace TestValuesTest;

public static class CotTest
{
    public static void cot_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COT_VALUES_TEST tests COT_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 January 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("COT_VALUES_TEST:");
        Console.WriteLine("   COT_VALUES stores values of the cotangent function.");
        Console.WriteLine("");
        Console.WriteLine("                X                     FX");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Cot.cot_values(ref n_data, ref x, ref fx);
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