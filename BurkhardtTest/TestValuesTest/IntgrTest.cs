using System;
using System.Globalization;
using Burkardt.Values;

namespace TestValuesTest;

public static class IntgrTest
{
    public static void int_values_test ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INT_VALUES_TEST tests INT_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 January 2015
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
        Console.WriteLine("INT_VALUES_TEST:");
        Console.WriteLine("  INT_VALUES stores values of the integer part of a real number.");
        Console.WriteLine("");
        Console.WriteLine("      X            INT(X)");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Intgr.int_values(ref n_data, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + fx.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }
}