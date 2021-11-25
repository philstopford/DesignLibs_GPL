using System;
using System.Globalization;
using Burkardt.Values;

namespace TestValuesTest;

public static class DawsonTest
{
    public static void dawson_values_test ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DAWSON_VALUES_TEST tests DAWSON_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("DAWSON_VALUES_TEST:");
        Console.WriteLine("  DAWSON_VALUES stores values of");
        Console.WriteLine("  Dawson's integral function.");
        Console.WriteLine("");
        Console.WriteLine("      X          DAWSON(X)");
        Console.WriteLine("");
        int n_data = 0;
        for ( ; ; )
        {
            Dawson.dawson_values ( ref n_data, ref x, ref fx );
            if ( n_data == 0 )
            {
                break;
            }
            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12)  + "  "
                              + fx.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }
}