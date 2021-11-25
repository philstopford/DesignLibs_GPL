using System;
using Burkardt.Values;

namespace TestValuesTest;

public static class ClausenTest
{
    public static void clausen_values_test ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CLAUSEN_VALUES_TEST tests CLAUSEN_VALUES.
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
        Console.WriteLine("CLAUSEN_VALUES_TEST:");
        Console.WriteLine("  CLAUSEN_VALUES stores values of ");
        Console.WriteLine("  Clausen's integral function.");
        Console.WriteLine("");
        Console.WriteLine("                X                     FX");
        Console.WriteLine("");
        int n_data = 0;
        for ( ; ; )
        {
            Clausen.clausen_values ( ref n_data, ref x, ref fx );
            if ( n_data == 0 )
            {
                break;
            }
            Console.WriteLine("  "
                              + x.ToString("0.################").PadLeft(24)  + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }
}