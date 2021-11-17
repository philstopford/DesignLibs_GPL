using System;
using Burkardt.Values;

namespace TestValuesTest;

public static class DilogarithmTest
{
    public static void dilogarithm_values_test ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DILOGARITHM_VALUES_TEST tests DILOGARITHM_VALUES.
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
        int n_data;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("DILOGARITHM_VALUES_TEST:");
        Console.WriteLine("  DILOGARITHM_VALUES stores values of");
        Console.WriteLine("  the dilogarithm function.");
        Console.WriteLine("");
        Console.WriteLine("      X          DILOGARITHM(X)");
        Console.WriteLine("");
        n_data = 0;
        for ( ; ; )
        {
            Dilogarithm.dilogarithm_values ( ref n_data, ref x, ref fx );
            if ( n_data == 0 )
            {
                break;
            }
            Console.WriteLine("  "
                              + x.ToString().PadLeft(12)  + "  "
                              + fx.ToString().PadLeft(12) + "");
        }
    }
}