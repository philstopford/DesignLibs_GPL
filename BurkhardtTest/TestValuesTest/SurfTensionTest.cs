using System;
using System.Globalization;
using Burkardt.Values;

namespace TestValuesTest;

public static class SurfTensionTest
{
    public static void surten_values_test ( )
        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    SURTEN_VALUES_TEST tests SURTEN_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double sigma = 0;
        double tc = 0;
        Console.WriteLine("");
        Console.WriteLine("SURTEN_VALUES_TEST:");
        Console.WriteLine("  SURTEN_VALUES stores values of");
        Console.WriteLine("  the surface tension of water");
        Console.WriteLine("  as a function of temperature.");
        Console.WriteLine("");
        Console.WriteLine("      T            SIGMA(T)");
        Console.WriteLine("");
        int n_data = 0;
        for ( ; ; )
        {
            SurfTension.surten_values ( ref n_data, ref tc, ref sigma );
            if ( n_data == 0 )
            {
                break;
            }
            Console.WriteLine( "  "
                               + tc.ToString(CultureInfo.InvariantCulture).PadLeft(12) + tc    + "  "
                               + sigma.ToString(CultureInfo.InvariantCulture).PadLeft(12) + sigma + "");
        }
    }
}