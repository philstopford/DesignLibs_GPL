using System;
using Burkardt.Values;

namespace TestValuesTest;

public static class AGMTest
{
    public static void agm_values_test ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    AGM_VALUES_TEST tests AGM_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 February 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a = 0;
        double b = 0;
        double fx = 0;

        Console.WriteLine("");
        Console.WriteLine("AGM_VALUES_TEST:");
        Console.WriteLine("  AGM_VALUES stores values of ");
        Console.WriteLine("  the arithmetic geometric mean function.");
        Console.WriteLine("");
        Console.WriteLine("           A          B              AGM(A,B)");
        Console.WriteLine("");
        int n_data = 0;
        for ( ; ; )
        {
            AGM.agm_values ( ref n_data, ref a, ref b, ref fx );
            if ( n_data == 0 )
            {
                break;
            }
            Console.WriteLine("  "
                              + a.ToString("0.######").PadLeft(14)  + "  "
                              + b.ToString("0.######").PadLeft(14)  + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }
}