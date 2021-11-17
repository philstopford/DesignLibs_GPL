using System;
using Burkardt.Values;

namespace TestValuesTest;

public static class ErrorFuncTest
{
    public static void erf_values_test ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ERF_VALUES_TEST tests ERF_VALUES.
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
        Console.WriteLine("ERF_VALUES_TEST:");
        Console.WriteLine("  ERF_VALUES stores values of");
        Console.WriteLine("  the error function ERF(X).");
        Console.WriteLine("");
        Console.WriteLine("      X          ERF(X)");
        Console.WriteLine("");
        n_data = 0;
        for ( ; ; )
        {
            ErrorFunc.erf_values ( ref n_data, ref x, ref fx );
            if ( n_data == 0 )
            {
                break;
            }
            Console.WriteLine("  "
                              + x.ToString().PadLeft(12) + "  "
                              + fx.ToString().PadLeft(12) + "");
        }
    }
    public static void erfc_values_test ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ERFC_VALUES_TEST tests ERFC_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 May 2007
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
        Console.WriteLine("ERFC_VALUES_TEST:");
        Console.WriteLine("  ERFC_VALUES stores values of");
        Console.WriteLine("  the complementary error function ERFC(X).");
        Console.WriteLine("");
        Console.WriteLine("      X          ERFC(X)");
        Console.WriteLine("");
        n_data = 0;
        for ( ; ; )
        {
            ErrorFunc.erfc_values ( ref n_data, ref x, ref fx );
            if ( n_data == 0 )
            {
                break;
            }
            Console.WriteLine("  "
                              + x.ToString().PadLeft(12) + "  "
                              + fx.ToString().PadLeft(12) + "");
        }
    }

}