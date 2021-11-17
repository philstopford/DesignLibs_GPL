using System;
using Burkardt.Values;

namespace TestValuesTest;

public static class TangentTest
{
    public static void tan_values_test ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TAN_VALUES_TEST tests TAN_VALUES.
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
        Console.WriteLine("TAN_VALUES_TEST:");
        Console.WriteLine("   TAN_VALUES stores values of the tangent function.");
        Console.WriteLine("");
        Console.WriteLine("                X                     FX");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Tangent.tan_values(ref n_data, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + x.ToString("0.################").PadLeft(24) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }
    public static void tanh_values_test ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TANH_VALUES_TEST tests TANH_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 June 2007
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
        Console.WriteLine("TANH_VALUES_TEST:");
        Console.WriteLine("   TANH_VALUES stores values of the hyperbolic tangent function.");
        Console.WriteLine("");
        Console.WriteLine("                X                     FX");
        Console.WriteLine("");
        n_data = 0;
        for ( ; ; )
        {
            Tangent.tanh_values ( ref n_data, ref x, ref fx );
            if ( n_data == 0 )
            {
                break;
            }
            Console.WriteLine("  "
                              + x.ToString("0.################").PadLeft(24) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }
}