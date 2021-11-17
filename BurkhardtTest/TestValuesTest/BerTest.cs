using System;
using Burkardt.Values;

namespace TestValuesTest;

public static class BerTest
{
    public static void ber0_values_test ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BER0_VALUES_TEST tests BER0_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 June 2006
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
        Console.WriteLine("BER0_VALUES_TEST:");
        Console.WriteLine("  BER0_VALUES stores values of ");
        Console.WriteLine("  the Kelvin function BER of order 0.");
        Console.WriteLine("");
        Console.WriteLine("                X                     FX");
        Console.WriteLine("");
        n_data = 0;
        for ( ; ; )
        {
            ber.ber0_values ( ref n_data, ref x, ref fx );
            if ( n_data == 0 )
            {
                break;
            }
            Console.WriteLine("  " + x.ToString("0.################").PadLeft(24)
                                   + "  " + fx.ToString("0.################").PadLeft(24) + "");
        }
    }
    public static void ber1_values_test ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BER1_VALUES_TEST tests BER1_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 June 2006
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
        Console.WriteLine("BER1_VALUES_TEST:");
        Console.WriteLine("  BER1_VALUES stores values of ");
        Console.WriteLine("  the Kelvin function BER of order 1.");
        Console.WriteLine("");
        Console.WriteLine("                X                     FX");
        Console.WriteLine("");
        n_data = 0;
        for ( ; ; )
        {
            ber.ber1_values ( ref n_data, ref x, ref fx );
            if ( n_data == 0 )
            {
                break;
            }
            Console.WriteLine("  " + x.ToString("0.################").PadLeft(24)
                                   + "  " + fx.ToString("0.################").PadLeft(24) + "");
        }
    }
}