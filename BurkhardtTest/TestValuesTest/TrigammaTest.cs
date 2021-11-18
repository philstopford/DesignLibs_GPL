using System;
using Burkardt.Values;

namespace TestValuesTest;

public class TrigammaTest
{
    public static void trigamma_values_test()
        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    TRIGAMMA_VALUES_TEST tests TRIGAMMA_VALUES.
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
        Console.WriteLine("TRIGAMMA_VALUES_TEST");
        Console.WriteLine("  TRIGAMMA_VALUES stores values of");
        Console.WriteLine("  the TriGamma function.");
        Console.WriteLine("");
        Console.WriteLine("      X            FX");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Trigamma.trigamma_values(ref n_data, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

}