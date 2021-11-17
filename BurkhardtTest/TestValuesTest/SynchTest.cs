using System;
using Burkardt.Values;

namespace TestValuesTest;

public class SynchTest
{
    public static void synch1_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SYNCH1_VALUES_TEST tests SYNCH1_VALUES.
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
        double fx = 0;
        int n_data;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("SYNCH1_VALUES_TEST:");
        Console.WriteLine("  SYNCH1_VALUES stores values of ");
        Console.WriteLine("  the Synchrotron function of order 1.");
        Console.WriteLine("");
        Console.WriteLine("                X                     FX");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Synch.synch1_values(ref n_data, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + x.ToString("0.################").PadLeft(24) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void synch2_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SYNCH2_VALUES_TEST tests SYNCH2_VALUES.
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
        double fx = 0;
        int n_data;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("SYNCH2_VALUES_TEST:");
        Console.WriteLine("  SYNCH2_VALUES stores values of ");
        Console.WriteLine("  the Synchrotron function of order 2.");
        Console.WriteLine("");
        Console.WriteLine("                X                     FX");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Synch.synch2_values(ref n_data, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + x.ToString("0.################").PadLeft(24) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }
}