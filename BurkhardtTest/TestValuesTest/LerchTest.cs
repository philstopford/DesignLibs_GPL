using System;
using Burkardt.Values;

namespace TestValuesTest;

public class LerchTest
{
    public static void lerch_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LERCH_VALUES_TEST tests LERCH_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 June 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a = 0;
        double fx = 0;
        int n_data;
        int s = 0;
        double z = 0;
        Console.WriteLine("");
        Console.WriteLine("LERCH_VALUES_TEST:");
        Console.WriteLine("  LERCH_VALUES returns values of");
        Console.WriteLine("  the Lerch transcendent function.");
        Console.WriteLine("");
        Console.WriteLine("      Z      S      A      Fx");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Lerch.lerch_values(ref n_data, ref z, ref s, ref a, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + z.ToString("0.################").PadLeft(24) + "  "
                              + s.ToString().PadLeft(6) + "  "
                              + a.ToString().PadLeft(12) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }
}