using System;
using Burkardt.Values;

namespace TestValuesTest;

public static class GudermannianTest
{
    public static void gud_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GUD_VALUES_TEST tests GUD_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 January 2007
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
        Console.WriteLine("GUD_VALUES_TEST:");
        Console.WriteLine("  GUD_VALUES stores values of");
        Console.WriteLine("  the Gudermannian function.");
        Console.WriteLine("");
        Console.WriteLine("      X            GUD(X)");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Gudermannian.gud_values(ref n_data, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }
}