using System;
using System.Globalization;
using Burkardt.Values;

namespace TestValuesTest;

public static class ThermCondTest
{
    public static void thercon_values_test()
        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    THERCON_VALUES_TEST tests THERCON_VALUES.
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
        double lambda = 0;
        double p = 0;
        double tc = 0;
        Console.WriteLine("");
        Console.WriteLine("THERCON_VALUES_TEST:");
        Console.WriteLine("  THERCON_VALUES stores values of");
        Console.WriteLine("  the thermal conductivity of water");
        Console.WriteLine("  as a function of temperature and pressure.");
        Console.WriteLine("");
        Console.WriteLine("      T            P            LAMBDA(T,P)");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            ThermCond.thercon_values(ref n_data, ref tc, ref p, ref lambda);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + tc.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + p.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + lambda.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

}