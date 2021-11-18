using System;
using Burkardt.IntegralNS;

namespace PolPakTest;

public static class sinpowerTest
{
    public static void sin_power_int_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SIN_POWER_INT_TEST tests SIN_POWER_INT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 June 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a = 0;
        double b = 0;
        double fx = 0;
        double fx2 = 0;
        int n = 0;
        int n_data;

        Console.WriteLine("");
        Console.WriteLine("SIN_POWER_INT_TEST:");
        Console.WriteLine("  SIN_POWER_INT computes the integral of the N-th power");
        Console.WriteLine("  of the sine function.");
        Console.WriteLine("");
        Console.WriteLine("         A         B       N        Exact    Computed");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            Burkardt.Values.Sine.sin_power_int_values(ref n_data, ref a, ref b, ref n, ref fx);

            if (n_data == 0)
            {
                break;
            }

            fx2 = SinPower.sin_power_int(a, b, n);

            Console.WriteLine("  "
                              + a.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + b.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + n.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + fx.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + fx2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

}