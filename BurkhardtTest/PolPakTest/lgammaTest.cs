using System;
using Burkardt;

namespace PolPakTest;

public static class lgammaTest
{
    public static void lgamma_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LGAMMA_TEST tests LGAMMA.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        double fx2 = 0;
        int n_data = 0;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("LGAMMA_TEST:");
        Console.WriteLine("  LGAMMA is a C math library function which evaluates");
        Console.WriteLine("  the logarithm of the Gamma function.");
        Console.WriteLine("");
        Console.WriteLine("     X       Exact F       LGAMMA(X)");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            Burkardt.Values.Gamma.gamma_log_values(ref n_data, ref x, ref fx);

            if (n_data == 0)
            {
                break;
            }

            fx2 = Helpers.LogGamma(x);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + fx.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                              + fx2.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");

        }

    }

}