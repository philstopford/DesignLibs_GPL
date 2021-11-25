using System;
using System.Globalization;
using Burkardt.Values;

namespace TestValuesTest;

public static class LaguerreTest
{
    public static void laguerre_associated_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGUERRE_ASSOCIATED_VALUES_TEST tests LAGUERRE_ASSOCIATED_VALUES.
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
        int m = 0;
        int n = 0;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("LAGUERRE_ASSOCIATED_VALUES_TEST:");
        Console.WriteLine("  LAGUERRE_ASSOCIATED_VALUES stores values of");
        Console.WriteLine("  the associated Laguerre polynomials.");
        Console.WriteLine("");
        Console.WriteLine("     N     M    X             L(N,M)(X)");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Laguerre.laguerre_associated_values(ref n_data, ref n, ref m, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + n.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + m.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void laguerre_general_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGUERRE_GENERAL_VALUES_TEST tests LAGUERRE_GENERAL_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    26 January 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a = 0;
        double fx = 0;
        int n = 0;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("LAGUERRE_GENERAL_VALUES_TEST:");
        Console.WriteLine("  LAGUERRE_GENERAL_VALUES stores values of");
        Console.WriteLine("  the generalized Laguerre function.");
        Console.WriteLine("");
        Console.WriteLine("     N     A    X             L(N,A)(X)");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Laguerre.laguerre_general_values(ref n_data, ref n, ref a, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + n.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + a.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

    public static void laguerre_polynomial_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGUERRE_POLYNOMIAL_VALUES_TEST tests LAGUERRE_POLYNOMIAL_VALUES.
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
        int n = 0;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("LAGUERRE_POLYNOMIAL_VALUES_TEST:");
        Console.WriteLine("  LAGUERRE_POLYNOMIAL_VALUES stores values of ");
        Console.WriteLine("  the Laguerre polynomials.");
        Console.WriteLine("");
        Console.WriteLine("     N     X            L(N)(X)");
        Console.WriteLine("");
        int n_data = 0;
        for (;;)
        {
            Laguerre.laguerre_polynomial_values(ref n_data, ref n, ref x, ref fx);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + n.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + fx.ToString("0.################").PadLeft(24) + "");
        }
    }

}