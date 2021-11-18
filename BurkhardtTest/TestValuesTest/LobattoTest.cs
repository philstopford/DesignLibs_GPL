using System;
using Burkardt.Values;

namespace TestValuesTest;

public static class LobattoTest
{
    public static void lobatto_polynomial_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOBATTO_POLYNOMIAL_VALUES_TEST tests LOBATTO_POLYNOMIAL_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 May 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        int n = 0;
        int n_data;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("LOBATTO_POLYNOMIAL_VALUES_TEST:");
        Console.WriteLine("  LOBATTO_POLYNOMIAL_VALUES stores values of ");
        Console.WriteLine("  the completed Lobatto polynomials.");
        Console.WriteLine("");
        Console.WriteLine("     N    X             Lo(N)(X)");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Lobatto.lobatto_polynomial_values(ref n_data, ref n, ref x, ref fx);
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

    public static void lobatto_polynomial_derivatives_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOBATTO_POLYNOMIAL_DERIVATIVES_TEST tests LOBATTO_POLYNOMIAL_DERIVATIVES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 November 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        int n = 0;
        int n_data;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("LOBATTO_POLYNOMIAL_DERIVATIVES_TEST:");
        Console.WriteLine("  LOBATTO_POLYNOMIAL_VALUES stores derivatives of ");
        Console.WriteLine("  the completed Lobatto polynomials.");
        Console.WriteLine("");
        Console.WriteLine("     N    X             Lo'(N)(X)");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Lobatto.lobatto_polynomial_derivatives(ref n_data, ref n, ref x, ref fx);
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