using System;
using Burkardt.Values;

namespace TestValuesTest;

public static class BernsteinTest
{
    public static void bernstein_poly_01_values_test ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BERNSTEIN_POLY_01_VALUES_TEST tests BERNSTEIN_POLY_01_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double b = 0;
        int k = 0;
        int n = 0;
        int n_data;
        double x = 0;
        Console.WriteLine("");
        Console.WriteLine("BERNSTEIN_POLY_01_VALUES_TEST:");
        Console.WriteLine("  BERNSTEIN_POLY_01_VALUES returns values of ");
        Console.WriteLine("  the Bernstein Polynomials.");
        Console.WriteLine("");
        Console.WriteLine("     N     K       X      BERNSTEIN(N,K)(X)");
        Console.WriteLine("");
        n_data = 0;
        for ( ; ; )
        {
            Bernstein.bernstein_poly_01_values ( ref n_data, ref n, ref k, ref x, ref b );
            if ( n_data == 0 )
            {
                break;
            }
            Console.WriteLine("  "
                              + n.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + k.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + b.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }
}