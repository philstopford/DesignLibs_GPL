using System;
using Burkardt;

namespace SubsetTestNS;

public static class InverseTest
{
    public static void inverse_mod_n_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INVERSE_MOD_N_TEST tests INVERSE_MOD_N.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 November 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int b;
        int n;
        int y;
        int z;

        Console.WriteLine("");
        Console.WriteLine("INVERSE_MOD_N_TEST");
        Console.WriteLine("  INVERSE_MOD_N seeks Y, the inverse of B mod N,");
        Console.WriteLine("  so that mod ( B * Y, N ) = 1, but returns 0");
        Console.WriteLine("  if the inverse does not exist.");

        Console.WriteLine("");
        Console.WriteLine("     B     N     Y     Z = ( ( B * Y ) % N )");
        Console.WriteLine("");

        for ( n = 1; n <= 10;  n++ )
        {
            for ( b = 1; b < n; b++ )
            {
                y = Helpers.inverse_mod_n ( b, n );
                z = b * y % n;
                Console.WriteLine("  " + b.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                       + "  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                       + "  " + y.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                       + "  " + z.ToString(CultureInfo.InvariantCulture).PadLeft(2) + "");
            }
        }
    }
}