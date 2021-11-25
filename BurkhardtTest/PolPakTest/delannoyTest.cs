using System;
using Burkardt.Sequence;
using Burkardt.Types;

namespace PolPakTest;

public static class delannoyTest
{
    public static void delannoy_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DELANNOY_TEST tests DELANNOY.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 June 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int m = 8;
        const int n = 8;

        Console.WriteLine("");
        Console.WriteLine("DELANNOY_TEST");
        Console.WriteLine("  DELANNOY computes the Delannoy numbers A(0:M,0:N).");
        Console.WriteLine("  A(M,N) counts the paths from (0,0) to (M,N).");
        Console.WriteLine("");

        int[] a = Delannoy.delannoy ( m, n );

        typeMethods.i4mat_print ( m + 1, n + 1, a, "  The Delannoy numbers:" );

    }

}