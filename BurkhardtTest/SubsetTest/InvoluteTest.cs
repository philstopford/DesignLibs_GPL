using System;
using Burkardt.Function;

namespace SubsetTestNS;

public static class InvoluteTest
{
    public static void involute_enum_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INVOLUTE_ENUM_TEST tests INVOLUTE_ENUM;
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
        const int N = 10;

        int i;
        int[] s = new int[N+1];

        Console.WriteLine("");
        Console.WriteLine("INVOLUTE_ENUM_TEST");
        Console.WriteLine("  INVOLUTE_ENUM counts involutions;");
        Console.WriteLine("");

        Involute.involute_enum ( N, ref s );

        Console.WriteLine("");
        Console.WriteLine("  N    # of involutions");
        Console.WriteLine("");

        for ( i = 0; i <= N; i++ )
        {
            Console.WriteLine(i.ToString().PadLeft(10)    + "  "
                                                          + s[i].ToString().PadLeft(10) + "");
        }

    }
}