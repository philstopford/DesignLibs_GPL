using System;
using Burkardt.RankingNS;
using Burkardt.Types;

namespace ComboTest;

internal partial class Program
{
    private static void stirling_numbers1_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    STIRLING_NUMBERS1_TEST tests STIRLING_NUMBERS1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 July 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int maxm = 6;
        int maxn = 6;
        int[] s;

        Console.WriteLine("");
        Console.WriteLine("STIRLING_NUMBERS1_TEST");
        Console.WriteLine("  STIRLING_NUMBERS1 computes a table of Stirling");
        Console.WriteLine("  numbers of the first kind.");

        s = Ranking.stirling_numbers1 ( maxm, maxn );

        typeMethods.i4mat_print ( maxm + 1, maxn + 1, s, "  Stirling number of first kind" ); 
    }

    private static void stirling_numbers2_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    STIRLING_NUMBERS2_TEST tests STIRLING_NUMBERS2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 July 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int maxm = 6;
        int maxn = 6;
        int[] s;

        Console.WriteLine("");
        Console.WriteLine("STIRLING_NUMBERS2_TEST");
        Console.WriteLine("  STIRLING_NUMBERS2 computes a table of Stirling");
        Console.WriteLine("  numbers of the second kind.");

        s = Ranking.stirling_numbers2 ( maxm, maxn );

        typeMethods.i4mat_print ( maxm + 1, maxn + 1, s, "  Stirling number of second kind" ); 
    }
}