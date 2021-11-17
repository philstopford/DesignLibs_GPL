using System;
using Burkardt.RankingNS;
using Burkardt.Types;

namespace ComboTest;

internal partial class Program
{
    private static void bell_numbers_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BELL_NUMBERS_TEST tests BELL_NUMBERS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 July 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] b;
        int bn = 0;
        int n = 0;
        int n_data;

        n_data = 0;

        Console.WriteLine("");
        Console.WriteLine("BELL_NUMBERS_TEST");
        Console.WriteLine("  BELL_NUMBERS computes Bell numbers.");
        Console.WriteLine("");
        Console.WriteLine("     N          BELL(N)      BELL_NUMBERS(N)");
        Console.WriteLine("");
        for ( ; ; )
        {
            Ranking.bell_values ( ref n_data, ref n, ref bn );

            if ( n_data == 0 )
            {
                break;
            }
            b = Ranking.bell_numbers ( n );
            Console.WriteLine("  " + n.ToString().PadLeft(8)
                                   + "  " + bn.ToString().PadLeft(12)
                                   + "  " + b[n].ToString().PadLeft(12) + "");
        }
    }
}