using System;
using Burkardt.RankingNS;
using Burkardt.Types;

namespace ComboTest;

internal static partial class Program
{
    private static void queens_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    QUEENS_TEST tests QUEENS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 July 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 8;

        int[] iarray = new int[N];
        int[] istack = new int[N*N];
        int k = 0;
        const int maxstack = N * N;
        int nstack = 0;

        Console.WriteLine("");
        Console.WriteLine("QUEENS_TEST");
        Console.WriteLine("  QUEENS produces nonattacking queens");
        Console.WriteLine("  on a chessboard using a backtrack search.");
        Console.WriteLine("");

        int indx = 0;

        for ( ; ; )
        {
            Ranking.backtrack ( N, ref iarray, ref indx, ref k, ref nstack, ref istack, maxstack );

            if ( indx == 1 )
            {
                typeMethods.i4vec_transpose_print ( N, iarray, "" );
            }
            else if ( indx == 2 )
            {
                Ranking.queens ( N, iarray, k, ref nstack, ref istack, maxstack );
            }
            else
            {
                break;
            }
        }
    }
}