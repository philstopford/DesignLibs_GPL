using System;
using Burkardt.RankingNS;
using Burkardt.Types;

namespace Burkardt.ComboTest
{
    partial class Program
    {
        static void backtrack_test ( )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BACKTRACK_TEST tests BACKTRACK.
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
            int N = 8;

            int[] iarray = new int[N];
            int indx;
            int[] istack = new int[N*N];
            int k = 0;
            int n = N;
            int maxstack = N * N;
            int nstack = 0;

            Console.WriteLine("");
            Console.WriteLine("BACKTRACK_TEST");
            Console.WriteLine("  BACKTRACK supervises a backtrack search.");
            Console.WriteLine("  Here we arrange non-attacking chess queens.");
            Console.WriteLine("");

            indx = 0;

            for ( ; ; )
            {
                Ranking.backtrack ( n, ref iarray, ref indx, ref k, ref nstack, ref istack, maxstack );

                if ( indx == 1 )
                {
                    typeMethods.i4vec_transpose_print ( n, iarray, "" );
                }
                else if ( indx == 2 )
                {
                    Ranking.queens ( n, iarray, k, ref nstack, ref istack, maxstack );
                }
                else
                {
                    break;
                }
            }
        }
    }
}