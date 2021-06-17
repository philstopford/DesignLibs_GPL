using System;
using Burkardt.RankingNS;
using Burkardt.Types;

namespace Burkardt.ComboTest
{
    partial class Program
    {
        static void marriage_test ( )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MARRIAGE_TEST tests MARRIAGE.
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
            int N = 5;

            int[] fiancee;
            int i;
            int n = N;
            int[] next;
            int[] prefer = {
                2, 1, 2, 1, 5, 
                5, 2, 3, 3, 3, 
                1, 3, 5, 2, 2, 
                3, 4, 4, 4, 1, 
                4, 5, 1, 5, 4 };
            int[] rank = {
                2, 4, 1, 4, 5, 
                4, 3, 3, 2, 2, 
                5, 5, 4, 1, 3, 
                3, 1, 2, 3, 1, 
                1, 2, 5, 5, 4 };

            Console.WriteLine("");
            Console.WriteLine("MARRIAGE_TEST");
            Console.WriteLine("  MARRIAGE arranges a set of stable marriages");
            Console.WriteLine("  given a set of preferences.");

            fiancee = new int[n];
            next = new int[n];

            Ranking.marriage ( n, prefer, rank, ref fiancee, ref next );

            Console.WriteLine("");
            Console.WriteLine("  Man, Wife's rank, Wife");
            Console.WriteLine("");
            for ( i = 1; i <= n; i++ )
            {
                Console.WriteLine("  " + i.ToString().PadLeft(4)
                    + "  " + next[i-1].ToString().PadLeft(4)
                    + "  " + prefer[i-1+(next[i-1]-1)*n].ToString().PadLeft(4) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("  Woman, Husband's rank, Husband");
            Console.WriteLine("");
            for ( i = 1; i <= n; i++ )
            {
                Console.WriteLine("  " + i.ToString().PadLeft(4)
                    + "  " + rank[i-1+(fiancee[i-1]-1)*n].ToString().PadLeft(4)
                    + "  " + fiancee[i-1].ToString().PadLeft(4) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("  Correct result:");
            Console.WriteLine("");
            Console.WriteLine("  M:W 1  2  3  4  5");
            Console.WriteLine("   1  +  .  .  .  .");
            Console.WriteLine("   2  .  .  .  +  .");
            Console.WriteLine("   3  .  .  .  .  +");
            Console.WriteLine("   4  .  .  +  .  .");
            Console.WriteLine("   5  .  +  .  .  .");
        }
    }
}