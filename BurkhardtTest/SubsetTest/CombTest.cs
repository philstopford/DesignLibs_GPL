﻿using System;
using Burkardt.SubsetNS;
using Burkardt.Types;

namespace SubsetTestNS;

public static class CombTest
{
    public static void comb_next_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COMB_NEXT_TEST tests COMB_NEXT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 April 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 5;

        int[] a = new int[N];
        int k;

        Console.WriteLine("");
        Console.WriteLine("COMB_NEXT_TEST");
        Console.WriteLine("  COMB_NEXT produces combinations.");

        for (k = 1; k <= N; k++)
        {
            Console.WriteLine("");
            Console.WriteLine("  Combinations of size K = " + k + "");
            Console.WriteLine("");

            bool done = true;

            for (;;)
            {
                Comb.comb_next(N, k, a, ref done);

                if (done)
                {
                    break;
                }

                string cout = "";
                int i;
                for (i = 0; i < k; i++)
                {
                    cout += "  " + a[i].ToString().PadLeft(4);
                }

                Console.WriteLine(cout);
            }
        }

    }

    public static void comb_row_next_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COMB_ROW_NEXT_TEST tests COMB_ROW_TEST.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 December 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n;
        const int N_MAX = 10;

        Console.WriteLine("");
        Console.WriteLine("COMB_ROW_NEXT_TEST");
        Console.WriteLine("  COMB_ROW_NEXT computes the next row of the Pascal triangle.");
        Console.WriteLine("");

        int[] c = new int[N_MAX + 1];

        for (n = 0; n <= N_MAX; n++)
        {
            Comb.comb_row_next(n, ref c);
            string cout = "  " + n.ToString().PadLeft(2) + "  ";
            int i;
            for (i = 0; i <= n; i++)
            {
                cout += c[i].ToString().PadLeft(5);
            }

            Console.WriteLine(cout);
        }
    }

    public static void comb_unrank_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COMB_UNRANK_TEST tests COMB_UNRANK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 October 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 5;

        int[] a = new int[N];
        int i;
        const int m = 10;
        int rank;
        string cout;

        int cnk = typeMethods.i4_choose(m, N);

        Console.WriteLine("");
        Console.WriteLine("COMB_UNRANK_TEST");
        Console.WriteLine("  COMB_UNRANK returns a combination of N things");
        Console.WriteLine("  out of M, given the lexicographic rank.");
        Console.WriteLine("");
        Console.WriteLine("  The total set size is M = " + m + "");
        Console.WriteLine("  The subset size is N =    " + N + "");
        Console.WriteLine("  The number of combinations of N out of M is " + cnk + "");
        Console.WriteLine("");
        Console.WriteLine("   Rank	  Combination");
        Console.WriteLine("");

        for (rank = 1; rank <= 3; rank++)
        {
            Comb.comb_unrank(m, N, rank, ref a);
            cout = "  "
                   + rank.ToString().PadLeft(3) + "  ";
            for (i = 0; i < N; i++)
            {
                cout += a[i].ToString().PadLeft(4) + "  ";
            }

            Console.WriteLine(cout);
        }

        for (rank = 6; rank <= 8; rank++)
        {
            Comb.comb_unrank(m, N, rank, ref a);
            cout = rank.ToString().PadLeft(3) + "  ";
            for (i = 0; i < N; i++)
            {
                cout += a[i].ToString().PadLeft(4) + "  ";
            }

            Console.WriteLine(cout);
        }

        for (rank = 250; rank <= 252; rank++)
        {
            Comb.comb_unrank(m, N, rank, ref a);
            cout = "  " + rank.ToString().PadLeft(3) + "  ";
            for (i = 0; i < N; i++)
            {
                cout += a[i].ToString().PadLeft(4) + "  ";
            }

            Console.WriteLine(cout);
        }
    }

}