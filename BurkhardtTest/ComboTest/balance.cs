using System;
using Burkardt.RankingNS;
using Burkardt.Types;

namespace ComboTest;

internal static partial class Program
{
    private static void bal_seq_check_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BAL_SEQ_CHECK_TEST tests BAL_SEQ_CHECK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int test;
        int[] t = new int[1];
        int[] t1 =  {
                0, 0, 1, 0, 1, 0, 0, 1, 1, 1
            }
            ;
        int[] t2 =  {
                1, 1, 0, 1, 0, 1, 1, 0, 0, 0
            }
            ;
        int[] t3 =  {
                0, 0, 1, 0, 1, 0, 0, 1, 0, 1
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("BAL_SEQ_CHECK TEST");
        Console.WriteLine("  BAL_SEQ_CHECK checks N and T(1:2*N).");
        Console.WriteLine("");
        Console.WriteLine("  Check?   N    T(1:2*N)");
        Console.WriteLine("");

        for (test = 1; test <= 3; test++)
        {
            const int n = 5;

            t = test switch
            {
                1 => typeMethods.i4vec_copy_new(2 * n, t1),
                2 => typeMethods.i4vec_copy_new(2 * n, t2),
                3 => typeMethods.i4vec_copy_new(2 * n, t3),
                _ => t
            };

            bool check = Ranking.bal_seq_check(n, t);
            string cout = "    "
                          + "  " + check.ToString().PadLeft(1)
                          + "  " + n.ToString().PadLeft(2);
            int i;
            for (i = 0; i < 2 * n; i++)
            {
                cout += "  " + t[i].ToString().PadLeft(2);
            }

            Console.WriteLine(cout);
        }
    }

    private static void bal_seq_enum_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BAL_SEQ_ENUM_TEST tests BAL_SEQ_ENUM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n;

        Console.WriteLine("");
        Console.WriteLine("BAL_SEQ_ENUM_TEST");
        Console.WriteLine("  BAL_SEQ_ENUM enumerates balanced sequences of N terms.");

        for (n = 0; n <= 10; n++)
        {
            int bal_seq_num = Ranking.bal_seq_enum(n);
            Console.WriteLine("  " + n.ToString().PadLeft(2)
                                   + "  " + bal_seq_num.ToString().PadLeft(6) + "");
        }
    }

    private static void bal_seq_rank_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BAL_SEQ_RANK_TEST tests BAL_SEQ_RANK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] t =  {
                0, 0, 1, 0, 1, 1, 0, 0, 1, 1
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("BAL_SEQ_RANK_TEST");
        Console.WriteLine("  BAL_SEQ_RANK ranks a balanced sequence of N items.");

        const int n = 5;
        int rank = Ranking.bal_seq_rank(n, t);

        typeMethods.i4vec_transpose_print(2 * n, t, "  Element to be ranked:");
        Console.WriteLine("");
        Console.WriteLine("  Rank is computed as: " + rank + "");
    }

    private static void bal_seq_successor_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BAL_SEQ_SUCCESSOR_TEST tests BAL_SEQ_SUCCESSOR.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("BAL_SEQ_SUCCESSOR_TEST:");
        Console.WriteLine("  BAL_SEQ_SUCCESSOR lists balanced sequences of N items, one at a time.");

        const int n = 5;
        int[] t = new int[2 * n];

        int rank = -1;

        for (;;)
        {
            int rank_old = rank;

            Ranking.bal_seq_successor(n, ref t, ref rank);

            if (rank <= rank_old)
            {
                break;
            }

            string cout = "  " + rank.ToString().PadLeft(4);
            int i;
            for (i = 0; i < 2 * n; i++)
            {
                cout += "  " + t[i].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);
        }
    }

    private static void bal_seq_to_tableau_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BAL_SEQ_TO_TABLEAU_TEST tests BAL_SEQ_TO_TABLEAU.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int n = 4;

        Console.WriteLine("");
        Console.WriteLine("BAL_SEQ_TO_TABLEAU_TEST");
        Console.WriteLine("  BAL_SEQ_TO_TABLEAU converts a balanced");
        Console.WriteLine("  sequence to a tableau;");
        //
        //  Pick a random balanced sequence.
        //
        const int rank = 7;

        int[] t = Ranking.bal_seq_unrank(rank, n);

        Console.WriteLine("");
        Console.WriteLine("  Random balanced sequence:");
        Console.WriteLine("");
        typeMethods.i4vec_transpose_print(2 * n, t, "");
        //
        //  Convert to a tableau.
        //
        int[] tab = Ranking.bal_seq_to_tableau(n, t);

        typeMethods.i4mat_print(2, n, tab, "  Corresponding tableau");
    }

    private static void bal_seq_unrank_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BAL_SEQ_UNRANK_TEST tests BAL_SEQ_UNRANK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("BAL_SEQ_UNRANK_TEST:");
        Console.WriteLine("  BAL_SEQ_UNRANK unranks a balanced sequence of N items.");

        const int rank = 21;
        const int n = 5;

        int[] t = Ranking.bal_seq_unrank(rank, n);
        Console.WriteLine("");
        Console.WriteLine("  The element of rank " + rank + "");
        Console.WriteLine("");
        typeMethods.i4vec_transpose_print(2 * n, t, "");
    }
}