using System;
using Burkardt.RankingNS;
using Burkardt.Types;

namespace ComboTest;

internal static partial class Program
{
    private static void ksubset_colex_check_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    KSUBSET_COLEX_CHECK_TEST tests KSUBSET_COLEX_CHECK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int k = 0;
        int n = 0;
        int[] s = new int[1];
        int[] s2 =  {
                5, 3, 2
            }
            ;
        int[] s3 =  {
                5, 2, 3
            }
            ;
        int[] s4 =  {
                7, 3, 2
            }
            ;
        int[] s5 =  {
                5, 3, 2
            }
            ;
        int test;

        Console.WriteLine("");
        Console.WriteLine("KSUBSET_COLEX_CHECK TEST");
        Console.WriteLine("  KSUBSET_COLEX_CHECK checks a K subset of an N set.");

        for (test = 1; test <= 7; test++)
        {
            switch (test)
            {
                case 1:
                    k = -1;
                    n = 5;
                    s = null;
                    break;
                case 2:
                    k = 3;
                    n = 0;
                    s = typeMethods.i4vec_copy_new(k, s2);
                    break;
                case 3:
                    k = 3;
                    n = 5;
                    s = typeMethods.i4vec_copy_new(k, s3);
                    break;
                case 4:
                    k = 3;
                    n = 5;
                    s = typeMethods.i4vec_copy_new(k, s4);
                    break;
                case 5:
                    k = 3;
                    n = 5;
                    s = typeMethods.i4vec_copy_new(k, s5);
                    break;
                case 6:
                    k = 0;
                    n = 5;
                    s = null;
                    break;
                case 7:
                    k = 0;
                    n = 0;
                    s = null;
                    break;
            }

            bool check = Ranking.ksubset_colex_check(k, n, s);
            typeMethods.i4vec_transpose_print(k, s, "  Subset:");
            Console.WriteLine("  N = " + n + ", K = " + k + "");
            Console.WriteLine("  Check = " + check + "");
        }
    }

    private static void ksubset_colex_rank_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    KSUBSET_COLEX_RANK_TEST tests KSUBSET_COLEX_RANK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] t =  {
                5, 3, 1
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("KSUBSET_COLEX_RANK_TEST");
        Console.WriteLine("  KSUBSET_COLEX_RANK ranks");
        Console.WriteLine("  K-subsets of an N set,");
        Console.WriteLine("  using the colexicographic ordering:");

        const int k = 3;
        const int n = 5;
        typeMethods.i4vec_transpose_print(k, t, "  Element to be ranked:");

        int rank = Ranking.ksubset_colex_rank(k, n, t);

        Console.WriteLine("");
        Console.WriteLine("  Rank is computed as " + rank + "");
    }

    private static void ksubset_colex_successor_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    KSUBSET_COLEX_SUCCESSOR_TEST tests KSUBSET_COLEX_SUCCESSOR.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("KSUBSET_COLEX_SUCCESSOR_TEST");
        Console.WriteLine("  KSUBSET_COLEX_SUCCESSOR lists");
        Console.WriteLine("  K-subsets of an N set using the colexicographic ordering:");

        const int k = 3;
        const int n = 5;
        int[] t = new int[k];

        int rank = -1;

        for (;;)
        {
            int rank_old = rank;

            Ranking.ksubset_colex_successor(k, n, ref t, ref rank);

            if (rank <= rank_old)
            {
                break;
            }

            string cout = "  " + rank.ToString().PadLeft(4);
            int i;
            for (i = 0; i < k; i++)
            {
                cout += "  " + t[i].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);
        }
    }

    private static void ksubset_colex_unrank_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    KSUBSET_COLEX_UNRANK_TEST tests KSUBSET_COLEX_UNRANK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("KSUBSET_COLEX_UNRANK_TEST");
        Console.WriteLine("  KSUBSET_COLEX_UNRANK unranks");
        Console.WriteLine("  K-subsets of an N set");
        Console.WriteLine("  using the colexicographic ordering.");

        const int rank = 5;
        const int k = 3;
        const int n = 5;

        int[] t = Ranking.ksubset_colex_unrank(rank, k, n);

        Console.WriteLine("");
        Console.WriteLine("  The element of rank " + rank + "");
        Console.WriteLine("");
        typeMethods.i4vec_transpose_print(k, t, "");
    }

    private static void ksubset_enum_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    KSUBSET_ENUM_TEST tests KSUBSET_ENUM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n;

        Console.WriteLine("");
        Console.WriteLine("KSUBSET_ENUM_TEST");
        Console.WriteLine("  KSUBSET_ENUM enumerates K subsets of an N set.");
        Console.WriteLine("");
        Console.WriteLine("      K:      0       1       2       3       4       5");
        Console.WriteLine("   N");

        for (n = 0; n <= 10; n++)
        {
            string cout = "  " + n.ToString().PadLeft(2)
                               + ":  ";
            int k;
            for (k = 0; k <= Math.Min(n, 5); k++)
            {
                cout += "  " + Ranking.ksubset_enum(k, n).ToString().PadLeft(6);
            }

            Console.WriteLine(cout);
        }
    }

    private static void ksubset_lex_check_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    KSUBSET_LEX_CHECK_TEST tests KSUBSET_LEX_CHECK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int k = 0;
        int n = 0;
        int[] s = new int[1];
        int[] s2 =  {
                2, 3, 5
            }
            ;
        int[] s3 =  {
                3, 2, 5
            }
            ;
        int[] s4 =  {
                2, 3, 7
            }
            ;
        int[] s5 =  {
                2, 3, 5
            }
            ;
        int test;

        Console.WriteLine("");
        Console.WriteLine("KSUBSET_LEX_CHECK TEST");
        Console.WriteLine("  KSUBSET_LEX_CHECK checks a K subset of an N set.");

        for (test = 1; test <= 7; test++)
        {
            switch (test)
            {
                case 1:
                    k = -1;
                    n = 5;
                    s = null;
                    break;
                case 2:
                    k = 3;
                    n = 0;
                    s = typeMethods.i4vec_copy_new(k, s2);
                    break;
                case 3:
                    k = 3;
                    n = 5;
                    s = typeMethods.i4vec_copy_new(k, s3);
                    break;
                case 4:
                    k = 3;
                    n = 5;
                    s = typeMethods.i4vec_copy_new(k, s4);
                    break;
                case 5:
                    k = 3;
                    n = 5;
                    s = typeMethods.i4vec_copy_new(k, s5);
                    break;
                case 6:
                    k = 0;
                    n = 5;
                    s = null;
                    break;
                case 7:
                    k = 0;
                    n = 0;
                    s = null;
                    break;
            }

            bool check = Ranking.ksubset_lex_check(k, n, s);
            typeMethods.i4vec_transpose_print(k, s, "  Subset:");
            Console.WriteLine("  N = " + n + ", K = " + k + "");
            Console.WriteLine("  Check = " + check + "");
        }
    }

    private static void ksubset_lex_rank_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    KSUBSET_LEX_RANK_TEST tests KSUBSET_LEX_RANK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] t =  {
                1, 4, 5
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("KSUBSET_LEX_RANK_TEST");
        Console.WriteLine("  KSUBSET_LEX_RANK ranks");
        Console.WriteLine("  K-subsets of an N set,");
        Console.WriteLine("  using the lexicographic ordering:");

        const int k = 3;
        const int n = 5;
        typeMethods.i4vec_transpose_print(k, t, "  Element to be ranked:");

        int rank = Ranking.ksubset_lex_rank(k, n, t);

        Console.WriteLine("");
        Console.WriteLine("  Rank is computed as " + rank + "");
    }

    private static void ksubset_lex_successor_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    KSUBSET_LEX_SUCCESSOR_TEST tests KSUBSET_LEX_SUCCESSOR.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("KSUBSET_LEX_SUCCESSOR_TEST");
        Console.WriteLine("  KSUBSET_LEX_SUCCESSOR lists");
        Console.WriteLine("  K-subsets of an N set using the lexicographic ordering:");

        const int k = 3;
        const int n = 5;
        int[] t = new int[k];

        int rank = -1;

        for (;;)
        {
            int rank_old = rank;

            Ranking.ksubset_lex_successor(k, n, ref t, ref rank);

            if (rank <= rank_old)
            {
                break;
            }

            string cout = "  " + rank.ToString().PadLeft(4);
            int i;
            for (i = 0; i < k; i++)
            {
                cout += "  " + t[i].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);
        }
    }

    private static void ksubset_lex_unrank_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    KSUBSET_LEX_UNRANK_TEST tests KSUBSET_LEX_UNRANK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("KSUBSET_LEX_UNRANK_TEST");
        Console.WriteLine("  KSUBSET_LEX_UNRANK unranks");
        Console.WriteLine("  K-subsets of an N set");
        Console.WriteLine("  using the lexicographic ordering.");

        const int rank = 5;
        const int k = 3;
        const int n = 5;

        int[] t = Ranking.ksubset_lex_unrank(rank, k, n);

        Console.WriteLine("");
        Console.WriteLine("  The element of rank " + rank + "");
        Console.WriteLine("");
        typeMethods.i4vec_transpose_print(k, t, "");
    }

    private static void ksubset_revdoor_rank_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    KSUBSET_REVDOOR_RANK_TEST tests KSUBSET_REVDOOR_RANK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] t =  {
                2, 4, 5
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("KSUBSET_REVDOOR_RANK_TEST");
        Console.WriteLine("  KSUBSET_REVDOOR_RANK ranks");
        Console.WriteLine("  K-subsets of an N set,");
        Console.WriteLine("  using the revolving door ordering:");

        const int k = 3;
        const int n = 5;
        typeMethods.i4vec_transpose_print(k, t, "  Element to be ranked:");

        int rank = Ranking.ksubset_revdoor_rank(k, n, t);

        Console.WriteLine("");
        Console.WriteLine("  Rank is computed as " + rank + "");
    }

    private static void ksubset_revdoor_successor_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    KSUBSET_REVDOOR_SUCCESSOR_TEST tests KSUBSET_REVDOOR_SUCCESSOR.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("KSUBSET_REVDOOR_SUCCESSOR_TEST");
        Console.WriteLine("  KSUBSET_REVDOOR_SUCCESSOR lists");
        Console.WriteLine("  K-subsets of an N set using the revolving door ordering:");

        const int k = 3;
        const int n = 5;
        int[] t = new int[k];

        int rank = -1;

        for (;;)
        {
            int rank_old = rank;

            Ranking.ksubset_revdoor_successor(k, n, ref t, ref rank);

            if (rank <= rank_old)
            {
                break;
            }

            string cout = "  " + rank.ToString().PadLeft(4);
            int i;
            for (i = 0; i < k; i++)
            {
                cout += "  " + t[i].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);
        }
    }

    private static void ksubset_revdoor_unrank_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    KSUBSET_REVDOOR_UNRANK_TEST tests KSUBSET_REVDOOR_UNRANK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("KSUBSET_REVDOOR_UNRANK_TEST");
        Console.WriteLine("  KSUBSET_REVDOOR_UNRANK unranks");
        Console.WriteLine("  K-subsets of an N set");
        Console.WriteLine("  using the revolving door ordering.");

        const int rank = 5;
        const int k = 3;
        const int n = 5;

        int[] t = Ranking.ksubset_revdoor_unrank(rank, k, n);

        Console.WriteLine("");
        Console.WriteLine("  The element of rank " + rank + "");
        Console.WriteLine("");
        typeMethods.i4vec_transpose_print(k, t, "");
    }
}