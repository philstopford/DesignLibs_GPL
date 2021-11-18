using System;
using Burkardt.RankingNS;
using Burkardt.SubsetNS;

namespace SubsetTestNS;

public static class SubsetTest
{
    public static void subset_by_size_next_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SUBSET_BY_SIZE_NEXT_TEST tests SUBSET_BY_SIZE_NEXT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 June 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] a;
        int i;
        int m;
        int m2;
        bool more;
        bool more2;
        int n;
        int rank;
        int subsize;

        Console.WriteLine("");
        Console.WriteLine("SUBSET_BY_SIZE_NEXT_TEST");
        Console.WriteLine("  SUBSET_BY_SIZE_NEXT generates all subsets of an N set.");
        Console.WriteLine("");

        n = 5;
        a = new int[n];
        subsize = 0;
        more = false;
        more2 = false;
        m = 0;
        m2 = 0;

        rank = 0;

        for (;;)
        {
            Subset.subset_by_size_next(n, ref a, ref subsize, ref more, ref more2, ref m, ref m2);

            rank += 1;

            string cout = rank.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  ";

            switch (subsize)
            {
                case > 0:
                {
                    for (i = 0; i < subsize; i++)
                    {
                        cout += a[i].ToString(CultureInfo.InvariantCulture).PadLeft(2) + "  ";
                    }

                    Console.WriteLine(cout);
                    break;
                }
                default:
                    Console.WriteLine("The empty set");
                    break;
            }

            if (!more)
            {
                break;
            }

        }
    }

    public static void subset_lex_next_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SUBSET_LEX_NEXT_TEST tests SUBSET_LEX_NEXT with size restrictions.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int NDIM = 3;

        int[] a = new int[NDIM];
        int i;
        int k;
        bool ltest;
        int n = 5;

        Console.WriteLine("");
        Console.WriteLine("SUBSET_LEX_NEXT_TEST");
        Console.WriteLine("  SUBSET_LEX_NEXT generates all subsets of an N set.");
        Console.WriteLine("  The user can impose a restriction on the");
        Console.WriteLine("  maximum size of the subsets.");
        Console.WriteLine("");
        Console.WriteLine("  Here, we require the subsets to be no larger");
        Console.WriteLine("  than NDIM = " + NDIM + "");

        k = 0;

        for (;;)
        {
            ltest = k == NDIM;

            Subset.subset_lex_next(n, ltest, NDIM, ref k, ref a);

            switch (k)
            {
                case > 0:
                {
                    string cout = "  ";
                    for (i = 0; i < k; i++)
                    {
                        cout += a[i].ToString(CultureInfo.InvariantCulture).PadLeft(2) + "  ";
                    }

                    Console.WriteLine(cout);
                    break;
                }
                default:
                    Console.WriteLine("  The empty set.");
                    break;
            }

            if (k == 0)
            {
                break;
            }

        }
    }

    public static void subset_gray_next_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SUBSET_GRAY_NEXT_TEST tests SUBSET_GRAY_NEXT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 5;

        int[] a = new int[N];
        int i;
        int j;
        int iadd = 0;
        bool more;
        int ncard = 0;

        Console.WriteLine("");
        Console.WriteLine("SUBSET_GRAY_NEXT_TEST");
        Console.WriteLine("  SUBSET_GRAY_NEXT generates all subsets of an N set");
        Console.WriteLine("  using the Gray code ordering:");
        Console.WriteLine("  0 0 1 0 1 means the subset contains 3 and 5.");
        Console.WriteLine("");
        Console.WriteLine("  Gray code");
        Console.WriteLine("");

        more = false;
        j = 0;

        for (;;)
        {
            Subset.subset_gray_next(N, ref a, ref more, ref ncard, ref iadd);

            j += 1;
            for (i = 0; i < N; i++)
            {
            }

            Console.WriteLine("");

            if (!more)
            {
                break;
            }

        }


    }

    public static void subset_random_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SUBSET_RANDOM_TEST tests SUBSET_RANDOM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 5;

        int[] a = new int[N];
        int i;
        int j;
        int seed;

        Console.WriteLine("");
        Console.WriteLine("SUBSET_RANDOM_TEST");
        Console.WriteLine("  SUBSET_RANDOM picks a subset at random.");
        Console.WriteLine("  The number of elements in the main set is " + N + "");
        Console.WriteLine("");

        seed = 123456789;

        for (j = 1; j <= 5; j++)
        {
            a = Subset.subset_random(N, ref seed);

            string cout = j.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "    ";
            for (i = 0; i < N; i++)
            {
                cout += a[i].ToString(CultureInfo.InvariantCulture).PadLeft(2);
            }

            Console.WriteLine(cout);

        }

    }

    public static void subset_gray_rank_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SUBSET_GRAY_RANK_TEST tests SUBSET_GRAY_RANK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 5;

        int[] a = { 1, 0, 1, 1, 0 };
        int i;
        int rank;

        Console.WriteLine("");
        Console.WriteLine("SUBSET_GRAY_RANK_TEST");
        Console.WriteLine("  SUBSET_GRAY_RANK returns rank of a subset of an N set");
        Console.WriteLine("  using the Gray code ordering.");
        Console.WriteLine("");
        Console.WriteLine("  For N = " + N + ", the subset is:");

        for (i = 0; i < N; i++)
        {
        }

        Console.WriteLine("");

        rank = Ranking.subset_gray_rank(N, a);

        Console.WriteLine("");
        Console.WriteLine("  The rank is " + rank + "");
    }

    public static void subset_gray_unrank_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SUBSET_GRAY_UNRANK_TEST tests SUBSET_GRAY_UNRANK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 5;

        int[] a = new int[N];
        int i;
        int rank;

        Console.WriteLine("");
        Console.WriteLine("SUBSET_GRAY_UNRANK_TEST");
        Console.WriteLine("  SUBSET_GRAY_UNRANK finds the subset of an N set");
        Console.WriteLine("  of a given rank under the Gray code ordering.");
        Console.WriteLine("");
        Console.WriteLine("  N is " + N + "");
        Console.WriteLine("");
        Console.WriteLine("  Rank   Subset");
        Console.WriteLine("");

        for (rank = 1; rank <= 10; rank++)
        {
            Ranking.subset_gray_unrank(rank, N, ref a);

            string cout = "  "
                          + rank.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  ";
            for (i = 0; i < N; i++)
            {
                cout += a[i].ToString(CultureInfo.InvariantCulture).PadLeft(2);
            }

            Console.WriteLine(cout);

        }
    }

}