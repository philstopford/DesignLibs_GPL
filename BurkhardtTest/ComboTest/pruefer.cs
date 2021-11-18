using System;
using Burkardt.RankingNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace ComboTest;

internal partial class Program
{
    private static void pruefer_check_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PRUEFER_CHECK_TEST tests PRUEFER_CHECK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 December 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        bool check;
        int i;
        int n = 0;
        int test;
        int[] p = new int[1];
        int[] p1 =  {
            }
            ;
        int[] p2 =  {
                1
            }
            ;
        int[] p3 =  {
                5, 2
            }
            ;
        int[] p4 =  {
                5, 1, 3
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("BAL_SEQ_CHECK TEST");
        Console.WriteLine("  BAL_SEQ_CHECK checks N and T(1:2*N).");
        Console.WriteLine("");
        Console.WriteLine("  Check?   N    T(1:2*N)");
        Console.WriteLine("");

        for (test = 1; test <= 4; test++)
        {
            switch (test)
            {
                case 1:
                    n = 2;
                    p = typeMethods.i4vec_copy_new(n - 2, p1);
                    break;
                case 2:
                    n = 3;
                    p = typeMethods.i4vec_copy_new(n - 2, p2);
                    break;
                case 3:
                    n = 4;
                    p = typeMethods.i4vec_copy_new(n - 2, p3);
                    break;
                case 4:
                    n = 5;
                    p = typeMethods.i4vec_copy_new(n - 2, p4);
                    break;
            }

            check = Ranking.pruefer_check(n, p);
            Console.WriteLine("    "
                              + "  " + check.ToString(CultureInfo.InvariantCulture).PadLeft(1)
                              + "  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(2));
            for (i = 0; i < n - 2; i++)
            {
            }

            Console.WriteLine("");
        }
    }

    private static void pruefer_enum_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PRUEFER_ENUM_TEST tests PRUEFER_ENUM.
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
        int n;
        int pruefer_num;

        Console.WriteLine("");
        Console.WriteLine("PRUEFER_ENUM_TEST");
        Console.WriteLine("  PRUEFER_ENUM enumerates trees on N nodes, using the Pruefer code.");
        Console.WriteLine("");
        Console.WriteLine("   N           #");
        Console.WriteLine("");

        for (n = 0; n <= 10; n++)
        {
            pruefer_num = Ranking.pruefer_enum(n);
            Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + pruefer_num.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }
    }

    private static void pruefer_rank_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PRUEFER_RANK_TEST tests PRUEFER_RANK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 December 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n = 4;
        int[] p =  {
                3, 1
            }
            ;
        int rank;

        Console.WriteLine("");
        Console.WriteLine("PRUEFER_RANK_TEST");
        Console.WriteLine("  PRUEFER_RANK ranks Pruefer codes.");

        typeMethods.i4vec_transpose_print(n - 2, p, "  Element to be ranked:");

        rank = Ranking.pruefer_rank(n, p);

        Console.WriteLine("");
        Console.WriteLine("  Rank is computed as " + rank + "");
    }

    private static void pruefer_successor_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PRUEFER_SUCCESSOR_TEST tests PRUEFER_SUCCESSOR.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 December 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int n = 4;
        int[] p = new int[1];
        int rank = 0;
        int rank_old;

        Console.WriteLine("");
        Console.WriteLine("PRUEFER_SUCCESSOR_TEST");
        Console.WriteLine("  PRUEFER_SUCCESSOR lists Pruefer codes.");

        p = new int[n - 2];

        rank = -1;

        for (;;)
        {
            rank_old = rank;

            Ranking.pruefer_successor(n, ref p, ref rank);

            if (rank <= rank_old)
            {
                break;
            }

            string cout = "  " + rank.ToString(CultureInfo.InvariantCulture).PadLeft(4);
            for (i = 0; i < n - 2; i++)
            {
                cout += "  " + p[i].ToString(CultureInfo.InvariantCulture).PadLeft(4);
            }

            Console.WriteLine(cout);
        }
    }

    private static void pruefer_to_tree_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PRUEFER_TO_TREE_TEST tests PRUEFER_TO_TREE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 December 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i4_hi;
        int i4_lo;
        int j;
        int n = 5;
        int[] p;
        int pruefer_num;
        int rank;
        int seed = 123456789;
        int[] t;
        int test;
        int test_num = 5;

        Console.WriteLine("");
        Console.WriteLine("PRUEFER_TO_TREE_TEST");
        Console.WriteLine("  PRUEFER_TO_TREE converts a Pruefer code to a tree;");

        pruefer_num = Ranking.pruefer_enum(n);

        i4_lo = 0;
        i4_hi = pruefer_num - 1;

        for (test = 1; test <= test_num; test++)
        {
            //
            //  Pick a "random" Pruefer code.
            //
            rank = UniformRNG.i4_uniform_ab(i4_lo, i4_hi, ref seed);

            p = Ranking.pruefer_unrank(rank, n);

            Console.WriteLine("");
            Console.WriteLine("  Random Pruefer code of rank " + rank + "");
            typeMethods.i4vec_transpose_print(n - 2, p, "");
            //
            //  Convert the Pruefer code to a tree.
            //
            t = Ranking.pruefer_to_tree_new(n, p);

            Console.WriteLine("");
            Console.WriteLine("  Edge list for the corresponding tree:");
            Console.WriteLine("");
            for (j = 0; j < n - 1; j++)
            {
                Console.WriteLine("  " + j.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                       + "  " + t[0 + j * 2].ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                       + "  " + t[1 + j * 2].ToString(CultureInfo.InvariantCulture).PadLeft(4) + "");
            }

            //
            //  Convert the tree to a Pruefer code.
            //

            p = Ranking.tree_to_pruefer(n, t);

            Console.WriteLine("");
            typeMethods.i4vec_transpose_print(n - 2, p, "  Pruefer code:");
        }
    }

    private static void pruefer_unrank_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PRUEFER_UNRANK_TEST tests PRUEFER_UNRANK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 December 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n;
        int[] p;
        int rank;

        Console.WriteLine("");
        Console.WriteLine("PRUEFER_UNRANK_TEST");
        Console.WriteLine("  PRUEFER_UNRANK unranks Pruefer codes.");

        rank = 8;
        n = 4;

        p = Ranking.pruefer_unrank(rank, n);

        typeMethods.i4vec_transpose_print(n - 2, p, "  The element of rank 8:");
    }
}