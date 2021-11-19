using System;
using Burkardt.RankingNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace ComboTest;

internal partial class Program
{
    private static void tree_check_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TREE_CHECK_TEST tests TREE_CHECK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 December 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        bool check;
        int n = 0;
        int test;
        int[] t = new int[1];
        int[] t2 =  {
                1, 2,
                2, 3
            }
            ;
        int[] t3 =  {
                1, 2,
                3, 4,
                4, 5,
                5, 3
            }
            ;
        int[] t4 =  {
                1, 3,
                2, 3,
                3, 4,
                4, 5,
                5, 6
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("TREE_CHECK TEST");
        Console.WriteLine("  TREE_CHECK checks a tree.");
        Console.WriteLine("");
        Console.WriteLine("  Check?");
        Console.WriteLine("");

        for (test = 1; test <= 4; test++)
        {
            switch (test)
            {
                case 1:
                    n = 0;
                    t = null;
                    break;
                case 2:
                    n = 3;
                    t = typeMethods.i4vec_copy_new(2 * (n - 1), t2);
                    break;
                case 3:
                    n = 5;
                    t = typeMethods.i4vec_copy_new(2 * (n - 1), t3);
                    break;
                case 4:
                    n = 6;
                    t = typeMethods.i4vec_copy_new(2 * (n - 1), t4);
                    break;
            }

            check = Ranking.tree_check(n, t);
            Console.WriteLine("       " + check.ToString().PadLeft(2) + "");
            typeMethods.i4mat_print(2, n - 1, t, "  Tree:");
        }
    }

    private static void tree_enum_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TREE_ENUM_TEST tests TREE_ENUM.
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
        int tree_num;

        Console.WriteLine("");
        Console.WriteLine("TREE_ENUM_TEST");
        Console.WriteLine("  TREE_ENUM enumerates trees on N nodes.");

        for (n = 0; n <= 10; n++)
        {
            tree_num = Ranking.tree_enum(n);
            Console.WriteLine("  " + n.ToString().PadLeft(2)
                                   + "  " + tree_num.ToString().PadLeft(6) + "");
        }
    }

    private static void tree_rank_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TREE_RANK_TEST tests TREE_RANK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 December 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n;
        int rank;
        int[] t =  {
                4, 3,
                3, 1,
                2, 1
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("TREE_RANK_TEST");
        Console.WriteLine("  TREE_RANK ranks trees.");


        n = 4;
        typeMethods.i4mat_print(2, n - 1, t, "  Tree to be ranked:");

        rank = Ranking.tree_rank(n, t);

        Console.WriteLine("");
        Console.WriteLine("  Rank is computed as " + rank + "");
    }

    private static void tree_successor_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TREE_SUCCESSOR_TEST tests TREE_SUCCESSOR.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 December 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int j;
        int n = 4;
        int rank;
        int rank_old;
        int[] t = new int[1];

        Console.WriteLine("");
        Console.WriteLine("TREE_SUCCESSOR_TEST");
        Console.WriteLine("  TREE_SUCCESSOR lists trees.");

        t = new int[2 * (n - 1)];

        rank = -1;

        for (;;)
        {
            rank_old = rank;

            Ranking.tree_successor(n, ref t, ref rank);

            if (rank <= rank_old)
            {
                break;
            }

            string cout = "  " + rank.ToString().PadLeft(4);
            for (j = 0; j < n - 1; j++)
            {
                cout += "  " + t[0 + j * 2].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);
            string cout2 = "  " + "    ";
            for (j = 0; j < n - 1; j++)
            {
                cout2 += "  " + t[1 + j * 2].ToString().PadLeft(4);
            }

            Console.WriteLine(cout2);
        }
    }

    private static void tree_to_pruefer_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TREE_TO_PRUEFER_TEST tests TREE_TO_PRUEFER.
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
        Console.WriteLine("TREE_TO_PRUEFER_TEST");
        Console.WriteLine("  TREE_TO_PRUEFER converts a tree to a Pruefer code.");

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
                Console.WriteLine("  " + j.ToString().PadLeft(2)
                                       + "  " + t[0 + j * 2].ToString().PadLeft(4)
                                       + "  " + t[1 + j * 2].ToString().PadLeft(4) + "");
            }

            //
            //  Convert the tree to a Pruefer code.
            //

            p = Ranking.tree_to_pruefer(n, t);

            Console.WriteLine("");
            typeMethods.i4vec_transpose_print(n - 2, p, "  Pruefer code:");
        }
    }

    private static void tree_unrank_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TREE_UNRANK_TEST tests TREE_UNRANK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 December 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n;
        int rank;
        int[] t;

        Console.WriteLine("");
        Console.WriteLine("TREE_UNRANK_TEST");
        Console.WriteLine("  TREE_UNRANK unranks trees.");

        rank = 8;
        n = 4;

        t = Ranking.tree_unrank(rank, n);

        typeMethods.i4mat_print(2, n - 1, t, "  The tree of rank 8:");
    }
}