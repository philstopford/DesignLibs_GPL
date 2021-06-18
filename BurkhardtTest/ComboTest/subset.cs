using System;
using Burkardt.RankingNS;
using Burkardt.Types;

namespace ComboTest
{
    partial class Program
    {
        static void subset_check_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SUBSET_CHECK_TEST tests SUBSET_CHECK.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    10 January 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            bool check;
            int n = 0;
            int[] s = new int[1];
            int[] s1 =  {
            }
            ;
            int[] s2 =  {
                1, 2, 0
            }
            ;
            int[] s3 =  {
                1, 0, 0, 1, 0
            }
            ;
            int test;

            Console.WriteLine("");
            Console.WriteLine("SUBSET_CHECK TEST");
            Console.WriteLine("  SUBSET_CHECK checks a subset.");

            for (test = 1; test <= 3; test++)
            {
                if (test == 1)
                {
                    n = 0;
                    s = typeMethods.i4vec_copy_new(n, s1);
                }
                else if (test == 2)
                {
                    n = 3;
                    s = typeMethods.i4vec_copy_new(n, s2);
                }
                else if (test == 3)
                {
                    n = 5;
                    s = typeMethods.i4vec_copy_new(n, s3);
                }

                check = Ranking.subset_check(n, s);
                typeMethods.i4vec_transpose_print(n, s, "  Subset:");
                Console.WriteLine("  Check = " + check + "");
            }
        }

        static void subset_colex_rank_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SUBSET_COLEX_RANK_TEST tests SUBSET_COLEX_RANK.
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
            int n = 5;
            int rank;
            int[] t =  {
                0, 1, 0, 1, 0
            }
            ;

            Console.WriteLine("");
            Console.WriteLine("SUBSET_COLEX_RANK_TEST");
            Console.WriteLine("  SUBSET_COLEX_RANK ranks subsets of a set,");
            Console.WriteLine("  using the colexicographic ordering.");

            typeMethods.i4vec_transpose_print(n, t, "  Element to be ranked:");

            rank = Ranking.subset_colex_rank(n, t);

            Console.WriteLine("");
            Console.WriteLine("  Rank is computed as " + rank + "");
        }

        static void subset_colex_successor_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SUBSET_COLEX_SUCCESSOR_TEST tests SUBSET_COLEX_SUCCESSOR.
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
            int n = 5;
            int rank;
            int rank_old;
            int[] t;

            Console.WriteLine("");
            Console.WriteLine("SUBSET_COLEX_SUCCESSOR_TEST");
            Console.WriteLine("  SUBSET_COLEX_SUCCESSOR lists subsets of a set,");
            Console.WriteLine("  using the colexicographic ordering.");

            t = new int[n];

            rank = -1;

            for (;;)
            {
                rank_old = rank;

                Ranking.subset_colex_successor(n, ref t, ref rank);

                if (rank <= rank_old)
                {
                    break;
                }

                string cout = "  " + rank.ToString().PadLeft(4);
                for (i = 0; i < n; i++)
                {
                    cout += "  " + t[i].ToString().PadLeft(4);
                }

                Console.WriteLine(cout);
            }
        }

        static void subset_colex_unrank_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SUBSET_COLEX_UNRANK_TEST tests SUBSET_COLEX_UNRANK.
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
            int rank;
            int[] t;

            Console.WriteLine("");
            Console.WriteLine("SUBSET_COLEX_UNRANK_TEST");
            Console.WriteLine("  SUBSET_COLEX_UNRANK unranks subsets of a set,");
            Console.WriteLine("  using the colexicographic ordering.");

            rank = 10;
            n = 5;

            t = Ranking.subset_colex_unrank(rank, n);

            typeMethods.i4vec_transpose_print(n, t, "  The element of rank 10:");
        }

        static void subset_complement_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SUBSET_COMPLEMENT_TEST tests SUBSET_COMPLEMENT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 January 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int n;
            int[] s1;
            int[] s2;
            int seed;

            Console.WriteLine("");
            Console.WriteLine("SUBSET_COMPLEMENT_TEST");
            Console.WriteLine("  SUBSET_COMPLEMENT returns the complement of a subset.");
            Console.WriteLine("");
            n = 5;
            seed = 123456789;

            s1 = Ranking.subset_random(n, ref seed);
            typeMethods.i4vec_transpose_print(n, s1, "  Subset S1:            ");

            s2 = Ranking.subset_complement(n, s1);
            typeMethods.i4vec_transpose_print(n, s2, "  S2 = complement of S1:");
        }

        static void subset_distance_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SUBSET_DISTANCE_TEST tests SUBSET_DISTANCE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    08 January 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int distance;
            int n;
            int[] s1;
            int[] s2;
            int seed;

            Console.WriteLine("");
            Console.WriteLine("SUBSET_DISTANCE_TEST");
            Console.WriteLine("  SUBSET_DISTANCE returns the distance between two subsets.");
            Console.WriteLine("");

            n = 5;
            seed = 123456789;

            s1 = Ranking.subset_random(n, ref seed);
            typeMethods.i4vec_transpose_print(n, s1, "  Subset S1:");

            s2 = Ranking.subset_random(n, ref seed);
            typeMethods.i4vec_transpose_print(n, s2, "  Subset S2:");

            distance = Ranking.subset_distance(n, s1, s2);
            Console.WriteLine("");
            Console.WriteLine("  Distance = " + distance + "");
        }

        static void subset_enum_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SUBSET_ENUM_TEST tests SUBSET_ENUM.
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
            int subset_num;

            Console.WriteLine("");
            Console.WriteLine("SUBSET_ENUM_TEST");
            Console.WriteLine("  SUBSET_ENUM enumerates subsets of a set of N items.");
            Console.WriteLine("");

            for (n = 0; n <= 10; n++)
            {
                subset_num = Ranking.subset_enum(n);
                Console.WriteLine("  " + n.ToString().PadLeft(2)
                    + "  " + subset_num.ToString().PadLeft(6) + "");
            }
        }

        static void subset_intersect_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SUBSET_INTERSECT_TEST tests SUBSET_INTERSECT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    09 January 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int n;
            int[] s1;
            int[] s2;
            int[] s3;
            int seed;

            Console.WriteLine("");
            Console.WriteLine("SUBSET_INTERSECT_TEST");
            Console.WriteLine("  SUBSET_INTERSECT returns the intersection of two subsets.");
            Console.WriteLine("");

            n = 7;
            seed = 123456789;

            s1 = Ranking.subset_random(n, ref seed);
            typeMethods.i4vec_transpose_print(n, s1, "  Subset S1:");

            s2 = Ranking.subset_random(n, ref seed);
            typeMethods.i4vec_transpose_print(n, s2, "  Subset S2:");

            s3 = Ranking.subset_intersect(n, s1, s2);
            typeMethods.i4vec_transpose_print(n, s3, "  Intersect:");
        }

        static void subset_lex_rank_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SUBSET_LEX_RANK_TEST tests SUBSET_LEX_RANK.
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
            int rank;
            int[] t =  {
                0, 1, 0, 1, 0
            }
            ;

            Console.WriteLine("");
            Console.WriteLine("SUBSET_LEX_RANK_TEST");
            Console.WriteLine("  SUBSET_LEX_RANK ranks subsets of a set,");
            Console.WriteLine("  using the lexicographic ordering.");

            n = 5;
            typeMethods.i4vec_transpose_print(n, t, "  Element to be ranked:");

            rank = Ranking.subset_lex_rank(n, t);

            Console.WriteLine("");
            Console.WriteLine("  Rank is computed as " + rank + "");
        }

        static void subset_lex_successor_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SUBSET_LEX_SUCCESSOR_TEST tests SUBSET_LEX_SUCCESSOR.
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
            int n = 5;
            int rank;
            int rank_old;
            int[] t;

            Console.WriteLine("");
            Console.WriteLine("SUBSET_LEX_SUCCESSOR_TEST");
            Console.WriteLine("  SUBSET_LEX_SUCCESSOR lists,");
            Console.WriteLine("  subsets of a set,");
            Console.WriteLine("  using the lexicographic ordering,");

            t = new int[n];

            rank = -1;

            for (;;)
            {
                rank_old = rank;

                Ranking.subset_lex_successor(n, ref t, ref rank);

                if (rank <= rank_old)
                {
                    break;
                }

                string cout = "  " + rank.ToString().PadLeft(4);
                for (i = 0; i < n; i++)
                {
                    cout += "  " + t[i].ToString().PadLeft(4);
                }

                Console.WriteLine(cout);
            }
        }

        static void subset_lex_unrank_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SUBSET_LEX_UNRANK_TEST tests SUBSET_LEX_UNRANK.
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
            int n = 5;
            int rank;
            int[] t;

            Console.WriteLine("");
            Console.WriteLine("SUBSET_LEX_UNRANK_TEST");
            Console.WriteLine("  SUBSET_LEX_UNRANK unranks");
            Console.WriteLine("  subsets of a set,");
            Console.WriteLine("  using the lexicographic ordering:");

            rank = 10;

            t = Ranking.subset_lex_unrank(rank, n);

            typeMethods.i4vec_transpose_print(n, t, "  The element of rank 10:");
        }

        static void subset_random_test()

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
            //    25 December 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i;
            int n;
            int[] s;
            int seed;

            Console.WriteLine("");
            Console.WriteLine("SUBSET_RANDOM_TEST");
            Console.WriteLine("  SUBSET_RANDOM returns a random subset.");

            n = 5;
            seed = 123456789;

            for (i = 0; i < 10; i++)
            {
                s = Ranking.subset_random(n, ref seed);
                typeMethods.i4vec_transpose_print(n, s, "  Subset:");
            }
        }

        static void subset_union_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SUBSET_UNION_TEST tests SUBSET_UNION.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    09 January 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int n;
            int[] s1;
            int[] s2;
            int[] s3;
            int seed;

            Console.WriteLine("");
            Console.WriteLine("SUBSET_UNION_TEST");
            Console.WriteLine("  SUBSET_UNION returns the union of two subsets.");
            Console.WriteLine("");

            n = 7;
            seed = 123456789;

            s1 = Ranking.subset_random(n, ref seed);
            typeMethods.i4vec_transpose_print(n, s1, "  Subset S1:");

            s2 = Ranking.subset_random(n, ref seed);
            typeMethods.i4vec_transpose_print(n, s2, "  Subset S2:");

            s3 = Ranking.subset_union(n, s1, s2);
            typeMethods.i4vec_transpose_print(n, s3, "  Union:    ");

        }

        static void subset_weight_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SUBSET_WEIGHT_TEST tests SUBSET_WEIGHT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 January 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int n;
            int[] s;
            int seed;
            int weight;

            Console.WriteLine("");
            Console.WriteLine("SUBSET_WEIGHT_TEST");
            Console.WriteLine("  SUBSET_WEIGHT returns the weight of a subset.");

            n = 5;
            seed = 123456789;

            s = Ranking.subset_random(n, ref seed);
            typeMethods.i4vec_transpose_print(n, s, "  Subset S:");

            weight = Ranking.subset_weight(n, s);
            Console.WriteLine("");
            Console.WriteLine("  Weight = " + weight + "");
        }

        static void subset_xor_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SUBSET_XOR_TEST tests SUBSET_XOR.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    09 January 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int n;
            int[] s1;
            int[] s2;
            int[] s3;
            int seed;

            Console.WriteLine("");
            Console.WriteLine("SUBSET_XOR_TEST");
            Console.WriteLine("  SUBSET_XOR returns the exclusive OR of two subsets.");
            Console.WriteLine("");

            n = 7;
            seed = 123456789;

            s1 = Ranking.subset_random(n, ref seed);
            typeMethods.i4vec_transpose_print(n, s1, "  Subset S1:");

            s2 = Ranking.subset_random(n, ref seed);
            typeMethods.i4vec_transpose_print(n, s2, "  Subset S2:");

            s3 = Ranking.subset_xor(n, s1, s2);
            typeMethods.i4vec_transpose_print(n, s3, "  XOR:      ");
        }

        static void subsetsumswap_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SUBSETSUMSWAP_TEST tests SUBSETSUM_SWAP.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    27 July 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 7;

            int[] a =  {
                12, 8, 11, 30, 8, 3, 7
            }
            ;
            int i;
            int[] index = new int[N];
            int n = N;
            int sum_achieved;
            int sum_desired = 17;

            Console.WriteLine("");
            Console.WriteLine("SUBSETSUMSWAP_TEST");
            Console.WriteLine("  SUBSETSUM_SWAP seeks a solution of the subset");
            Console.WriteLine("  sum problem using pair swapping.");
            Console.WriteLine("");
            Console.WriteLine("  The desired sum is " + sum_desired + "");

            sum_achieved = Ranking.subsetsum_swap(n, ref a, sum_desired, ref index);

            Console.WriteLine("");
            Console.WriteLine("    A(I), INDEX(I)");
            Console.WriteLine("");

            for (i = 0; i < n; i++)
            {
                Console.WriteLine("  " + a[i].ToString().PadLeft(5)
                    + "  " + index[i].ToString().PadLeft(5) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("  The achieved sum is " + sum_achieved + "");
        }
    }
}