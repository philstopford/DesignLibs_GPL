using System;
using Burkardt.RankingNS;
using Burkardt.Types;

namespace ComboTest
{
    partial class Program
    {
        static void npart_enum_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NPART_ENUM_TEST tests NPART_ENUM.
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
            int part_num;

            Console.WriteLine("");
            Console.WriteLine("NPART_ENUM_TEST");
            Console.WriteLine("  NPART_ENUM enumerates partitions of N into PART_NUM parts.");
            Console.WriteLine("");
            Console.WriteLine("   PART_NUM:  1       2       3       4       5       6");
            Console.WriteLine("   N");

            for (n = 0; n <= 10; n++)
            {
                string cout = "  " + n.ToString().PadLeft(2) + ":  ";
                for (part_num = 1; part_num <= Math.Min(n, 6); part_num++)
                {
                    cout += "  " + Ranking.npart_enum(n, part_num).ToString().PadLeft(6);
                }

                Console.WriteLine(cout);
            }
        }

        static void npart_rsf_lex_rank_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NPART_RSF_LEX_RANK_TEST tests NPART_RSF_LEX_RANK.
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
            int npart;
            int rank;
            int[] t =  {
                1, 5, 6
            }
            ;

            n = 12;
            npart = 3;

            Console.WriteLine("");
            Console.WriteLine("NPART_RSF_LEX_RANK_TEST:");
            Console.WriteLine("  NPART_RSF_LEX_RANK ranks");
            Console.WriteLine("  partitions of N with NPART parts");
            Console.WriteLine("  in reverse standard form:");

            typeMethods.i4vec_transpose_print(npart, t, "  Element to be ranked:");

            rank = Ranking.npart_rsf_lex_rank(n, npart, t);

            Console.WriteLine("");
            Console.WriteLine("  Rank is computed as " + rank + "");
        }

        static void npart_rsf_lex_successor_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NPART_RSF_LEX_SUCCESSOR_TEST tests NPART_RSF_LEX_SUCCESSOR.
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
            int i;
            int n;
            int npart = 3;
            int rank;
            int rank_old;
            int[] t;

            n = 12;

            Console.WriteLine("");
            Console.WriteLine("NPART_RSF_LEX_SUCCESSOR_TEST");
            Console.WriteLine("  NPART_RSF_LEX_SUCCESSOR lists");
            Console.WriteLine("  partitions of N with NPART parts");
            Console.WriteLine("  in reverse standard form:");

            t = new int[npart];

            rank = -1;

            for (;;)
            {
                rank_old = rank;

                Ranking.npart_rsf_lex_successor(n, npart, ref t, ref rank);

                if (rank <= rank_old)
                {
                    break;
                }

                string cout = "  " + rank.ToString().PadLeft(4);
                for (i = 0; i < npart; i++)
                {
                    cout += "  " + t[i].ToString().PadLeft(4);
                }

                Console.WriteLine(cout);
            }
        }

        static void npart_rsf_lex_unrank_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NPART_RSF_LEX_UNRANK_TEST tests NPART_RSF_LEX_UNRANK.
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
            int npart;
            int rank;
            int[] t;


            Console.WriteLine("");
            Console.WriteLine("NPART_RSF_LEX_UNRANK_TEST");
            Console.WriteLine("  NPART_RSF_LEX_UNRANK unranks");
            Console.WriteLine("  partitions of N with NPART parts");
            Console.WriteLine("  in reverse standard form:");

            rank = 4;
            n = 12;
            npart = 3;

            t = Ranking.npart_rsf_lex_unrank(rank, n, npart);

            typeMethods.i4vec_transpose_print(npart, t, "  The element of rank 4:");
        }

        static void npart_rsf_lex_random_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NPART_RSF_LEX_RANDOM_TEST tests NPART_RSF_LEX_RANDOM;
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
            int i;
            int n = 12;
            int npart = 3;
            int seed = 123456789;
            int[] t;

            Console.WriteLine("");
            Console.WriteLine("NPART_RSF_LEX_RANDOM_TEST");
            Console.WriteLine("  NPART_RSF_LEX_RANDOM produces random examples of");
            Console.WriteLine("  partitions of N with NPART parts");
            Console.WriteLine("  in reverse standard form.");
            Console.WriteLine("");

            for (i = 1; i <= 10; i++)
            {
                t = Ranking.npart_rsf_lex_random(n, npart, ref seed);
                typeMethods.i4vec_transpose_print(npart, t, "");
            }
        }

        static void npart_sf_lex_successor_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NPART_SF_LEX_SUCCESSOR_TEST tests NPART_SF_LEX_SUCCESSOR;
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
            int i;
            int n = 12;
            int npart = 3;
            int rank;
            int rank_old;
            int[] t;

            Console.WriteLine("");
            Console.WriteLine("NPART_SF_LEX_SUCCESSOR_TEST");
            Console.WriteLine("  NPART_SF_LEX_SUCCESSOR lists");
            Console.WriteLine("  Partitions of N with NPART parts.");

            t = new int[npart];

            rank = -1;

            for (;;)
            {
                rank_old = rank;

                Ranking.npart_sf_lex_successor(n, npart, ref t, ref rank);

                if (rank <= rank_old)
                {
                    break;
                }

                string cout = "  " + rank.ToString().PadLeft(4);
                for (i = 0; i < npart; i++)
                {
                    cout += "  " + t[i].ToString().PadLeft(4);
                }

                Console.WriteLine("");
            }
        }

        static void npart_table_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NPART_TABLE_TEST tests NPART_TABLE.
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
            int i;
            int j;
            int maxn = 10;
            int maxpart = 5;
            int[] p;

            Console.WriteLine("");
            Console.WriteLine("NPART_TABLE_TEST");
            Console.WriteLine("  NPART_TABLE tabulates partitions");
            Console.WriteLine("  of N with NPART parts;");

            p = Ranking.npart_table(maxn, maxpart);

            Console.WriteLine("");
            Console.WriteLine("   I     1      2      3      4      5");
            Console.WriteLine("");

            for (i = 0; i <= maxn; i++)
            {
                string cout = "  " + i.ToString().PadLeft(2);
                for (j = 0; j <= maxpart; j++)
                {
                    cout += "  " + p[i + j * (maxn + 1)].ToString().PadLeft(4);
                }

                Console.WriteLine(cout);
            }
        }
    }
}