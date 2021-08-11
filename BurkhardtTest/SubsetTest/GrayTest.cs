using System;
using Burkardt.RankingNS;

namespace SubsetTest
{
    public static class GrayTest
    {
        public static void gray_next_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GRAY_NEXT_TEST tests GRAY_NEXT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 June 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int[] a = new int[4];
            int change;
            int[] g = new int[4];
            int i;
            int k;
            int n = 4;

            Console.WriteLine("");
            Console.WriteLine("GRAY_NEXT_TEST");
            Console.WriteLine("  GRAY_NEXT returns the index of the single item");
            Console.WriteLine("  to be changed in order to get the next Gray code.");

            Console.WriteLine("");
            Console.WriteLine("   K  Switch  Gray Code");
            Console.WriteLine("");

            change = -n;
            k = 0;

            for (;;)
            {
                Ranking.gray_next(n, ref change, ref k, ref a);

                if (change == -n)
                {
                    break;
                }
                else if (change == 0)
                {
                    for (i = 0; i < n; i++)
                    {
                        g[i] = 0;
                    }
                }
                else
                {
                    g[Math.Abs(change) - 1] = 1 - g[Math.Abs(change) - 1];
                }

                string cout = "  "
                              + k.ToString().PadLeft(2) + "  "
                              + change.ToString().PadLeft(6) + "  ";
                for (i = 0; i < n; i++)
                {
                    cout += g[i].ToString().PadLeft(1);
                }

                Console.WriteLine(cout);
            }
        }

        public static void gray_rank_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GRAY_RANK_TEST tests GRAY_RANK.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 January 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int gray;
            int rank;
            int rank2;

            Console.WriteLine("");
            Console.WriteLine("GRAY_RANK_TEST");
            Console.WriteLine("  GRAY_RANK ranks a Gray code;");
            Console.WriteLine("");
            Console.WriteLine("    R  =                        RANK");
            Console.WriteLine("    G  =            GRAY_UNRANK(RANK)");
            Console.WriteLine("    R2 =  GRAY_RANK(GRAY_UNRANK(RANK))");
            Console.WriteLine("");
            Console.WriteLine("    R    G    R2");
            Console.WriteLine("");

            for (rank = 0; rank <= 24; rank++)
            {
                gray = Ranking.gray_unrank(rank);

                rank2 = Ranking.gray_rank(gray);

                Console.WriteLine(rank.ToString().PadLeft(9) + "  "
                                                             + gray.ToString().PadLeft(9) + "  "
                                                             + rank2.ToString().PadLeft(9) + "");
            }
        }

        public static void gray_rank2_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GRAY_RANK2_TEST tests GRAY_RANK2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 January 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int gray;
            int rank;
            int rank2;

            Console.WriteLine("");
            Console.WriteLine("GRAY_RANK2_TEST");
            Console.WriteLine("  GRAY_RANK2 ranks a Gray code;");
            Console.WriteLine("");
            Console.WriteLine("    R  =                          RANK");
            Console.WriteLine("    G  =             GRAY_UNRANK2(RANK)");
            Console.WriteLine("    R2 =  GRAY_RANK2(GRAY_UNRANK2(RANK))");
            Console.WriteLine("");
            Console.WriteLine("    R    G    R2");
            Console.WriteLine("");

            for (rank = 0; rank <= 24; rank++)
            {
                gray = Ranking.gray_unrank2(rank);

                rank2 = Ranking.gray_rank2(gray);

                Console.WriteLine(rank.ToString().PadLeft(9) + "  "
                                                             + gray.ToString().PadLeft(9) + "  "
                                                             + rank2.ToString().PadLeft(9) + "");
            }
        }

        public static void gray_unrank_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GRAY_UNRANK_TEST tests GRAY_UNRANK.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 January 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int gray;
            int rank;
            int rank2;

            Console.WriteLine("");
            Console.WriteLine("GRAY_UNRANK_TEST");
            Console.WriteLine("  GRAY_UNRANK unranks a Gray code.");
            Console.WriteLine("");
            Console.WriteLine("    R  =                        RANK");
            Console.WriteLine("    G  =            GRAY_UNRANK(RANK)");
            Console.WriteLine("    R2 =  GRAY_RANK(GRAY_UNRANK(RANK))");
            Console.WriteLine("");
            Console.WriteLine("    R    G    R2");
            Console.WriteLine("");

            for (rank = 0; rank <= 24; rank++)
            {
                gray = Ranking.gray_unrank(rank);

                rank2 = Ranking.gray_rank(gray);

                Console.WriteLine(rank.ToString().PadLeft(9) + "  "
                                                             + gray.ToString().PadLeft(9) + "  "
                                                             + rank2.ToString().PadLeft(9) + "");
            }
        }

        public static void gray_unrank2_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GRAY_UNRANK2_TEST tests GRAY_UNRANK2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 January 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int gray;
            int rank;
            int rank2;

            Console.WriteLine("");
            Console.WriteLine("GRAY_UNRANK2_TEST");
            Console.WriteLine("  GRAY_UNRANK2 unranks a Gray code.");
            Console.WriteLine("");
            Console.WriteLine("    R  =                          RANK");
            Console.WriteLine("    G  =             GRAY_UNRANK2(RANK)");
            Console.WriteLine("    R2 =  GRAY_RANK2(GRAY_UNRANK2(RANK))");
            Console.WriteLine("");
            Console.WriteLine("    R    G    R2");
            Console.WriteLine("");

            for (rank = 0; rank <= 24; rank++)
            {
                gray = Ranking.gray_unrank2(rank);

                rank2 = Ranking.gray_rank2(gray);

                Console.WriteLine(rank.ToString().PadLeft(9) + "  "
                                                             + gray.ToString().PadLeft(9) + "  "
                                                             + rank2.ToString().PadLeft(9) + "");
            }
        }

    }
}