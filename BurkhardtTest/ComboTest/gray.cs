using System;
using Burkardt.RankingNS;
using Burkardt.Types;

namespace ComboTest
{
    partial class Program
    {
        static void gray_code_check_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GRAY_CODE_CHECK_TEST tests GRAY_CODE_CHECK.
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
            bool check;
            int i;
            int n;
            int test;
            int[] t = new int[1];
            int[] t3 =  {
                1, 1, 1, 1, 1
            }
            ;

            Console.WriteLine("");
            Console.WriteLine("GRAY_CODE_CHECK TEST");
            Console.WriteLine("  GRAY_CODE_CHECK checks N and T(1:N).");
            Console.WriteLine("");
            Console.WriteLine("  Check?   N    T(1:N)");
            Console.WriteLine("");

            for (test = 1; test <= 3; test++)
            {
                n = 5;

                if (test == 1)
                {
                    t = typeMethods.i4vec_copy_new(n, t3);
                }
                else if (test == 2)
                {
                    t = typeMethods.i4vec_copy_new(n, t3);
                }
                else if (test == 3)
                {
                    t = typeMethods.i4vec_copy_new(n, t3);
                }

                check = typeMethods.gray_code_check(n, t3);
                string cout = "      " + check.ToString().PadLeft(1)
                                  + "  " + n.ToString().PadLeft(2)
                                  + ":  ";
                for (i = 0; i < n; i++)
                {
                    cout += "  " + t[i].ToString().PadLeft(2);
                }

                Console.WriteLine(cout);
            }
        }

        static void gray_code_enum_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GRAY_CODE_ENUM_TEST tests GRAY_CODE_ENUM.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    22 November 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int n;
            int ngray;

            Console.WriteLine("");
            Console.WriteLine("GRAY_CODE_ENUM_TEST");
            Console.WriteLine("  GRAY_CODE_ENUM enumerates Gray codes with N elements.");
            Console.WriteLine("");
            Console.WriteLine("   N   Enum(N)");
            Console.WriteLine("");
            for (n = 0; n <= 10; n++)
            {
                ngray = typeMethods.gray_code_enum(n);
                Console.WriteLine("  " + n.ToString().PadLeft(2)
                                  + "  " + ngray.ToString().PadLeft(6) + "");
            }
        }

        static void gray_code_rank_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GRAY_CODE_RANK_TEST tests GRAY_CODE_RANK.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    22 November 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int n = 5;
            int rank;
            int[] t =  {
                1, 1, 0, 0, 0
            }
            ;

            Console.WriteLine("");
            Console.WriteLine("GRAY_CODE_RANK_TEST");
            Console.WriteLine("  GRAY_CODE_RANK ranks a Gray code.");

            rank = Ranking.gray_code_rank(n, t);

            typeMethods.i4vec_transpose_print(n, t, "  Element to be ranked:");

            Console.WriteLine("");
            Console.WriteLine("  Computed rank: " + rank + "");
        }

        static void gray_code_successor_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GRAY_CODE_SUCCESSOR_TEST tests GRAY_CODE_SUCCESSOR.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    22 November 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i;
            int n;
            int rank;
            int rank_old;
            int[] t;

            Console.WriteLine("");
            Console.WriteLine("GRAY_CODE_SUCCESSOR_TEST");
            Console.WriteLine("  GRAY_CODE_SUCCESSOR lists Gray codes one by one.");

            n = 5;
            t = new int[n];
            rank = -1;

            for (;;)
            {
                rank_old = rank;

                typeMethods.gray_code_successor(n, ref t, ref rank);

                if (rank <= rank_old)
                {
                    break;
                }

                string cout = "  " + rank.ToString().PadLeft(4);
                for (i = 0; i < n; i++)
                {
                    cout += "  " + t[i].ToString().PadLeft(4);
                }

                Console.WriteLine("");
            }
        }

        static void gray_code_unrank_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GRAY_CODE_UNRANK_TEST tests GRAY_CODE_UNRANK.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    22 November 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i;
            int n;
            int ngray;
            int rank;
            int[] t;

            Console.WriteLine("");
            Console.WriteLine("GRAY_CODE_UNRANK_TEST");
            Console.WriteLine("  GRAY_CODE_UNRANK unranks a Gray code.");

            n = 5;
            ngray = typeMethods.gray_code_enum(n);
            rank = ngray / 2;

            t =Ranking. gray_code_unrank(rank, n);

            Console.WriteLine("");
            Console.WriteLine("  The element of rank " + rank + "");
            Console.WriteLine("");
            string cout = "";
            for (i = 0; i < n; i++)
            {
                cout += "  " + t[i].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);
       }
    }
}