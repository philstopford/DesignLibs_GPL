using System;
using Burkardt.RankingNS;
using Burkardt.Types;

namespace Burkardt.ComboTest
{
    partial class Program
    {
        static void partition_greedy_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PARTITION_GREEDY_TEST tests PARTITION_GREEDY.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 July 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 10;

            int[] a1 =  {
                2, 10, 3, 8, 5, 7, 9, 5, 3, 2
            }
            ;
            int[] a2 =  {
                771, 121, 281, 854, 885, 734, 486, 1003, 83, 62
            }
            ;
            int i;
            int[] indx;
            int n = N;
            int[] sums = new int[2];

            Console.WriteLine("");
            Console.WriteLine("PARTITION_GREEDY_TEST");
            Console.WriteLine("  PARTITION_GREEDY partitions an integer vector into");
            Console.WriteLine("  two subsets with nearly equal sum.");
            Console.WriteLine("");

            indx = Ranking.partition_greedy(n, a1);

            Console.WriteLine("");
            Console.WriteLine("");
            Console.WriteLine("Data set #1 partitioned:");
            Console.WriteLine("");
            sums[0] = 0;
            sums[1] = 0;

            for (i = 0; i < n; i++)
            {
                if (indx[i] == 1)
                {
                    sums[0] = sums[0] + a1[i];
                    Console.WriteLine("  " + a1[i].ToString().PadLeft(4) + "");
                }
                else
                {
                    sums[1] = sums[1] + a1[i];
                    Console.WriteLine("  " + "    "
                        + "  " + a1[i].ToString().PadLeft(4) + "");
                }
            }

            Console.WriteLine("");
            Console.WriteLine("Sums:");
            Console.WriteLine("");
            Console.WriteLine("  " + sums[0].ToString().PadLeft(4)
                + "  " + sums[1].ToString().PadLeft(4) + "");
            
            indx = Ranking.partition_greedy(n, a2);

            Console.WriteLine("");
            Console.WriteLine("");
            Console.WriteLine("Data set #2 partitioned:");
            Console.WriteLine("");

            sums[0] = 0;
            sums[1] = 0;

            for (i = 0; i < n; i++)
            {
                if (indx[i] == 1)
                {
                    sums[0] = sums[0] + a2[i];
                    Console.WriteLine("  " + a2[i].ToString().PadLeft(4) + "");
                }
                else
                {
                    sums[1] = sums[1] + a2[i];
                    Console.WriteLine("  " + "    "
                        + "  " + a2[i].ToString().PadLeft(4) + "");
                }
            }

            Console.WriteLine("");
            Console.WriteLine("Sums:");
            Console.WriteLine("");
            Console.WriteLine("  " + sums[0].ToString().PadLeft(4)
                + "  " + sums[1].ToString().PadLeft(4) + "");
        }

        static void partn_enum_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PARTN_ENUM_TEST tests PARTN_ENUM.
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
            int nmax;

            Console.WriteLine("");
            Console.WriteLine("PARTN_ENUM_TEST");
            Console.WriteLine("  PARTN_ENUM enumerates partitions of N with maximum part NMAX.");
            Console.WriteLine("");
            Console.WriteLine("   NMAX:      1       2       3       4       5       6");
            Console.WriteLine("   N");

            for (n = 0; n <= 10; n++)
            {
                string cout = "  " + n.ToString().PadLeft(2) + ":  ";
                for (nmax = 1; nmax <= Math.Min(n, 6); nmax++)
                {
                    cout += "  " + Ranking.partn_enum(n, nmax).ToString().PadLeft(6);
                }

                Console.WriteLine(cout);
            }
        }

        static void partn_sf_check_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PARTN_SF_CHECK_TEST tests PARTN_SF_CHECK.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 January 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int[] a = new int[1];
            int[] a1 =  {
                6, 4, 4, 1
            }
            ;
            int[] a2 =  {
                6, 4, 4, 1
            }
            ;
            int[] a3 =  {
                6, 6, 6, -3
            }
            ;
            int[] a4 =  {
                8, 4, 2, 1
            }
            ;
            int[] a5 =  {
                1, 4, 4, 6
            }
            ;
            int[] a6 =  {
                6, 5, 4, 1
            }
            ;
            int[] a7 =  {
                6, 4, 4, 1
            }
            ;
            bool check;
            int n = 0;
            int nmax = 0;
            int npart = 0;
            int test;

            Console.WriteLine("");
            Console.WriteLine("PARTN_SF_CHECK TEST");
            Console.WriteLine("  PARTN_SF_CHECK checks a standard form partition");
            Console.WriteLine("  of N with largest entry NMAX.");

            for (test = 1; test <= 7; test++)
            {
                if (test == 1)
                {
                    n = 0;
                    nmax = 6;
                    npart = 4;
                    // a = typeMethods.i4vec_copy_new(n, a1);
                    a = a1;
                }
                else if (test == 2)
                {
                    n = 15;
                    nmax = 6;
                    npart = 0;
                    // a = typeMethods.i4vec_copy_new(n, a2);
                    a = a2;
                }
                else if (test == 3)
                {
                    n = 15;
                    nmax = 6;
                    npart = 4;
                    // a = typeMethods.i4vec_copy_new(n, a3);
                    a = a3;
                }
                else if (test == 4)
                {
                    n = 15;
                    nmax = 6;
                    npart = 4;
                    // a = typeMethods.i4vec_copy_new(n, a4);
                    a = a4;
                }
                else if (test == 5)
                {
                    n = 15;
                    nmax = 6;
                    npart = 4;
                    // a = typeMethods.i4vec_copy_new(n, a5);
                    a = a5;
                }
                else if (test == 6)
                {
                    n = 15;
                    nmax = 6;
                    npart = 4;
                    // a = typeMethods.i4vec_copy_new(n, a6);
                    a = a6;
                }
                else if (test == 7)
                {
                    n = 15;
                    nmax = 6;
                    npart = 4;
                    // a = typeMethods.i4vec_copy_new(n, a7);
                    a = a7;
                }

                Console.WriteLine("");
                Console.WriteLine("  Partition in SF form.");
                Console.WriteLine("  Partition of N = " + n + "");
                Console.WriteLine("  Maximum entry NMAX = " + nmax + "");
                Console.WriteLine("  Number of parts NPART = " + npart + "");
                typeMethods.i4vec_transpose_print(npart, a, "");
                check = Ranking.partn_sf_check(n, nmax, npart, a);
                Console.WriteLine("  Check = " + check + "");
            }
        }

        static void partn_successor_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PARTN_SUCCESSOR_TEST tests PARTN_SUCCESSOR.
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
            int[] b;
            int i;
            int n;
            int nmax;
            int npart = 0;
            int npart2 = 0;
            int rank;
            int rank_old;
            int[] t;

            Console.WriteLine("");
            Console.WriteLine("PARTN_SUCCESSOR_TEST");
            Console.WriteLine("  PARTN_SUCCESSOR lists partitions of N with maximum element NMAX:");
            Console.WriteLine("");

            n = 11;
            nmax = 4;
            t = new int[n];

            rank = -1;

            for (;;)
            {
                rank_old = rank;

                Ranking.partn_successor(n, nmax, ref npart, ref t, ref rank);

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

            //
            //  List conjugates.
            //
            Console.WriteLine("");
            Console.WriteLine("  Repeat, but list RSF conjugated partitions.");
            Console.WriteLine("");

            rank = -1;

            for (;;)
            {
                rank_old = rank;

                Ranking.partn_successor(n, nmax, ref npart, ref t, ref rank);

                if (rank <= rank_old)
                {
                    break;
                }

                b = Ranking.part_sf_conjugate(n, npart, t, ref npart2);

                typeMethods.i4vec_reverse(npart2, ref b);

                string cout = "  " + rank.ToString().PadLeft(4);
                for (i = 0; i < npart2; i++)
                {
                    cout += "  " + b[i].ToString().PadLeft(4);
                }

                Console.WriteLine(cout);
            }
        }
    }
}