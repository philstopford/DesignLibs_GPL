﻿using System;
using Burkardt.RankingNS;
using Burkardt.Types;

namespace ComboTest;

internal static partial class Program
{
    private static void partition_greedy_test()

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
        const int N = 10;

        int[] a1 =  {
                2, 10, 3, 8, 5, 7, 9, 5, 3, 2
            }
            ;
        int[] a2 =  {
                771, 121, 281, 854, 885, 734, 486, 1003, 83, 62
            }
            ;
        int i;
        int[] sums = new int[2];

        Console.WriteLine("");
        Console.WriteLine("PARTITION_GREEDY_TEST");
        Console.WriteLine("  PARTITION_GREEDY partitions an integer vector into");
        Console.WriteLine("  two subsets with nearly equal sum.");
        Console.WriteLine("");

        int[] indx = Ranking.partition_greedy(N, ref a1);

        Console.WriteLine("");
        Console.WriteLine("");
        Console.WriteLine("Data set #1 partitioned:");
        Console.WriteLine("");
        sums[0] = 0;
        sums[1] = 0;

        for (i = 0; i < N; i++)
        {
            switch (indx[i])
            {
                case 1:
                    sums[0] += a1[i];
                    Console.WriteLine("  " + a1[i].ToString().PadLeft(4) + "");
                    break;
                default:
                    sums[1] += a1[i];
                    Console.WriteLine("  " + "    "
                                           + "  " + a1[i].ToString().PadLeft(4) + "");
                    break;
            }
        }

        Console.WriteLine("");
        Console.WriteLine("Sums:");
        Console.WriteLine("");
        Console.WriteLine("  " + sums[0].ToString().PadLeft(4)
                               + "  " + sums[1].ToString().PadLeft(4) + "");
            
        indx = Ranking.partition_greedy(N, ref a2);

        Console.WriteLine("");
        Console.WriteLine("");
        Console.WriteLine("Data set #2 partitioned:");
        Console.WriteLine("");

        sums[0] = 0;
        sums[1] = 0;

        for (i = 0; i < N; i++)
        {
            switch (indx[i])
            {
                case 1:
                    sums[0] += a2[i];
                    Console.WriteLine("  " + a2[i].ToString().PadLeft(4) + "");
                    break;
                default:
                    sums[1] += a2[i];
                    Console.WriteLine("  " + "    "
                                           + "  " + a2[i].ToString().PadLeft(4) + "");
                    break;
            }
        }

        Console.WriteLine("");
        Console.WriteLine("Sums:");
        Console.WriteLine("");
        Console.WriteLine("  " + sums[0].ToString().PadLeft(4)
                               + "  " + sums[1].ToString().PadLeft(4) + "");
    }

    private static void partn_enum_test()

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

        Console.WriteLine("");
        Console.WriteLine("PARTN_ENUM_TEST");
        Console.WriteLine("  PARTN_ENUM enumerates partitions of N with maximum part NMAX.");
        Console.WriteLine("");
        Console.WriteLine("   NMAX:      1       2       3       4       5       6");
        Console.WriteLine("   N");

        for (n = 0; n <= 10; n++)
        {
            string cout = "  " + n.ToString().PadLeft(2) + ":  ";
            int nmax;
            for (nmax = 1; nmax <= Math.Min(n, 6); nmax++)
            {
                cout += "  " + Ranking.partn_enum(n, nmax).ToString().PadLeft(6);
            }

            Console.WriteLine(cout);
        }
    }

    private static void partn_sf_check_test()

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
            switch (test)
            {
                case 1:
                    n = 0;
                    nmax = 6;
                    npart = 4;
                    // a = typeMethods.i4vec_copy_new(n, a1);
                    a = a1;
                    break;
                case 2:
                    n = 15;
                    nmax = 6;
                    npart = 0;
                    // a = typeMethods.i4vec_copy_new(n, a2);
                    a = a2;
                    break;
                case 3:
                    n = 15;
                    nmax = 6;
                    npart = 4;
                    // a = typeMethods.i4vec_copy_new(n, a3);
                    a = a3;
                    break;
                case 4:
                    n = 15;
                    nmax = 6;
                    npart = 4;
                    // a = typeMethods.i4vec_copy_new(n, a4);
                    a = a4;
                    break;
                case 5:
                    n = 15;
                    nmax = 6;
                    npart = 4;
                    // a = typeMethods.i4vec_copy_new(n, a5);
                    a = a5;
                    break;
                case 6:
                    n = 15;
                    nmax = 6;
                    npart = 4;
                    // a = typeMethods.i4vec_copy_new(n, a6);
                    a = a6;
                    break;
                case 7:
                    n = 15;
                    nmax = 6;
                    npart = 4;
                    // a = typeMethods.i4vec_copy_new(n, a7);
                    a = a7;
                    break;
            }

            Console.WriteLine("");
            Console.WriteLine("  Partition in SF form.");
            Console.WriteLine("  Partition of N = " + n + "");
            Console.WriteLine("  Maximum entry NMAX = " + nmax + "");
            Console.WriteLine("  Number of parts NPART = " + npart + "");
            typeMethods.i4vec_transpose_print(npart, a, "");
            bool check = Ranking.partn_sf_check(n, nmax, npart, a);
            Console.WriteLine("  Check = " + check + "");
        }
    }

    private static void partn_successor_test()

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
        int i;
        int npart = 0;
        int npart2 = 0;
        int rank_old;

        Console.WriteLine("");
        Console.WriteLine("PARTN_SUCCESSOR_TEST");
        Console.WriteLine("  PARTN_SUCCESSOR lists partitions of N with maximum element NMAX:");
        Console.WriteLine("");

        const int n = 11;
        const int nmax = 4;
        int[] t = new int[n];

        int rank = -1;

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

            Console.WriteLine(cout);
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

            int[] b = Ranking.part_sf_conjugate(n, npart, t, ref npart2);

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