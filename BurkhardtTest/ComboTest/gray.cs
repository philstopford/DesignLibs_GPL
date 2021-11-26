﻿using System;
using Burkardt.RankingNS;
using Burkardt.Types;

namespace ComboTest;

internal static partial class Program
{
    private static void gray_code_check_test()

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
            const int n = 5;

            switch (test)
            {
                case 1:
                case 2:
                case 3:
                    t = typeMethods.i4vec_copy_new(n, t3);
                    break;
            }

            bool check = typeMethods.gray_code_check(n, t3);
            string cout = "      " + check.ToString().PadLeft(1)
                                   + "  " + n.ToString().PadLeft(2)
                                   + ":  ";
            int i;
            for (i = 0; i < n; i++)
            {
                cout += "  " + t[i].ToString().PadLeft(2);
            }

            Console.WriteLine(cout);
        }
    }

    private static void gray_code_enum_test()

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

        Console.WriteLine("");
        Console.WriteLine("GRAY_CODE_ENUM_TEST");
        Console.WriteLine("  GRAY_CODE_ENUM enumerates Gray codes with N elements.");
        Console.WriteLine("");
        Console.WriteLine("   N   Enum(N)");
        Console.WriteLine("");
        for (n = 0; n <= 10; n++)
        {
            int ngray = typeMethods.gray_code_enum(n);
            Console.WriteLine("  " + n.ToString().PadLeft(2)
                                   + "  " + ngray.ToString().PadLeft(6) + "");
        }
    }

    private static void gray_code_rank_test()

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
        const int n = 5;
        int[] t =  {
                1, 1, 0, 0, 0
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("GRAY_CODE_RANK_TEST");
        Console.WriteLine("  GRAY_CODE_RANK ranks a Gray code.");

        int rank = Ranking.gray_code_rank(n, t);

        typeMethods.i4vec_transpose_print(n, t, "  Element to be ranked:");

        Console.WriteLine("");
        Console.WriteLine("  Computed rank: " + rank + "");
    }

    private static void gray_code_successor_test()

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
        Console.WriteLine("");
        Console.WriteLine("GRAY_CODE_SUCCESSOR_TEST");
        Console.WriteLine("  GRAY_CODE_SUCCESSOR lists Gray codes one by one.");

        const int n = 5;
        int[] t = new int[n];
        int rank = -1;

        for (;;)
        {
            int rank_old = rank;

            typeMethods.gray_code_successor(n, ref t, ref rank);

            if (rank <= rank_old)
            {
                break;
            }

            string cout = "  " + rank.ToString().PadLeft(4);
            int i;
            for (i = 0; i < n; i++)
            {
                cout += "  " + t[i].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);
        }
    }

    private static void gray_code_unrank_test()

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

        Console.WriteLine("");
        Console.WriteLine("GRAY_CODE_UNRANK_TEST");
        Console.WriteLine("  GRAY_CODE_UNRANK unranks a Gray code.");

        const int n = 5;
        int ngray = typeMethods.gray_code_enum(n);
        int rank = ngray / 2;

        int[] t = Ranking. gray_code_unrank(rank, n);

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