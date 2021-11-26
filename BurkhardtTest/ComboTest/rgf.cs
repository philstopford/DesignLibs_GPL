using System;
using Burkardt.RankingNS;
using Burkardt.Types;

namespace ComboTest;

internal static partial class Program
{
    private static void rgf_check_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RGF_CHECK_TEST tests RGF_CHECK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] f = new int[1];
        int[] f1 =  {
                -1
            }
            ;
        int[] f2 =  {
                0, 1, 2, 3, 4, 5, 6
            }
            ;
        int[] f3 =  {
                1, 3, 5, 8, 9, 10, 12
            }
            ;
        int[] f4 =  {
                1, 2, 3, 1, 4, 5, 4
            }
            ;
        int m = 0;
        int test;

        Console.WriteLine("");
        Console.WriteLine("RGF_CHECK TEST");
        Console.WriteLine("  RGF_CHECK checks a restricted growth function.");

        for (test = 1; test <= 4; test++)
        {
            switch (test)
            {
                case 1:
                    m = 1;
                    //f = typeMethods.i4vec_copy_new(m, f1);
                    f = f1;
                    break;
                case 2:
                    m = 7;
                    // f = typeMethods.i4vec_copy_new(m, f2);
                    f = f2;
                    break;
                case 3:
                    m = 7;
                    // f = typeMethods.i4vec_copy_new(m, f3);
                    f = f3;
                    break;
                case 4:
                    m = 7;
                    // f = typeMethods.i4vec_copy_new(m, f4);
                    f = f4;
                    break;
            }

            Console.WriteLine("");
            typeMethods.i4vec_transpose_print(m, f, "  RGF:");
            bool check = Ranking.rgf_check(m, f);
            Console.WriteLine("  CHECK = " + check + "");
        }
    }

    private static void rgf_enum_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RGF_ENUM_TEST tests RGF_ENUM.
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
        Console.WriteLine("RGF_ENUM_TEST");
        Console.WriteLine("  RGF_ENUM enumerates restricted growth functions.");
        Console.WriteLine("");
        Console.WriteLine("   N   #");
        Console.WriteLine("");

        for (n = 0; n <= 10; n++)
        {
            int rgf_num = Ranking.rgf_enum(n);
            Console.WriteLine("  " + n.ToString().PadLeft(2)
                                   + "  " + rgf_num.ToString().PadLeft(6) + "");
        }
    }

    private static void rgf_g_table_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RGF_G_TABLE_TEST tests RGF_G_TABLE.
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
        int i;
        const int m = 6;

        Console.WriteLine("");
        Console.WriteLine("RGF_G_TABLE_TEST");
        Console.WriteLine("  RGF_G_TABLE tabulates generalized restricted");
        Console.WriteLine("  growth functions.");
        Console.WriteLine("");

        int[] d = Ranking.rgf_g_table(m);

        for (i = 0; i <= m; i++)
        {
            string cout = "";
            int j;
            for (j = 0; j <= m - i; j++)
            {
                cout += "  " + d[i + j * (m + 1)].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);
        }
    }

    private static void rgf_rank_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RGF_RANK_TEST tests RGF_RANK.
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
        int[] f =  {
                1, 2, 1, 3
            }
            ;
        const int m = 4;

        Console.WriteLine("");
        Console.WriteLine("RGF_RANK_TEST");
        Console.WriteLine("  RGF_RANK ranks restricted growth functions:");

        typeMethods.i4vec_transpose_print(m, f, "  Element to be ranked:");

        int rank = Ranking.rgf_rank(m, f);

        Console.WriteLine("");
        Console.WriteLine("  Rank is computed as " + rank + "");
    }

    private static void rgf_successor_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RGF_SUCCESSOR_TEST tests RGF_SUCCESSOR.
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
        const int m = 4;

        Console.WriteLine("");
        Console.WriteLine("RGF_SUCCESSOR_TEST");
        Console.WriteLine("  RGF_SUCCESSOR lists restricted growth functions:");

        int[] f = new int[m];

        int rank = -1;

        for (;;)
        {
            int rank_old = rank;

            Ranking.rgf_successor(m, ref f, ref rank);

            if (rank <= rank_old)
            {
                break;
            }

            string cout = "  " + rank.ToString().PadLeft(4);
            int i;
            for (i = 0; i < m; i++)
            {
                cout += "  " + f[i].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);
        }
    }

    private static void rgf_to_setpart_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RGF_TO_SETPART_TEST tests RGF_TO_SETPART.
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
        int[] f =  {
                1, 1, 1, 1, 1, 2, 1, 3
            }
            ;
        const int m = 8;
        int nsub = 0;

        Console.WriteLine("");
        Console.WriteLine("RGF_TO_SETPART_TEST");
        Console.WriteLine("  RGF_TO_SETPART converts a restricted growth function");
        Console.WriteLine("  to a set partition;");

        typeMethods.i4vec_transpose_print(m, f, "  Restricted growth function:");
        //
        //  Convert the RGF to a set partition.
        //
        int[] s = new int[m];
        int[] index = new int[m];

        Ranking.rgf_to_setpart(m, f, ref nsub, ref s, ref index);

        Console.WriteLine("");
        Console.WriteLine("  Corresponding set partition");
        Console.WriteLine("");
        int jlo = 1;
        for (i = 1; i <= nsub; i++)
        {
            string cout = "";
            int j;
            for (j = jlo; j <= index[i - 1]; j++)
            {
                cout += "  " + s[j - 1].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);
            jlo = index[i - 1] + 1;
        }
    }

    private static void rgf_unrank_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RGF_UNRANK_TEST tests RGF_UNRANK.
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
        Console.WriteLine("");
        Console.WriteLine("RGF_UNRANK_TEST");
        Console.WriteLine("  RGF_UNRANK unranks restricted growth functions:");

        const int rank = 7;
        int m;
        m = 4;

        int[] f = Ranking.rgf_unrank(rank, m);

        Console.WriteLine("");
        Console.WriteLine("  The element of rank " + rank + "");
        Console.WriteLine("");
        typeMethods.i4vec_transpose_print(m, f, "");
    }
}