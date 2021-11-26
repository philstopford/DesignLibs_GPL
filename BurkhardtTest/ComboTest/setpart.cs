using System;
using Burkardt.RankingNS;
using Burkardt.Types;

namespace ComboTest;

internal static partial class Program
{
    private static void setpart_check_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SETPART_CHECK_TEST tests SETPART_CHECK.
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
        int[] index = new int[1];
        int[] index1 =  {
                2, 5, 8
            }
            ;
        int[] index2 =  {
                2, 5, 8
            }
            ;
        int[] index3 =  {
                2, 8, 5
            }
            ;
        int[] index4 =  {
                2, 5, 8
            }
            ;
        int[] index5 =  {
                2, 5, 8
            }
            ;
        int[] index6 =  {
                2, 5, 8
            }
            ;
        int m = 0;
        int nsub = 0;
        int[] s = new int[1];
        int[] s1 =  {
                3, 6, 1, 4, 7, 2, 5, 8
            }
            ;
        int[] s2 =  {
                3, 6, 1, 4, 7, 2, 5, 8
            }
            ;
        int[] s3 =  {
                3, 6, 1, 4, 7, 2, 5, 8
            }
            ;
        int[] s4 =  {
                3, 6, 1, 4, 9, 2, 5, 8
            }
            ;
        int[] s5 =  {
                3, 6, 1, 4, 6, 2, 5, 8
            }
            ;
        int[] s6 =  {
                3, 6, 1, 4, 7, 2, 5, 8
            }
            ;
        int test;

        Console.WriteLine("");
        Console.WriteLine("SETPART_CHECK TEST");
        Console.WriteLine("  SETPART_CHECK checks a set partition.");

        for (test = 1; test <= 6; test++)
        {
            switch (test)
            {
                case 1:
                    m = 0;
                    nsub = 3;
                    s = typeMethods.i4vec_copy_new(8, s1);
                    index = typeMethods.i4vec_copy_new(nsub, index1);
                    break;
                case 2:
                    m = 8;
                    nsub = 0;
                    s = typeMethods.i4vec_copy_new(m, s2);
                    index = typeMethods.i4vec_copy_new(3, index2);
                    break;
                case 3:
                    m = 8;
                    nsub = 3;
                    s = typeMethods.i4vec_copy_new(m, s3);
                    index = typeMethods.i4vec_copy_new(nsub, index3);
                    break;
                case 4:
                    m = 8;
                    nsub = 3;
                    s = typeMethods.i4vec_copy_new(m, s4);
                    index = typeMethods.i4vec_copy_new(nsub, index4);
                    break;
                case 5:
                    m = 8;
                    nsub = 3;
                    s = typeMethods.i4vec_copy_new(m, s5);
                    index = typeMethods.i4vec_copy_new(nsub, index5);
                    break;
                case 6:
                    m = 8;
                    nsub = 3;
                    s = typeMethods.i4vec_copy_new(m, s6);
                    index = typeMethods.i4vec_copy_new(nsub, index6);
                    break;
            }

            Console.WriteLine("");
            Console.WriteLine("  The set partition:");
            Console.WriteLine("  M = " + m + "");
            Console.WriteLine("  NSUB = " + nsub + "");
            Console.WriteLine("");
            int jlo = 0;
            int i;
            for (i = 0; i < nsub; i++)
            {
                string cout = "";
                int j;
                for (j = jlo; j <= index[i] - 1; j++)
                {
                    cout += "  " + s[j];
                }

                Console.WriteLine(cout);
                jlo = index[i];
            }

            bool check = Ranking.setpart_check(m, nsub, s, index);
            Console.WriteLine("  CHECK = " + check + "");
        }
    }

    private static void setpart_enum_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SETPART_ENUM_TEST tests SETPART_ENUM.
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
        int n;

        Console.WriteLine("");
        Console.WriteLine("SETPART_ENUM_TEST");
        Console.WriteLine("  SETPART_ENUM enumerates set partitions.");
        Console.WriteLine("");
        //
        //  Enumerate.
        //
        for (n = 1; n <= 6; n++)
        {
            int npart = Ranking.setpart_enum(n);
            Console.WriteLine("  " + n.ToString().PadLeft(4)
                                   + "  " + npart.ToString().PadLeft(4) + "");
        }
    }

    private static void setpart_to_rgf_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SETPART_TO_RGF_TEST tests SETPART_TO_RGF.
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
        int[] index =  {
                6, 7, 8
            }
            ;
        const int m = 8;
        const int nsub = 3;
        int[] s =  {
                1, 2, 3, 4, 5, 7, 6, 8
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("SETPART_TO_RGF_TEST");
        Console.WriteLine("  SETPART_TO_RGF converts a set partition to a");
        Console.WriteLine("  restricted growth function.");

        Console.WriteLine("");
        Console.WriteLine("  Set partition");
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

        //
        //  Convert the set partition to an RGF.
        //
        int[] f = Ranking.setpart_to_rgf(m, nsub, s, index);

        typeMethods.i4vec_transpose_print(m, f, "  Corresponding RGF:");
    }
}