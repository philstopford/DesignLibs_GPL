using System;
using Burkardt;
using Burkardt.RankingNS;
using Burkardt.Types;

namespace ComboTest
{
    partial class Program
    {
        static void cycle_check_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CYCLE_CHECK_TEST tests CYCLE_CHECK.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 January 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            bool check;
            int i;
            int[] index = new int[1];
            int[] index1 =  {
                1, 4, 3
            }
            ;
            int[] index3 =  {
                1, 4, 2
            }
            ;
            int[] index4 =  {
                1, 4, 3
            }
            ;
            int[] index5 =  {
                1, 4, 3
            }
            ;
            int[] index6 =  {
                1, 4, 3
            }
            ;
            int j;
            int jlo;
            int n = 0;
            int ncycle = 0;
            int[] t = new int[1];
            int[] t1 =  {
                5, 1, 3, 8, 6, 2, 4, 7
            }
            ;
            int[] t3 =  {
                5, 1, 3, 8, 6, 2, 4, 7
            }
            ;
            int[] t4 =  {
                5, 1, 3, 12, 6, 2, 4, 7
            }
            ;
            int[] t5 =  {
                5, 1, 3, 8, 5, 2, 4, 7
            }
            ;
            int[] t6 =  {
                5, 1, 3, 8, 6, 2, 4, 7
            }
            ;
            int test;

            Console.WriteLine("");
            Console.WriteLine("CYCLE_CHECK TEST");
            Console.WriteLine("  CYCLE_CHECK checks a permutation in cycle form.");

            for (test = 1; test <= 6; test++)
            {
                if (test == 1)
                {
                    n = 0;
                    ncycle = 3;
                    t = typeMethods.i4vec_copy_new(8, t1);
                    index = typeMethods.i4vec_copy_new(ncycle, index1);
                }
                else if (test == 2)
                {
                    n = 8;
                    ncycle = 0;
                    t = typeMethods.i4vec_copy_new(n, t1);
                    index = typeMethods.i4vec_copy_new(3, index1);
                }
                else if (test == 3)
                {
                    n = 8;
                    ncycle = 3;
                    t = typeMethods.i4vec_copy_new(n, t3);
                    index = typeMethods.i4vec_copy_new(ncycle, index3);
                }
                else if (test == 4)
                {
                    n = 8;
                    ncycle = 3;
                    t = typeMethods.i4vec_copy_new(n, t4);
                    index = typeMethods.i4vec_copy_new(ncycle, index4);
                }
                else if (test == 5)
                {
                    n = 8;
                    ncycle = 3;
                    t = typeMethods.i4vec_copy_new(n, t5);
                    index = typeMethods.i4vec_copy_new(ncycle, index5);
                }
                else if (test == 6)
                {
                    n = 8;
                    ncycle = 3;
                    t = typeMethods.i4vec_copy_new(n, t6);
                    index = typeMethods.i4vec_copy_new(ncycle, index6);
                }

                Console.WriteLine("");
                Console.WriteLine("  Permutation in cycle form:");
                Console.WriteLine("  Number of cycles is " + ncycle + "");
                Console.WriteLine("");
                jlo = 0;
                for (i = 0; i < ncycle; i++)
                {
                    string cout = "    ";
                    for (j = jlo; j < jlo + index[i]; j++)
                    {
                        cout += "  " + t[j];
                    }

                    Console.WriteLine(cout);
                    jlo = jlo + index[i];
                }

                check = Ranking.cycle_check(n, ncycle, t, index);
                Console.WriteLine("  Check = " + check + "");
            }
        }

        static void cycle_to_perm_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CYCLE_TO_PERM_TEST tests CYCLE_TO_PERM.
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
            int j;
            int jlo;
            int[] index =  {
                5, 1, 1
            }
            ;
            int n = 7;
            int ncycle;
            int[] p;
            int[] t =  {
                4, 2, 5, 3, 1, 6, 7
            }
            ;

            Console.WriteLine("");
            Console.WriteLine("CYCLE_TO_PERM_TEST");
            Console.WriteLine("  CYCLE_TO_PERM converts a permutation from");
            Console.WriteLine("  cycle to array form;");

            n = 7;
            ncycle = 3;

            Console.WriteLine("");
            Console.WriteLine("  Cycle form:");
            Console.WriteLine("  Number of cycles is " + ncycle + "");
            Console.WriteLine("");
            jlo = 0;
            for (i = 1; i <= ncycle; i++)
            {
                string cout = "";
                for (j = jlo + 1; j <= jlo + index[i - 1]; j++)
                {
                    cout += "  " + t[j - 1].ToString().PadLeft(4);
                }

                Console.WriteLine(cout);
                jlo = jlo + index[i - 1];
            }

            p = Ranking.cycle_to_perm(n, ncycle, t, index);

            Permutation.perm_print(n, p, "  Corresponding permutation:");
        }
    }
}