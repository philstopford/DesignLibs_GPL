using System;
using Burkardt;
using Burkardt.Types;

namespace SubsetTestNS
{
    public static class DerangeTest
    {
        public static void derange_enum_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DERANGE_ENUM_TEST tests DERANGE_ENUM.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 10;

            int i;

            Console.WriteLine("");
            Console.WriteLine("DERANGE_ENUM_TEST");
            Console.WriteLine("  DERANGE_ENUM counts derangements;");
            Console.WriteLine("");
            Console.WriteLine("       N    # of derangements");
            Console.WriteLine("");

            for (i = 0; i <= N; i++)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(8)
                                  + "  " + Derange.derange_enum(i).ToString().PadLeft(8) + "");
            }
            
        }

        public static void derange_enum2_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DERANGE_ENUM2_TEST tests DERANGE_ENUM2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 10;

            int[] d = new int[N + 1];
            int i;

            Console.WriteLine("");
            Console.WriteLine("DERANGE_ENUM2_TEST");
            Console.WriteLine("  DERANGE_ENUM2 counts derangements.");
            Console.WriteLine("");
            Console.WriteLine("       N    # of derangements");
            Console.WriteLine("");

            Derange.derange_enum2(N, ref d);

            for (i = 0; i <= N; i++)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(8)
                                  + "  " + d[i] + "");
            }
        }

        public static void derange_enum3_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DERANGE_ENUM3_TEST tests DERANGE_ENUM3.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 10;

            int i;

            Console.WriteLine("");
            Console.WriteLine("DERANGE_ENUM3_TEST");
            Console.WriteLine("  DERANGE_ENUM3 counts derangements.");
            Console.WriteLine("");
            Console.WriteLine("       N    # of derangements");
            Console.WriteLine("");

            for (i = 0; i <= N; i++)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(8)
                                  + "  " + Derange.derange_enum3(i) + "");
            }
        }

        public static void derange0_back_next_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DERANGE0_BACK_NEXT_TEST tests DERANGE0_BACK_NEXT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 June 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 5;

            int[] a = new int[N];
            int i;
            bool more;
            int rank;
            Derange.DerangeBackData data = new Derange.DerangeBackData();

            Console.WriteLine("");
            Console.WriteLine("DERANGE0_BACK_NEXT_TEST");
            Console.WriteLine("  DERANGE0_BACK_NEXT generates derangements");
            Console.WriteLine("  using backtracking.");
            Console.WriteLine("");

            more = false;
            rank = 0;

            for (;;)
            {
                Derange.derange0_back_next(ref data, N, ref a, ref more);

                if (!more)
                {
                    break;
                }

                rank = rank + 1;

                string cout = rank.ToString().PadLeft(4) + "    ";
                for (i = 0; i < N; i++)
                {
                    cout += a[i].ToString().PadLeft(4) + "  ";
                }

                Console.WriteLine(cout);

            }
        }

        public static void derange0_check_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DERANGE0_CHECK_TEST tests DERANGE0_CHECK.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 June 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int[] a = new int[5];
            int[] a_test =  {
                1, 2, 3, 4, 0,
                1, 4, 2, 0, 3,
                1, 2, 3, 0, 3,
                -1, 2, 3, 4, 0,
                0, 3, 8, 1, 2
            }
            ;
            bool check;
            int i;
            int j;
            int n = 5;
            int n_test = 5;

            Console.WriteLine("");
            Console.WriteLine("DERANGE0_CHECK_TEST");
            Console.WriteLine("  DERANGE0_CHECK checks whether a vector of N objects");
            Console.WriteLine("  is a derangement of (0,...,N-1).");

            for (j = 0; j < n_test; j++)
            {
                for (i = 0; i < n; i++)
                {
                    a[i] = a_test[i + j * n];
                }

                typeMethods.i4vec_transpose_print(n, a, "  Potential derangement:");
                check = Derange.derange0_check(n, a);
                Console.WriteLine("  CHECK = " + check + "");
            }
        }

        public static void derange0_weed_next_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DERANGE0_WEED_NEXT_TEST tests DERANGE0_WEED_NEXT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 June 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int[] a;
            int i;
            int maxder;
            bool more;
            int n;
            int numder;
            int rank;

            n = 5;
            a = new int[n];
            more = false;
            maxder = 0;
            numder = 0;

            Console.WriteLine("");
            Console.WriteLine("DERANGE0_WEED_NEXT_TEST");
            Console.WriteLine("  DERANGE0_WEED_NEXT generates derangements");
            Console.WriteLine("  by generating ALL permutations, and weeding out");
            Console.WriteLine("  the ones that are not derangements.");
            Console.WriteLine("");

            rank = 0;

            for (;;)
            {
                Derange.derange0_weed_next(n, a, ref more, ref maxder, ref numder);

                rank = rank + 1;

                string cout = rank.ToString().PadLeft(4) + ":   ";
                for (i = 0; i < n; i++)
                {
                    cout += a[i].ToString().PadLeft(4) + "  ";
                }

                Console.WriteLine(cout);

                if (!more)
                {
                    break;
                }

            }
        }

    }
}