using System;
using Burkardt;

namespace SubsetTestNS
{
    public static class YoungTableauTest
    {
        public static void ytb_enum_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    YTB_ENUM_TEST tests YTB_ENUM.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    11 October 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int n;

            Console.WriteLine("");
            Console.WriteLine("YTB_ENUM_TEST");
            Console.WriteLine("  YTB_ENUM counts Young tableau.");
            Console.WriteLine("");
            Console.WriteLine("   N  YTB_ENUM(N)");
            Console.WriteLine("");

            for (n = 0; n <= 10; n++)
            {
                Console.WriteLine(n.ToString().PadLeft(4) + "  "
                                  + YoungTableau.ytb_enum(n).ToString().PadLeft(10) + "");
            }
        }

        public static void ytb_next_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    YTB_NEXT_TEST tests YTB_NEXT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    11 October 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 6;

            int i;
            int[] a = new int[N];
            int[] lambda =  {
                3, 2, 1, 0, 0, 0
            }
            ;
            bool more;

            for (i = 0; i < N; i++)
            {
                a[i] = 0;
            }

            Console.WriteLine("");
            Console.WriteLine("YTB_NEXT_TEST");
            Console.WriteLine("  YTB_NEXT generates Young tableaus.");
            Console.WriteLine("");

            more = false;

            i = 0;

            for (;;)
            {
                i = i + 1;
                YoungTableau.ytb_next(N, ref lambda, ref a, ref more);

                YoungTableau.ytb_print(N, a, " ");

                if (!more || 100 < i)
                {
                    break;
                }

            }
        }

        public static void ytb_random_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    YTB_RANDOM_TEST tests YTB_RANDOM.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    11 October 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 6;

            int[] a = new int[N];
            int i;
            int[] lambda =  {
                3, 2, 1, 0, 0, 0
            }
            ;
            int seed;

            Console.WriteLine("");
            Console.WriteLine("YTB_RANDOM_TEST");
            Console.WriteLine("  YTB_RANDOM generates a random Young tableau");

            seed = 123456789;

            for (i = 1; i <= 5; i++)
            {
                YoungTableau.ytb_random(N, lambda, ref seed, ref a);

                YoungTableau.ytb_print(N, a, " ");

            }
        }

    }
}