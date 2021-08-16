using System;
using Burkardt.Function;

namespace PolPakTest
{
    public static class zeckendorfTest
    {
        public static void zeckendorf_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ZECKENDORF_TEST tests ZECKENDORF.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 May 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int M_MAX = 20;

            int[] i_list = new int[M_MAX];
            int j = 0;
            int[] f_list = new int[M_MAX];
            int m = 0;
            int n = 0;

            Console.WriteLine("");
            Console.WriteLine("ZECKENDORF_TEST");
            Console.WriteLine("  ZECKENDORF computes the Zeckendorf decomposition of");
            Console.WriteLine("  an integer N into nonconsecutive Fibonacci numbers.");
            Console.WriteLine("");
            Console.WriteLine("   N Sum M Parts");
            Console.WriteLine("");

            for (n = 1; n <= 100; n++)
            {
                Zeckendorf.zeckendorf(n, M_MAX, ref m, ref i_list, ref f_list);

                string cout = n.ToString().PadLeft(4) + "  ";
                for (j = 0; j < m; j++)
                {
                    cout += f_list[j].ToString().PadLeft(4) + "  ";
                }

                Console.WriteLine(cout);

            }

        }
    }
}