using System;
using Burkardt;

namespace KnapsackTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KNAPSACK_01_TEST tests the KNAPSACK_01 library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    22 August 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("KNAPSACK_01_TEST");
            Console.WriteLine("  Test the KNAPSACK_01 library.");

            test01();

            Console.WriteLine("");
            Console.WriteLine("KNAPSACK_01_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 seeks a solution of the 0/1 Knapsack problem.
            //
            //  Discussion:
            //
            //    In the 0/1 knapsack problem, a knapsack of capacity C is given,
            //    as well as N items, with the I-th item of weight W(I).
            //
            //    A selection is "acceptable" if the total weight is no greater than C.
            //
            //    It is desired to find an optimal acceptable selection, that is,
            //    an acceptable selection such that there is no acceptable selection
            //    of greater weight.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 August 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int c;
            int i;
            int n = 6;
            int[] s;
            int t;
            int[] w =
            {
                16, 17, 23, 24, 39, 40
            };

            c = 100;

            Console.WriteLine("");
            Console.WriteLine("TEST01:");
            Console.WriteLine("  Knapsack maximum capacity is " + c + "");
            Console.WriteLine("  Come as close as possible to filling the knapsack.");

            s = Knapsack.knapsack_01(n, w, c);

            Console.WriteLine("");
            Console.WriteLine("   # 0/1  Weight");
            Console.WriteLine("");
            for (i = 0; i < n; i++)
            {
                Console.WriteLine(i.ToString().PadLeft(4) + "  "
                                                          + s[i].ToString().PadLeft(1) + "  "
                                                          + w[i].ToString().PadLeft(4) + "");
            }

            t = 0;
            for (i = 0; i < n; i++)
            {
                t = t + s[i] * w[i];
            }

            Console.WriteLine("");
            Console.WriteLine("  Total:   " + t + "");
        }
    }
}