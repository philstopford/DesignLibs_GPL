using System;
using Burkardt.Sequence;

namespace PolPakTest
{
    public static class fibonacciTest
    {
        public static void fibonacci_direct_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FIBONACCI_DIRECT_TEST tests FIBONACCI_DIRECT.
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
            int f;
            int i;
            int n = 20;

            Console.WriteLine("");
            Console.WriteLine("FIBONACCI_DIRECT_TEST");
            Console.WriteLine("  FIBONACCI_DIRECT evalutes a Fibonacci number directly.");
            Console.WriteLine("");

            for (i = 1; i <= n; i++)
            {
                f = Fibonacci.fibonacci_direct(i);

                Console.WriteLine("  "
                                  + i.ToString().PadLeft(6) + "  "
                                  + f.ToString().PadLeft(10) + "");
            }

        }

        public static void fibonacci_floor_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FIBONACCI_FLOOR_TEST tests FIBONACCI_FLOOR.
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
            int f = 0;
            int i = 0;
            int n = 0;

            Console.WriteLine("");
            Console.WriteLine("FIBONACCI_FLOOR_TEST");
            Console.WriteLine("  FIBONACCI_FLOOR computes the largest Fibonacci number");
            Console.WriteLine("  less than or equal to a given positive integer.");
            Console.WriteLine("");
            Console.WriteLine("     N  Fibonacci  Index");
            Console.WriteLine("");

            for (n = 1; n <= 20; n++)
            {
                Fibonacci.fibonacci_floor(n, ref f, ref i);

                Console.WriteLine("  "
                                  + n.ToString().PadLeft(6) + "  "
                                  + f.ToString().PadLeft(6) + "  "
                                  + i.ToString().PadLeft(6) + "");
            }

        }

        public static void fibonacci_recursive_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FIBONACCI_RECURSIVE_TEST tests FIBONACCI_RECURSIVE.
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
            int N = 20;

            int[] f = new int[N];
            int i;

            Console.WriteLine("");
            Console.WriteLine("FIBONACCI_RECURSIVE_TEST");
            Console.WriteLine("  FIBONACCI_RECURSIVE computes the Fibonacci sequence.");
            Console.WriteLine("");

            Fibonacci.fibonacci_recursive(N, ref f);

            for (i = 1; i <= N; i++)
            {
                Console.WriteLine("  "
                                  + i.ToString().PadLeft(6) + "  "
                                  + f[i - 1].ToString().PadLeft(10) + "");
            }

        }

    }
}