using System;
using Burkardt.TestValues;

namespace TestValuesTest
{
    public static class CollatzTest
    {
        public static void collatz_count_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    COLLATZ_COUNT_VALUES_TEST tests COLLATZ_COUNT_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 March 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int count = 0;
            int n = 0;
            int n_data;
            Console.WriteLine("");
            Console.WriteLine("COLLATZ_COUNT_VALUES_TEST:");
            Console.WriteLine("  COLLATZ_COUNT_VALUES returns values of");
            Console.WriteLine("  the length of the Collatz sequence that");
            Console.WriteLine("  starts at N.");
            Console.WriteLine("");
            Console.WriteLine("         N      COLLATZ_COUNT(N)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Collatz.collatz_count_values(ref n_data, ref n, ref count);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  " + n.ToString().PadLeft(8)
                                       + "  " + count.ToString().PadLeft(12) + "");
            }
        }
    }
}