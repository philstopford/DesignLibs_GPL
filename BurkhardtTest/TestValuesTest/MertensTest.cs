using System;
using TestValues;

namespace TestValuesTest
{
    public static class MertensTest
    {
        public static void mertens_values_test()
            //****************************************************************************80
            //
            //  Purpose: 
            //
            //    MERTENS_VALUES_TEST tests MERTENS_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 October 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int fn = 0;
            int n = 0;
            int n_data;
            Console.WriteLine("");
            Console.WriteLine("MERTENS_VALUES_TEST:");
            Console.WriteLine("  MERTENS_VALUES returns values of");
            Console.WriteLine("  the Mertens function.");
            Console.WriteLine("");
            Console.WriteLine("     N         MERTENS(N)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Mertens.mertens_values(ref n_data, ref n, ref fn);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + n.ToString().PadLeft(8) + "  "
                                  + fn.ToString().PadLeft(12) + "");
            }
        }
    }
}