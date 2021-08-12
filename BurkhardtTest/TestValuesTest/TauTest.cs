using System;
using Burkardt.TestValues;

namespace TestValuesTest
{
    public class TauTest
    {
        public static void tau_values_test()
            //****************************************************************************80
            //
            //  Purpose: 
            //
            //    TAU_VALUES_TEST tests TAU_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 February 2007
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
            Console.WriteLine("TAU_VALUES_TEST:");
            Console.WriteLine("  TAU_VALUES returns values of");
            Console.WriteLine("  the TAU function.");
            Console.WriteLine("");
            Console.WriteLine("     N         TAU(N)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Tau.tau_values(ref n_data, ref n, ref fn);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + n.ToString().PadLeft(6) + "  "
                                  + fn.ToString().PadLeft(12) + "");
            }
        }
    }
}