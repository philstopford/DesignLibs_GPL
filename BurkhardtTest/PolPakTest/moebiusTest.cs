using System;
using Burkardt.Function;

namespace PolPakTest
{
    public static class moebiusTest
    {
        public static void moebius_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MOEBIUS_TEST tests MOEBIUS.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 June 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int c = 0;
            int n = 0;
            int n_data;

            Console.WriteLine("");
            Console.WriteLine("MOEBIUS_TEST");
            Console.WriteLine("  MOEBIUS computes the Moebius function.");
            Console.WriteLine("");
            Console.WriteLine("      N   Exact   MOEBIUS(N)");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                Burkardt.TestValues.Moebius.moebius_values(ref n_data, ref n, ref c);

                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + n.ToString().PadLeft(8) + "  "
                                  + c.ToString().PadLeft(10) + "  "
                                  + Moebius.moebius(n).ToString().PadLeft(10) + "");
            }

        }

    }
}