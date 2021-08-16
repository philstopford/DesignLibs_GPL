using System;
using Burkardt.Function;

namespace PolPakTest
{
    public static class omegaTest
    {
        public static void omega_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    OMEGA_TEST tests OMEGA.
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
            Console.WriteLine("OMEGA_TEST");
            Console.WriteLine("  OMEGA computes the OMEGA function.");
            Console.WriteLine("");
            Console.WriteLine("          N   Exact   OMEGA(N)");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                Burkardt.TestValues.Omega.omega_values(ref n_data, ref n, ref c);

                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + n.ToString().PadLeft(12) + "  "
                                  + c.ToString().PadLeft(10) + "  "
                                  + Omega.omega(n).ToString().PadLeft(10) + "");

            }

        }

    }
}