using System;
using Burkardt.TestValues;

namespace TestValuesTest
{
    public class SigmaTest
    {
        public static void sigma_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIGMA_VALUES_TEST tests SIGMA_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 February 2007
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
            Console.WriteLine("SIGMA_VALUES_TEST:");
            Console.WriteLine("  SIGMA_VALUES returns values of");
            Console.WriteLine("  the SIGMA function.");
            Console.WriteLine("");
            Console.WriteLine("       N         SIGMA(N)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Sigma.sigma_values(ref n_data, ref n, ref fn);
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