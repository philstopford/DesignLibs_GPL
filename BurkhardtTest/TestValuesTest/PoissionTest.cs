using System;
using Burkardt.Values;

namespace TestValuesTest
{
    public class PoissionTest
    {
        public static void poisson_cdf_values_test()
            //****************************************************************************80
            //
            //  Purpose: 
            //
            //    POISSON_CDF_VALUES_TEST tests POISSON_CDF_VALUES.
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
            double a = 0;
            double fx = 0;
            int n_data;
            int x = 0;
            Console.WriteLine("");
            Console.WriteLine("POISSON_CDF_VALUES_TEST:");
            Console.WriteLine("  POISSON_CDF_VALUES returns values of");
            Console.WriteLine("  the Poisson Cumulative Density Function.");
            Console.WriteLine("");
            Console.WriteLine("      A     X       CDF(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Poisson.poisson_cdf_values(ref n_data, ref a, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + a.ToString().PadLeft(8) + a + "  "
                                  + x.ToString().PadLeft(4) + x + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }
    }
}