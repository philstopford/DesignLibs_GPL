using System;
using Burkardt.Values;

namespace TestValuesTest
{
    public static class LogisticTest
    {
        public static void logistic_cdf_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LOGISTIC_CDF_VALUES_TEST tests LOGISTIC_CDF_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    13 June 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double beta = 0;
            double fx = 0;
            double mu = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("LOGISTIC_CDF_VALUES_TEST:");
            Console.WriteLine("  LOGISTIC_CDF_VALUES returns values of ");
            Console.WriteLine("  the Logistic Cumulative Density Function.");
            Console.WriteLine("");
            Console.WriteLine("     Mu      Beta         X   CDF(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Logistic.logistic_cdf_values(ref n_data, ref mu, ref beta, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + mu.ToString().PadLeft(8) + "  "
                                  + beta.ToString().PadLeft(8) + "  "
                                  + x.ToString().PadLeft(8) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }
    }
}