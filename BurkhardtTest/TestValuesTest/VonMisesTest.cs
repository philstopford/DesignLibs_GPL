using System;
using TestValues;

namespace TestValuesTest
{
    public class VonMisesTest
    {
        public static void von_mises_cdf_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    VON_MISES_CDF_VALUES_TEST tests VON_MISES_CDF_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 June 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a = 0;
            double b = 0;
            double fx = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("VON_MISES_CDF_VALUES_TEST:");
            Console.WriteLine("  VON_MISES_CDF_VALUES stores values of");
            Console.WriteLine("  the von Mises CDF.");
            Console.WriteLine("");
            Console.WriteLine("      A            B            X            CDF(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                VonMises.von_mises_cdf_values(ref n_data, ref a, ref b, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + a.ToString().PadLeft(12) + "  "
                                  + b.ToString().PadLeft(12) + "  "
                                  + x.ToString().PadLeft(12) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }
    }
}