using System;
using Burkardt.Values;

namespace TestValuesTest
{
    public class StromgenTest
    {
        public static void stromgen_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    STROMGEN_VALUES_TEST tests STROMGEN_VALUES.
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
            double fx = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("STROMGEN_VALUES_TEST:");
            Console.WriteLine("  STROMGEN_VALUES stores values of ");
            Console.WriteLine("  the Stromgen function.");
            Console.WriteLine("");
            Console.WriteLine("                X                     FX");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Stromgen.stromgen_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString("0.################").PadLeft(24) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }
    }
}