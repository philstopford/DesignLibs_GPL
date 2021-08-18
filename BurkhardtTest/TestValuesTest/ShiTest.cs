using System;
using Burkardt.Values;

namespace TestValuesTest
{
    public class ShiTest
    {
        public static void shi_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SHI_VALUES_TEST tests SHI_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    11 June 2007
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
            Console.WriteLine("SHI_VALUES_TEST:");
            Console.WriteLine("  SHI_VALUES stores values of");
            Console.WriteLine("  the hyperbolic sine integral function.");
            Console.WriteLine("");
            Console.WriteLine("      X            SHI(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Shi.shi_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString().PadLeft(12) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }
    }
}