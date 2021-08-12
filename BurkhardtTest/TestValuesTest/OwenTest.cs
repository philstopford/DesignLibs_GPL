using System;
using Burkardt.TestValues;

namespace TestValuesTest
{
    public class OwenTest
    {
        public static void owen_values_test()
            //****************************************************************************80
            //
            //  Purpose: 
            //
            //    OWEN_VALUES_TEST tests OWEN_VALUES.
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
            double h = 0;
            int n_data;
            double t = 0;
            Console.WriteLine("");
            Console.WriteLine("OWEN_VALUES_TEST");
            Console.WriteLine("  OWEN_VALUES stores values of");
            Console.WriteLine("  Owen's T function.");
            Console.WriteLine("");
            Console.WriteLine("          H            A            T");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Owen.owen_values(ref n_data, ref h, ref a, ref t);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + h.ToString().PadLeft(12) + "  "
                                  + a.ToString().PadLeft(12) + "  "
                                  + t.ToString("0.################").PadLeft(24) + "");
            }
        }
    }
}