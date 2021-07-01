using System;
using TestValues;

namespace TestValuesTest
{
    public static class ExpIntegralTest
    {
        public static void e1_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    E1_VALUES_TEST tests E1_VALUES.
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
            double fx = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("E1_VALUES_TEST:");
            Console.WriteLine("  E1_VALUES stores values of");
            Console.WriteLine("  the exponential integral function E1(X).");
            Console.WriteLine("");
            Console.WriteLine("      X          E1(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                ExpIntegral.e1_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString().PadLeft(12) + "  "
                                  + fx.ToString().PadLeft(12) + "");
            }
        }

        public static void ei_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    EI_VALUES_TEST tests EI_VALUES.
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
            double fx = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("EI_VALUES_TEST:");
            Console.WriteLine("  EI_VALUES stores values of");
            Console.WriteLine("  the exponential integral function EI(X).");
            Console.WriteLine("");
            Console.WriteLine("      X          EI(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                ExpIntegral.ei_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString().PadLeft(12) + "  "
                                  + fx.ToString().PadLeft(12) + "");
            }
        }
    }
}