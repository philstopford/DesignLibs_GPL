using System;
using TestValues;

namespace TestValuesTest
{
    public class ITest
    {

        public static void i0ml0_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I0ML0_VALUES_TEST tests I0ML0_VALUES.
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
            double fx = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("I0ML0_VALUES_TEST:");
            Console.WriteLine("  I0ML0_VALUES stores values of ");
            Console.WriteLine("  the I0-L0 function.");
            Console.WriteLine("");
            Console.WriteLine("                X                     FX");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                I.i0ml0_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString("0.################").PadLeft(24) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void i1ml1_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I1ML1_VALUES_TEST tests I1ML1_VALUES.
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
            double fx = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("I1ML1_VALUES_TEST:");
            Console.WriteLine("  I1ML1_VALUES stores values of ");
            Console.WriteLine("  the I1-L1 function.");
            Console.WriteLine("");
            Console.WriteLine("                X                     FX");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                I.i1ml1_values(ref n_data, ref x, ref fx);
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