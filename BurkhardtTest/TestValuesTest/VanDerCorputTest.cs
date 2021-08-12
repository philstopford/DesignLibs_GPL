using System;
using Burkardt.TestValues;

namespace TestValuesTest
{
    public class VanDerCorputTest
    {
        public static void van_der_corput_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    VAN_DER_CORPUT_VALUES_TEST tests VAN_DER_CORPUT_VALUES.
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
            int base_ = 0;
            int n_data;
            int seed = 0;
            double value = 0;
            Console.WriteLine("");
            Console.WriteLine("VAN_DER_CORPUT_VALUES_TEST:");
            Console.WriteLine("  VAN_DER_CORPUT_VALUES stores values of");
            Console.WriteLine("  the van der Corput sequence in a given base.");
            Console.WriteLine("");
            Console.WriteLine("      BASE      SEED    VDC(BASE,SEED)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                VanDerCorput.van_der_corput_values(ref n_data, ref base_, ref seed, ref value);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + base_.ToString().PadLeft(8) + "  "
                                  + seed.ToString().PadLeft(8) + "  "
                                  + value.ToString().PadLeft(14) + "");
            }
        }
    }
}