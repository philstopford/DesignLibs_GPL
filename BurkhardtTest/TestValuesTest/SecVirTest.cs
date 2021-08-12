using System;
using Burkardt.TestValues;

namespace TestValuesTest
{
    public class SecVirTest
    {
        public static void secvir_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SECVIR_VALUES_TEST tests SECVIR_VALUES.
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
            int n_data;
            double tc = 0;
            double vir = 0;
            Console.WriteLine("");
            Console.WriteLine("SECVIR_VALUES_TEST:");
            Console.WriteLine("  SECVIR_VALUES stores values of");
            Console.WriteLine("  the second virial coefficient of water");
            Console.WriteLine("  as a function of temperature.");
            Console.WriteLine("");
            Console.WriteLine("      T            VIR(T)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                SecVir.secvir_values(ref n_data, ref tc, ref vir);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + tc.ToString().PadLeft(12) + "  "
                                  + vir.ToString().PadLeft(12) + "");
            }
        }
    }
}