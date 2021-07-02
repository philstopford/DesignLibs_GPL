using System;
using System.Numerics;
using TestValues;

namespace TestValuesTest
{
    public class WrightOmegaTest
    {
        public static void wright_omega_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    WRIGHT_OMEGA_VALUES_TEST tests WRIGHT_OMEGA_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 May 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Complex fz = new Complex();
            int n_data;
            Complex z = new Complex();
            Console.WriteLine("");
            Console.WriteLine("WRIGHT_OMEGA_VALUES_TEST:");
            Console.WriteLine("  WRIGHT_OMEGA_VALUES stores values of ");
            Console.WriteLine("  the Wright Omega function.");
            Console.WriteLine("");
            Console.WriteLine("                Z                     FZ");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                WrightOmega.wright_omega_values(ref n_data, ref z, ref fz);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + z.Real.ToString("0.######").PadLeft(14) + "  "
                                  + z.Imaginary.ToString("0.######").PadLeft(14) + "  "
                                  + fz.Real.ToString("0.################").PadLeft(24) + "  "
                                  + fz.Imaginary.ToString("0.################").PadLeft(24) + "");
            }
        }

    }
}