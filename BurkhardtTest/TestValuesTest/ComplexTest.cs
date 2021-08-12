using System;
using System.Numerics;
using Burkardt.TestValues;

namespace TestValuesTest
{
    public static class ComplexTest
    {
        public static void c8_log_values_test ( )
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8_LOG_VALUES_TEST tests C8_LOG_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    04 January 2019
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
            Console.WriteLine("C8_LOG_VALUES_TEST:");
            Console.WriteLine("  C8_LOG_VALUES stores values of ");
            Console.WriteLine("  the logarithm of a complex value.");
            Console.WriteLine("");
            Console.WriteLine("                Z                     FZ");
            Console.WriteLine("");
            n_data = 0;
            for ( ; ; )
            {
                Cmplex.c8_log_values ( ref n_data, ref z, ref fz );
                if ( n_data == 0 )
                {
                    break;
                }
                Console.WriteLine("  "
                    + z.Real.ToString("0.######").PadLeft(14) + "  "
                    + z.Imaginary.ToString("0.######").PadLeft(14) + "  "
                    + fz.Real.ToString("0.################").PadLeft(24) + "  "
                    + fz.Imaginary.ToString("################").PadLeft(24) + "");
            }
        }
    }
}