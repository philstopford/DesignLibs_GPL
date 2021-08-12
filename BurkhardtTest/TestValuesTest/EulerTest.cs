using System;
using Burkardt.TestValues;

namespace TestValuesTest
{
    public static class EulerTest
    {
        public static void euler_number_values_test ( )
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    EULER_NUMBER_VALUES_TEST tests EULER_NUMBER_VALUES.
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
            int c = 0;
            int n = 0;
            int n_data;
            Console.WriteLine("");
            Console.WriteLine("EULER_NUMBER_VALUES_TEST:");
            Console.WriteLine("  EULER_NUMBER_VALUES returns values of");
            Console.WriteLine("  the Euler numbers.");
            Console.WriteLine("");
            Console.WriteLine("     N        EULER_NUMBER(N)");
            Console.WriteLine("");
            n_data = 0;
            for ( ; ; )
            {
                Euler.euler_number_values ( ref n_data, ref n, ref c );
                if ( n_data == 0 )
                {
                    break;
                }
                Console.WriteLine("  "
                                  + n.ToString().PadLeft(6) + "  "
                    + c.ToString().PadLeft(10) + "");
            }
        }
        public static void euler_poly_values_test ( )
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    EULER_POLY_VALUES_TEST tests EULER_POLY_VALUES.
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
            int n = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("EULER_POLY_VALUES_TEST:");
            Console.WriteLine("  EULER_POLY_VALUES returns values of");
            Console.WriteLine("  the Euler numbers.");
            Console.WriteLine("");
            Console.WriteLine("     N     X       EULER_POLY(N)(X)");
            Console.WriteLine("");
            n_data = 0;
            for ( ; ; )
            {
                Euler.euler_poly_values ( ref n_data, ref n, ref x, ref fx );
                if ( n_data == 0 )
                {
                    break;
                }
                Console.WriteLine("  "
                                  + n.ToString().PadLeft(6) + "  "
                    + x.ToString().PadLeft(8) + "  "
                    + fx.ToString().PadLeft(16) + "");
            }
        }
    }
}