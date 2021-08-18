using System;
using Burkardt.Values;

namespace TestValuesTest
{
    public class PiTest
    {
        public static void pi_values_test ( )
            //****************************************************************************80
            //
            //  Purpose: 
            //
            //    PI_VALUES_TEST tests PI_VALUES.
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
            int fn = 0;
            int n = 0;
            int n_data;
            Console.WriteLine("");
            Console.WriteLine("PI_VALUES_TEST:");
            Console.WriteLine("  PI_VALUES returns values of");
            Console.WriteLine("  the PI function.");
            Console.WriteLine("");
            Console.WriteLine("     N         PI(N)");
            Console.WriteLine("");
            n_data = 0;
            for ( ; ; )
            {
                Pi.pi_values ( ref n_data, ref n, ref fn );
                if ( n_data == 0 )
                {
                    break;
                }
                Console.WriteLine("  "
                                  + n.ToString().PadLeft(12) + n  + "  "
                    + fn.ToString().PadLeft(10) + fn + "");
            }
        }
    }
}