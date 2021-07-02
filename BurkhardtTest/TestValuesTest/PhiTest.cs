using System;
using TestValues;

namespace TestValuesTest
{
    public class PhiTest
    {
        public static void phi_values_test ( )
            //****************************************************************************80
            //
            //  Purpose: 
            //
            //    PHI_VALUES_TEST tests PHI_VALUES.
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
            Console.WriteLine("PHI_VALUES_TEST:");
            Console.WriteLine("  PHI_VALUES returns values of");
            Console.WriteLine("  the PHI function.");
            Console.WriteLine("");
            Console.WriteLine("     N         PHI(N)");
            Console.WriteLine("");
            n_data = 0;
            for ( ; ; )
            {
                Phi.phi_values ( ref n_data, ref n, ref fn );
                if ( n_data == 0 )
                {
                    break;
                }
                Console.WriteLine("  "
                                  + n.ToString().PadLeft(6)  + "  "
                    + fn.ToString().PadLeft(10) + "");
            }
        }
    }
}