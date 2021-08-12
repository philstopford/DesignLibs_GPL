using System;
using Burkardt.TestValues;

namespace TestValuesTest
{
    public class PrandtlTest
    {
        public static void prandtl_values_test ( )
            //****************************************************************************80
            //
            //  Purpose: 
            //
            //    PRANDTL_VALUES_TEST tests PRANDTL_VALUES.
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
            double p = 0;
            double pr = 0;
            double tc = 0;
            Console.WriteLine("");
            Console.WriteLine("PRANDTL_VALUES_TEST:");
            Console.WriteLine("  PRANDTL_VALUES stores values of");
            Console.WriteLine("  the Prandtl number of water");
            Console.WriteLine("  as a function of temperature and pressure.");
            Console.WriteLine("");
            Console.WriteLine("      T            P            Pr(T,P)");
            Console.WriteLine("");
            n_data = 0;
            for ( ; ; )
            {
                Prandtl.prandtl_values ( ref n_data, ref tc, ref p, ref pr );
                if ( n_data == 0 )
                {
                    break;
                }
                Console.WriteLine("  "
                                  + tc.ToString().PadLeft(12) + "  "
                    + p.ToString().PadLeft(12)  + "  "
                    + pr.ToString().PadLeft(12) + "");
            }
        }
    }
}