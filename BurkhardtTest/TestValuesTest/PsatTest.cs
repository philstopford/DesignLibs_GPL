using System;
using TestValues;

namespace TestValuesTest
{
    public class PsatTest
    {
        public static void psat_values_test ( )
            //****************************************************************************80
            //
            //  Purpose: 
            //
            //    PSAT_VALUES_TEST tests PSAT_VALUES.
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
            double psat = 0;
            double tc = 0;
            Console.WriteLine("");
            Console.WriteLine("PSAT_VALUES_TEST:");
            Console.WriteLine("  PSAT_VALUES stores values of");
            Console.WriteLine("  the saturation pressure of water");
            Console.WriteLine("  as a function of temperature.");
            Console.WriteLine("");
            Console.WriteLine("      T            PSAT(T)");
            Console.WriteLine("");
            n_data = 0;
            for ( ; ; )
            {
                Psat.psat_values ( ref n_data, ref tc, ref psat );
                if ( n_data == 0 )
                {
                    break;
                }
                Console.WriteLine( "  "
                    + tc.ToString("0.################").PadLeft(24) +  "  "
                    + psat.ToString("0.################").PadLeft(24) + "");
            }
        }
    }
}