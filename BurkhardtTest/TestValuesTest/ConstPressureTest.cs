using System;
using Burkardt.Values;

namespace TestValuesTest
{
    public static class ConstPressureTest
    {
        public static void cp_values_test ( )
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CP_VALUES_TEST tests CP_VALUES.
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
            double cp = 0;
            int n_data;
            double p = 0;
            double tc = 0;
            Console.WriteLine("");
            Console.WriteLine("CP_VALUES_TEST:");
            Console.WriteLine("  CP_VALUES stores values of");
            Console.WriteLine("  the specific heat CP");
            Console.WriteLine("  as a function of temperature and pressure.");
            Console.WriteLine("");
            Console.WriteLine("      T            P            CP(T,P)");
            Console.WriteLine("");
            n_data = 0;
            for ( ; ; )
            {
                ConstPressure.cp_values ( ref n_data, ref tc, ref p, ref cp );
                if ( n_data == 0 )
                {
                    break;
                }
                Console.WriteLine("  "
                                  + tc.ToString().PadLeft(12) + "  "
                    + p.ToString().PadLeft(12) + "  "
                    + cp.ToString().PadLeft(12) + "");
            }
        }

    }
}