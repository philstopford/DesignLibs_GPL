using System;
using Burkardt.Values;

namespace TestValuesTest
{
    public static class DielectricTest
    {
        public static void dielectric_values_test ( )
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DIELECTRIC_VALUES_TEST tests DIELECTRIC_VALUES.
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
            double eps = 0;
            int n_data;
            double p = 0;
            double tc = 0;
            Console.WriteLine("");
            Console.WriteLine("DIELECTRIC_VALUES_TEST:");
            Console.WriteLine("  DIELECTRIC_VALUES stores values of");
            Console.WriteLine("  the dielectric function.");
            Console.WriteLine("");
            Console.WriteLine("      T           P            EPS(T,P)");
            Console.WriteLine("");
            n_data = 0;
            for ( ; ; )
            {
                Dielectrics.dielectric_values ( ref n_data, ref tc, ref p, ref eps );
                if ( n_data == 0 )
                {
                    break;
                }
                Console.WriteLine("  "
                                  + tc.ToString().PadLeft(12)  + "  "
                    + p.ToString().PadLeft(12) + "  "
                    + eps.ToString().PadLeft(12) + "");
            }
        }
    }
}