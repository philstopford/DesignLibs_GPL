using System;
using Burkardt.Values;

namespace TestValuesTest
{
    public static class CubeRootTest
    {
        public static void cbrt_values_test ( )
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CBRT_VALUES_TEST tests CBRT_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    22 June 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fx = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("CBRT_VALUES_TEST:");
            Console.WriteLine("  CBRT_VALUES stores values of the cube root function.");
            Console.WriteLine("");
            Console.WriteLine("      X            CBRT(X)");
            Console.WriteLine("");
            n_data = 0;
            for ( ; ; )
            {
                CubeRoot.cbrt_values ( ref n_data, ref x, ref fx );
                if ( n_data == 0 )
                {
                    break;
                }
                Console.WriteLine("  "
                                  + x.ToString().PadLeft(12) + "  "
                    + fx.ToString().PadLeft(12) + "");
            }
        }
    }
}