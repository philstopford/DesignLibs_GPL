using System;
using Burkardt.Values;

namespace TestValuesTest
{
    public static class GegenbauerTest
    {
        public static void gegenbauer_poly_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GEGENBAUER_POLY_VALUES_TEST tests GEGENBAUER_POLY_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 January 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a = 0;
            double fx = 0;
            int n = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("GEGENBAUER_POLY_VALUES_TEST:");
            Console.WriteLine("  GEGENBAUER_POLY_VALUES returns values of");
            Console.WriteLine("  the Gegenbauer polynomials.");
            Console.WriteLine("");
            Console.WriteLine("       N       A       X       G(N,A)(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Gegenbauer.gegenbauer_poly_values(ref n_data, ref n, ref a, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + n.ToString().PadLeft(6) + "  "
                                  + a.ToString().PadLeft(10) + "  "
                                  + x.ToString().PadLeft(10) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }
    }
}