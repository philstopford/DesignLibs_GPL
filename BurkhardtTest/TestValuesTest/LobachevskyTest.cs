using System;
using Burkardt.Values;

namespace TestValuesTest
{
    public static class LobachevskyTest
    {
        public static void lobachevsky_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LOBACHEVSKY_VALUES_TEST tests LOBACHEVSKY_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    13 June 2007
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
            Console.WriteLine("LOBACHEVSKY_VALUES_TEST:");
            Console.WriteLine("  LOBACHEVSKY_VALUES stores values of ");
            Console.WriteLine("  the Lobachevsky function.");
            Console.WriteLine("");
            Console.WriteLine("                X                     FX");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Lobachevsky.lobachevsky_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString().PadLeft(12) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }
    }
}