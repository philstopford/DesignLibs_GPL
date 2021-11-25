using System;
using System.Globalization;
using Burkardt.Values;

namespace TestValuesTest
{
    public static class SqrtTest
    {
        public static void sqrt_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SQRT_VALUES_TEST tests SQRT_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fx = 0;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("SQRT_VALUES_TEST:");
            Console.WriteLine("  SQRT_VALUES returns some exact values.");
            Console.WriteLine("");
            Console.WriteLine("     X       Fx");
            Console.WriteLine("");
            int n_data = 0;
            for (;;)
            {
                Sqrt.sqrt_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }
    }
}