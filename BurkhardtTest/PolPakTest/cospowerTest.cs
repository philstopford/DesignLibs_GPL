using System;
using Burkardt.IntegralNS;

namespace PolPakTest
{
    public static class cospowerTest
    {

        public static void cos_power_int_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    COS_POWER_INT_TEST tests COS_POWER_INT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    31 March 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a = 0;
            double b = 0;
            double fx = 0;
            double fx2 = 0;
            int n = 0;
            int n_data;

            Console.WriteLine("");
            Console.WriteLine("COS_POWER_INT_TEST:");
            Console.WriteLine("  COS_POWER_INT computes the integral of the N-th power");
            Console.WriteLine("  of the cosine function.");
            Console.WriteLine("");
            Console.WriteLine("         A         B       N        Exact    Computed");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                Burkardt.Values.Cosine.cos_power_int_values(ref n_data, ref a, ref b, ref n, ref fx);

                if (n_data == 0)
                {
                    break;
                }

                fx2 = CosPower.cos_power_int(a, b, n);

                Console.WriteLine("  "
                                  + a.ToString().PadLeft(8) + "  "
                                  + b.ToString().PadLeft(8) + "  "
                                  + n.ToString().PadLeft(6) + "  "
                                  + fx.ToString().PadLeft(12) + "  "
                                  + fx2.ToString().PadLeft(12) + "");
            }
        }

    }
}