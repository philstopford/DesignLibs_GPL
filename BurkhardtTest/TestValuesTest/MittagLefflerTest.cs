using System;
using Burkardt.TestValues;

namespace TestValuesTest
{
    public class MittagLefflerTest
    {
        public static void mittag_leffler_ea_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MITTAG_LEFFLER_EA_VALUES_TEST tests MITTAG_LEFFLER_EA_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 February 2017
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a = 0;
            double fx = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("MITTAG_LEFFLER_EA_VALUES_TEST:");
            Console.WriteLine("   MITTAG_LEFFLER_EA_VALUES stores values of");
            Console.WriteLine("   the Mittag-Leffler one-parameter function E(A;X).");
            Console.WriteLine("");
            Console.WriteLine("      A            X            E(A;X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                MittagLeffler.mittag_leffler_ea_values(ref n_data, ref a, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + a.ToString().PadLeft(12) + "  "
                                  + x.ToString().PadLeft(12) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void mittag_leffler_eab_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MITTAG_LEFFLER_EAB_VALUES_TEST tests MITTAG_LEFFLER_EAB_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    04 February 2017
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a = 0;
            double b = 0;
            double fx = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("MITTAG_LEFFLER_EAB_VALUES_TEST:");
            Console.WriteLine("   MITTAG_LEFFLER_EAB_VALUES stores values of");
            Console.WriteLine("   the Mittag-Leffler two-parameter function E(A,B;X).");
            Console.WriteLine("");
            Console.WriteLine("      A            B             X            E(A,B;X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                MittagLeffler.mittag_leffler_eab_values(ref n_data, ref a, ref b, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  " + a.ToString().PadLeft(12)
                                       + "  " + b.ToString().PadLeft(12)
                                       + "  " + x.ToString().PadLeft(12)
                                       + "  " + fx.ToString("0.################").PadLeft(24) + "");
            }
        }
    }
}