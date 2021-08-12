using System;
using Burkardt.TestValues;

namespace TestValuesTest
{
    public static class KelvinTest
    {
        public static void kei0_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KEI0_VALUES_TEST tests KEI0_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 June 2006
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
            Console.WriteLine("KEI0_VALUES_TEST:");
            Console.WriteLine("  KEI0_VALUES stores values of ");
            Console.WriteLine("  the Kelvin function KEI of order 0.");
            Console.WriteLine("");
            Console.WriteLine("                X                     FX");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Kelvin.kei0_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString("0.################").PadLeft(24) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void kei1_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KEI1_VALUES_TEST tests KEI1_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 June 2006
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
            Console.WriteLine("KEI1_VALUES_TEST:");
            Console.WriteLine("  KEI1_VALUES stores values of ");
            Console.WriteLine("  the Kelvin function KEI of order 1.");
            Console.WriteLine("");
            Console.WriteLine("                X                     FX");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Kelvin.kei1_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString("0.################").PadLeft(24) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void ker0_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KER0_VALUES_TEST tests KER0_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 June 2006
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
            Console.WriteLine("KER0_VALUES_TEST:");
            Console.WriteLine("  KER0_VALUES stores values of ");
            Console.WriteLine("  the Kelvin function KER of order 0.");
            Console.WriteLine("");
            Console.WriteLine("                X                     FX");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Kelvin.ker0_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString("0.################").PadLeft(24) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void ker1_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KER1_VALUES_TEST tests KER1_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 June 2006
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
            Console.WriteLine("KER1_VALUES_TEST:");
            Console.WriteLine("  KER1_VALUES stores values of ");
            Console.WriteLine("  the Kelvin function KER of order 1.");
            Console.WriteLine("");
            Console.WriteLine("                X                     FX");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Kelvin.ker1_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString("0.################").PadLeft(24) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

    }
}