using System;
using TestValues;

namespace TestValuesTest
{
    public class BeiBellTest
    {
        public static void bei0_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BEI0_VALUES_TEST tests BEI0_VALUES.
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
            Console.WriteLine("BEI0_VALUES_TEST:");
            Console.WriteLine("  BEI0_VALUES stores values of ");
            Console.WriteLine("  the Kelvin function BEI of order 0.");
            Console.WriteLine("");
            Console.WriteLine("                X                     FX");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                bei.bei0_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString("0.################").PadLeft(24) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void bei1_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BEI1_VALUES_TEST tests BEI1_VALUES.
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
            Console.WriteLine("BEI1_VALUES_TEST:");
            Console.WriteLine("  BEI1_VALUES stores values of ");
            Console.WriteLine("  the Kelvin function BEI of order 1.");
            Console.WriteLine("");
            Console.WriteLine("                X                     FX");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                bei.bei1_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString("0.################").PadLeft(24) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void bell_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BELL_VALUES_TEST tests BELL_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int c = 0;
            int n = 0;
            int n_data;

            Console.WriteLine("");
            Console.WriteLine("BELL_VALUES_TEST:");
            Console.WriteLine("  BELL_VALUES returns values of ");
            Console.WriteLine("  the Bell numbers.");
            Console.WriteLine("");
            Console.WriteLine("     N        BELL(N)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Bell.bell_values(ref n_data, ref n, ref c);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + n.ToString().PadLeft(6) + "  "
                                  + c.ToString().PadLeft(10) + "");
            }
        }
    }
}