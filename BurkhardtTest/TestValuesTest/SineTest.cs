using System;
using Burkardt.TestValues;

namespace TestValuesTest
{
    public static class SineTest
    {
        public static void si_values_test ( )
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SI_VALUES_TEST tests SI_VALUES.
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
            double fx = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("SI_VALUES_TEST:");
            Console.WriteLine("  SI_VALUES stores values of");
            Console.WriteLine("  the sine integral function.");
            Console.WriteLine("");
            Console.WriteLine("      X            SI(X)");
            Console.WriteLine("");
            n_data = 0;
            for ( ; ; )
            {
                Sine.si_values ( ref n_data, ref x, ref fx );
                if ( n_data == 0 )
                {
                    break;
                }
                Console.WriteLine("  "
                                  + x.ToString("0.################").PadLeft(24) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }
        public static void sin_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIN_VALUES_TEST tests SIN_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 June 2007
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
            Console.WriteLine("SIN_VALUES_TEST:");
            Console.WriteLine("   SIN_VALUES stores values of the sine function.");
            Console.WriteLine("");
            Console.WriteLine("                X                     FX");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Sine.sin_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString("0.################").PadLeft(24) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void sin_degree_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIN_DEGREE_VALUES_TEST tests SIN_DEGREE_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    22 January 2015
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
            Console.WriteLine("SIN_DEGREE_VALUES_TEST:");
            Console.WriteLine("   SIN_DEGREE_VALUES stores values of the sine function.");
            Console.WriteLine("");
            Console.WriteLine("                X                     FX");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Sine.sin_degree_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString("0.################").PadLeft(24) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void sin_power_int_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIN_POWER_INT_VALUES_TEST tests SIN_POWER_INT_VALUES.
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
            double a = 0;
            double b = 0;
            double fx = 0;
            int n = 0;
            int n_data;
            Console.WriteLine("");
            Console.WriteLine("SIN_POWER_INT_VALUES_TEST:");
            Console.WriteLine("  SIN_POWER_INT_VALUES returns values of");
            Console.WriteLine("  the integral of the N-th power of the sine function.");
            Console.WriteLine("");
            Console.WriteLine("         A         B       N        FX");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Sine.sin_power_int_values(ref n_data, ref a, ref b, ref n, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + a.ToString().PadLeft(8) + a + "  "
                                  + b.ToString().PadLeft(8) + b + "  "
                                  + n.ToString().PadLeft(6) + n + "  "
                                  + fx.ToString("0.################").PadLeft(24) + fx + "");
            }
        }

        public static void sinh_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SINH_VALUES_TEST tests SINH_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 June 2007
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
            Console.WriteLine("SINH_VALUES_TEST:");
            Console.WriteLine("   SINH_VALUES stores values of the hyperbolic sine function.");
            Console.WriteLine("");
            Console.WriteLine("                X                     FX");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Sine.sinh_values(ref n_data, ref x, ref fx);
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