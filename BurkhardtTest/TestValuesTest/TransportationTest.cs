using System;
using Burkardt.Values;

namespace TestValuesTest
{
    public static class TransportationTest
    {
        public static void tran02_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRAN02_VALUES_TEST tests TRAN02_VALUES.
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
            Console.WriteLine("TRAN02_VALUES_TEST:");
            Console.WriteLine("  TRAN02_VALUES stores values of ");
            Console.WriteLine("  the Transport function of order 2.");
            Console.WriteLine("");
            Console.WriteLine("                X                     FX");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Transportation.tran02_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString("0.################").PadLeft(24) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void tran03_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRAN03_VALUES_TEST tests TRAN03_VALUES.
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
            Console.WriteLine("TRAN03_VALUES_TEST:");
            Console.WriteLine("  TRAN03_VALUES stores values of ");
            Console.WriteLine("  the Transport function of order 3.");
            Console.WriteLine("");
            Console.WriteLine("                X                     FX");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Transportation.tran03_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString("0.################").PadLeft(24) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void tran04_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRAN04_VALUES_TEST tests TRAN04_VALUES.
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
            Console.WriteLine("TRAN04_VALUES_TEST:");
            Console.WriteLine("  TRAN04_VALUES stores values of ");
            Console.WriteLine("  the Transport function of order 4.");
            Console.WriteLine("");
            Console.WriteLine("                X                     FX");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Transportation.tran04_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString("0.################").PadLeft(24) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void tran05_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRAN05_VALUES_TEST tests TRAN05_VALUES.
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
            Console.WriteLine("TRAN05_VALUES_TEST:");
            Console.WriteLine("  TRAN05_VALUES stores values of ");
            Console.WriteLine("  the Transport function of order 5.");
            Console.WriteLine("");
            Console.WriteLine("                X                     FX");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Transportation.tran05_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString("0.################").PadLeft(24) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void tran06_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRAN06_VALUES_TEST tests TRAN06_VALUES.
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
            Console.WriteLine("TRAN06_VALUES_TEST:");
            Console.WriteLine("  TRAN06_VALUES stores values of ");
            Console.WriteLine("  the Transport function of order 6.");
            Console.WriteLine("");
            Console.WriteLine("                X                     FX");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Transportation.tran06_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString("0.################").PadLeft(24) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void tran07_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRAN07_VALUES_TEST tests TRAN07_VALUES.
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
            Console.WriteLine("TRAN07_VALUES_TEST:");
            Console.WriteLine("  TRAN07_VALUES stores values of ");
            Console.WriteLine("  the Transport function of order 7.");
            Console.WriteLine("");
            Console.WriteLine("                X                     FX");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Transportation.tran07_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString("0.################").PadLeft(24) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void tran08_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRAN08_VALUES_TEST tests TRAN08_VALUES.
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
            Console.WriteLine("TRAN08_VALUES_TEST:");
            Console.WriteLine("  TRAN08_VALUES stores values of ");
            Console.WriteLine("  the Transport function of order 8.");
            Console.WriteLine("");
            Console.WriteLine("                X                     FX");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Transportation.tran08_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString("0.################").PadLeft(24) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void tran09_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRAN09_VALUES_TEST tests TRAN09_VALUES.
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
            Console.WriteLine("TRAN09_VALUES_TEST:");
            Console.WriteLine("  TRAN09_VALUES stores values of ");
            Console.WriteLine("  the Transport function of order 9.");
            Console.WriteLine("");
            Console.WriteLine("                X                     FX");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Transportation.tran09_values(ref n_data, ref x, ref fx);
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