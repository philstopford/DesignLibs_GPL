using System;
using TestValues;

namespace TestValuesTest
{
    public static class ArcTest
    {
        public static void arccos_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ARCCOS_VALUES_TEST tests ARCCOS_VALUES.
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
            Console.WriteLine("ARCCOS_VALUES_TEST:");
            Console.WriteLine("  ARCCOS_VALUES stores values of the arc cosine function.");
            Console.WriteLine("");
            Console.WriteLine("                X                     FX");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                arc.arccos_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString("0.################").PadLeft(24) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void arccosh_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ARCCOSH_VALUES_TEST tests ARCCOSH_VALUES.
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
            Console.WriteLine("ARCCOSH_VALUES_TEST:");
            Console.WriteLine(
                "  ARCCOSH_VALUES stores values of the hyperbolic arc cosine function.");
            Console.WriteLine("");
            Console.WriteLine("                X                     FX");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                arc.arccosh_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString("0.################").PadLeft(24) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void arcsin_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ARCSIN_VALUES_TEST tests ARCSIN_VALUES.
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
            Console.WriteLine("ARCSIN_VALUES_TEST:");
            Console.WriteLine("  ARCSIN_VALUES stores values of the arc sine function.");
            Console.WriteLine("");
            Console.WriteLine("                X                     FX");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                arc.arcsin_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString("0.################").PadLeft(24) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void arcsinh_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ARCSINH_VALUES_TEST tests ARCSINH_VALUES.
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
            Console.WriteLine("ARCSINH_VALUES_TEST:");
            Console.WriteLine(
                "  ARCSINH_VALUES stores values of the hyperbolic arc sine function.");
            Console.WriteLine("");
            Console.WriteLine("                X                     FX");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                arc.arcsinh_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString("0.################").PadLeft(24) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void arctan_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ARCTAN_VALUES_TEST tests ARCTAN_VALUES.
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
            Console.WriteLine("ARCTAN_VALUES_TEST:");
            Console.WriteLine("  ARCTAN_VALUES stores values of the arc tangent function.");
            Console.WriteLine("");
            Console.WriteLine("                X                     FX");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                arc.arctan_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString("0.################").PadLeft(24) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void arctan_int_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ARCTAN_INT_VALUES_TEST tests ARCTAN_INT_VALUES.
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
            double bip = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("ARCTAN_INT_VALUES_TEST:");
            Console.WriteLine("  ARCTAN_INT_VALUES stores values of ");
            Console.WriteLine("  the arctangent integral.");
            Console.WriteLine("");
            Console.WriteLine("                X                     FX");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                arc.arctan_int_values(ref n_data, ref x, ref bip);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString("0.################").PadLeft(24) + "  "
                                  + bip.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void arctanh_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ARCTANH_VALUES_TEST tests ARCTANH_VALUES.
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
            Console.WriteLine("ARCTANH_VALUES_TEST:");
            Console.WriteLine(
                "  ARCTANH_VALUES stores values of the hyperbolic arc tangent function.");
            Console.WriteLine("");
            Console.WriteLine("                X                     FX");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                arc.arctanh_values(ref n_data, ref x, ref fx);
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