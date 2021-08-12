using System;
using Burkardt.TestValues;

namespace TestValuesTest
{
    public class ChebyshevTest
    {
        public static void cheby_t_poly_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CHEBY_T_POLY_VALUES_TEST tests CHEBY_T_POLY_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 March 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fx = 0;
            int n = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("CHEBY_T_POLY_VALUES_TEST:");
            Console.WriteLine("  CHEBY_T_POLY_VALUES returns values of");
            Console.WriteLine("  the Chebyshev T polynomials.");
            Console.WriteLine("");
            Console.WriteLine("     N       X      T(N)(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Chebyshev.cheby_t_poly_values(ref n_data, ref n, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + n.ToString().PadLeft(6) + "  "
                                  + x.ToString().PadLeft(8) + "  "
                                  + fx.ToString().PadLeft(12) + "");
            }
        }

        public static void cheby_t01_poly_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CHEBY_T01_POLY_VALUES_TEST tests CHEBY_T01_POLY_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 July 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fx = 0;
            int n = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("CHEBY_T01_POLY_VALUES_TEST:");
            Console.WriteLine("  CHEBY_T01_POLY_VALUES returns values of");
            Console.WriteLine("  the shifted Chebyshev T polynomials.");
            Console.WriteLine("");
            Console.WriteLine("     N       X      T01(N)(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Chebyshev.cheby_t01_poly_values(ref n_data, ref n, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + n.ToString().PadLeft(6) + "  "
                                  + x.ToString().PadLeft(8) + "  "
                                  + fx.ToString().PadLeft(12) + "");
            }
        }

        public static void cheby_u_poly_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CHEBY_U_POLY_VALUES_TEST tests CHEBY_U_POLY_VALUES.
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
            int n = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("CHEBY_U_POLY_VALUES_TEST:");
            Console.WriteLine("  CHEBY_U_POLY_VALUES returns values of");
            Console.WriteLine("  the Chebyshev U polynomials.");
            Console.WriteLine("");
            Console.WriteLine("     N       X      U(N)(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Chebyshev.cheby_u_poly_values(ref n_data, ref n, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + n.ToString().PadLeft(6) + "  "
                                  + x.ToString().PadLeft(8) + "  "
                                  + fx.ToString().PadLeft(12) + "");
            }
        }

        public static void cheby_u01_poly_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CHEBY_U01_POLY_VALUES_TEST tests CHEBY_U01_POLY_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 July 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fx = 0;
            int n = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("CHEBY_U01_POLY_VALUES_TEST:");
            Console.WriteLine("  CHEBY_U01_POLY_VALUES returns values of");
            Console.WriteLine("  the shifted Chebyshev U polynomials.");
            Console.WriteLine("");
            Console.WriteLine("     N       X      U01(N)(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Chebyshev.cheby_u01_poly_values(ref n_data, ref n, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + n.ToString().PadLeft(6) + "  "
                                  + x.ToString().PadLeft(8) + "  "
                                  + fx.ToString().PadLeft(12) + "");
            }
        }

        public static void cheby_v_poly_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CHEBY_V_POLY_VALUES_TEST tests CHEBY_V_POLY_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 January 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fx = 0;
            int n = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("CHEBY_V_POLY_VALUES_TEST:");
            Console.WriteLine("  CHEBY_V_POLY_VALUES returns values of");
            Console.WriteLine("  the Chebyshev V polynomials.");
            Console.WriteLine("");
            Console.WriteLine("     N       X      V(N)(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Chebyshev.cheby_v_poly_values(ref n_data, ref n, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + n.ToString().PadLeft(6) + "  "
                                  + x.ToString().PadLeft(8) + "  "
                                  + fx.ToString().PadLeft(12) + "");
            }
        }

        public static void cheby_v01_poly_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CHEBY_V01_POLY_VALUES_TEST tests CHEBY_V01_POLY_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 July 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fx = 0;
            int n = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("CHEBY_V01_POLY_VALUES_TEST:");
            Console.WriteLine("  CHEBY_V01_POLY_VALUES returns values of");
            Console.WriteLine("  the shifted Chebyshev V polynomials.");
            Console.WriteLine("");
            Console.WriteLine("     N       X      V01(N)(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Chebyshev.cheby_v01_poly_values(ref n_data, ref n, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + n.ToString().PadLeft(6) + "  "
                                  + x.ToString().PadLeft(8) + "  "
                                  + fx.ToString().PadLeft(12) + "");
            }
        }

        public static void cheby_w_poly_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CHEBY_W_POLY_VALUES_TEST tests CHEBY_W_POLY_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 January 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fx = 0;
            int n = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("CHEBY_W_POLY_VALUES_TEST:");
            Console.WriteLine("  CHEBY_W_POLY_VALUES returns values of");
            Console.WriteLine("  the Chebyshev W polynomials.");
            Console.WriteLine("");
            Console.WriteLine("     N       X      W(N)(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Chebyshev.cheby_w_poly_values(ref n_data, ref n, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + n.ToString().PadLeft(6) + "  "
                                  + x.ToString().PadLeft(8) + "  "
                                  + fx.ToString().PadLeft(12) + "");
            }
        }

        public static void cheby_w01_poly_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CHEBY_W01_POLY_VALUES_TEST tests CHEBY_W01_POLY_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 July 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fx = 0;
            int n = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("CHEBY_W01_POLY_VALUES_TEST:");
            Console.WriteLine("  CHEBY_W01_POLY_VALUES returns values of");
            Console.WriteLine("  the shifted Chebyshev W polynomials.");
            Console.WriteLine("");
            Console.WriteLine("     N       X      W01(N)(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Chebyshev.cheby_w01_poly_values(ref n_data, ref n, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + n.ToString().PadLeft(6) + "  "
                                  + x.ToString().PadLeft(8) + "  "
                                  + fx.ToString().PadLeft(12) + "");
            }
        }

    }
}