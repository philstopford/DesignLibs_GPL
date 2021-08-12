using System;
using Burkardt.TestValues;

namespace TestValuesTest
{
    public static class BesselTest
    {
        public static void bessel_i0_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_I0_VALUES_TEST tests BESSEL_I0_VALUES.
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
            double fx = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("BESSEL_I0_VALUES_TEST:");
            Console.WriteLine("  BESSEL_I0_VALUES stores values of ");
            Console.WriteLine("  the Bessel I0 function.");
            Console.WriteLine("");
            Console.WriteLine("      X         I0(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Bessel.bessel_i0_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                     + x.ToString().PadLeft(12) + "  "
                     + fx.ToString().PadLeft(12) + "");
            }
        }

        public static void bessel_i0_int_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_I0_INT_VALUES_TEST tests BESSEL_I0_INT_VALUES.
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
            Console.WriteLine("BESSEL_I0_INT_VALUES_TEST:");
            Console.WriteLine("  BESSEL_I0_INT_VALUES stores values of ");
            Console.WriteLine("  the integral of the Bessel I0 function.");
            Console.WriteLine("");
            Console.WriteLine("                X                     FX");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Bessel.bessel_i0_int_values(ref n_data, ref x, ref bip);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString("0.################").PadLeft(24) + "  "
                                  + bip.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void bessel_i0_spherical_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_I0_SPHERICAL_VALUES_TEST tests BESSEL_I0_SPHERICAL_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    06 January 2007
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
            Console.WriteLine("BESSEL_I0_SPHERICAL_VALUES_TEST:");
            Console.WriteLine("  BESSEL_I0_SPHERICAL_VALUES stores values of");
            Console.WriteLine("  the spherical Bessel i0 function.");
            Console.WriteLine("");
            Console.WriteLine("      X            i0(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Bessel.bessel_i0_spherical_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                     + x.ToString().PadLeft(12) + "  "
                     + fx.ToString().PadLeft(12) + "");
            }
        }

        public static void bessel_i1_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_I1_VALUES_TEST tests BESSEL_I1_VALUES.
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
            double fx = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("BESSEL_I1_VALUES_TEST:");
            Console.WriteLine("  BESSEL_I1_VALUES stores values of ");
            Console.WriteLine("  the Bessel I1 function.");
            Console.WriteLine("");
            Console.WriteLine("      X         I1(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Bessel.bessel_i1_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                     + x.ToString().PadLeft(12) + "  "
                     + fx.ToString().PadLeft(12) + "");
            }
        }

        public static void bessel_i1_spherical_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_I1_SPHERICAL_VALUES_TEST tests BESSEL_I1_SPHERICAL_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    06 January 2007
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
            Console.WriteLine("BESSEL_I1_SPHERICAL_VALUES_TEST:");
            Console.WriteLine("  BESSEL_I1_SPHERICAL_VALUES stores values of");
            Console.WriteLine("  the spherical Bessel i1 function.");
            Console.WriteLine("");
            Console.WriteLine("      X            i1(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Bessel.bessel_i1_spherical_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                     + x.ToString().PadLeft(12) + "  "
                     + fx.ToString().PadLeft(12) + "");
            }
        }

        public static void bessel_in_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_IN_VALUES_TEST tests BESSEL_IN_VALUES.
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
            double fx = 0;
            int n = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("BESSEL_IN_VALUES_TEST:");
            Console.WriteLine("  BESSEL_IN_VALUES stores values of ");
            Console.WriteLine("  the Bessel In function.");
            Console.WriteLine("");
            Console.WriteLine("      N     X         IN(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Bessel.bessel_in_values(ref n_data, ref n, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                     + n.ToString().PadLeft(6) + "  "
                     + x.ToString().PadLeft(12) + "  "
                     + fx.ToString().PadLeft(12) + "");
            }
        }

        public static void bessel_ix_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_IX_VALUES_TEST tests BESSEL_IX_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    18 January 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fx = 0;
            double nu = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("BESSEL_IX_VALUES_TEST:");
            Console.WriteLine("  BESSEL_IX_VALUES stores values of ");
            Console.WriteLine("  the Bessel In function for NONINTEGER order.");
            Console.WriteLine("");
            Console.WriteLine("      NU    X         IN(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Bessel.bessel_ix_values(ref n_data, ref nu, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                     + nu.ToString().PadLeft(12) + "  "
                     + x.ToString().PadLeft(12) + "  "
                     + fx.ToString().PadLeft(12) + "");
            }
        }

        public static void bessel_j_spherical_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_J_SPHERICAL_VALUES_TEST tests BESSEL_J_SPHERICAL_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    16 January 2016
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
            Console.WriteLine("BESSEL_J_SPHERICAL_VALUES_TEST:");
            Console.WriteLine("  BESSEL_J_SPHERICAL_VALUES stores values of");
            Console.WriteLine("  the spherical Bessel jn function.");
            Console.WriteLine("");
            Console.WriteLine("      N      X            jn(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Bessel.bessel_j_spherical_values(ref n_data, ref n, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                     + n.ToString().PadLeft(4) + "  "
                     + x.ToString().PadLeft(12) + "  "
                     + fx.ToString().PadLeft(12) + "");
            }
        }

        public static void bessel_j0_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_J0_VALUES_TEST tests BESSEL_J0_VALUES.
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
            double fx = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("BESSEL_J0_VALUES_TEST:");
            Console.WriteLine("  BESSEL_J0_VALUES stores values of ");
            Console.WriteLine("  the Bessel J0 function.");
            Console.WriteLine("");
            Console.WriteLine("      X         J0(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Bessel.bessel_j0_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                     + x.ToString().PadLeft(12) + "  "
                     + fx.ToString().PadLeft(12) + "");
            }
        }

        public static void bessel_j0_int_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_J0_INT_VALUES_TEST tests BESSEL_J0_INT_VALUES.
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
            Console.WriteLine("BESSEL_J0_INT_VALUES_TEST:");
            Console.WriteLine("  BESSEL_J0_INT_VALUES stores values of ");
            Console.WriteLine("  the integral of the Bessel J0 function.");
            Console.WriteLine("");
            Console.WriteLine("                X                     FX");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Bessel.bessel_j0_int_values(ref n_data, ref x, ref bip);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString("0.################").PadLeft(24) + "  "
                                  + bip.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void bessel_j0_spherical_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_J0_SPHERICAL_VALUES_TEST tests BESSEL_J0_SPHERICAL_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    06 January 2007
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
            Console.WriteLine("BESSEL_J0_SPHERICAL_VALUES_TEST:");
            Console.WriteLine("  BESSEL_J0_SPHERICAL_VALUES stores values of");
            Console.WriteLine("  the spherical Bessel j0 function.");
            Console.WriteLine("");
            Console.WriteLine("      X            j0(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Bessel.bessel_j0_spherical_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                     + x.ToString().PadLeft(12) + "  "
                     + fx.ToString().PadLeft(12) + "");
            }
        }

        public static void bessel_j0_zero_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_J0_ZERO_VALUES_TEST tests BESSEL_J0_ZERO_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 January 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fx = 0;
            int k = 0;
            int n_data;
            Console.WriteLine("");
            Console.WriteLine("BESSEL_J0_ZERO_VALUES_TEST:");
            Console.WriteLine("  BESSEL_J0_ZERO_VALUES stores values of zeros of");
            Console.WriteLine("  the Bessel J0 function.");
            Console.WriteLine("");
            Console.WriteLine("       K         X(K)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Bessel.bessel_j0_zero_values(ref n_data, ref k, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  " + k.ToString().PadLeft(6)
                                  + "  " + fx.ToString().PadLeft(24) + "");
            }
        }

        public static void bessel_j1_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_J1_VALUES_TEST tests BESSEL_J1_VALUES.
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
            double fx = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("BESSEL_J1_VALUES_TEST:");
            Console.WriteLine("  BESSEL_J1_VALUES stores values of ");
            Console.WriteLine("  the Bessel J1 function.");
            Console.WriteLine("");
            Console.WriteLine("      X         J1(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Bessel.bessel_j1_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                     + x.ToString().PadLeft(12) + "  "
                     + fx.ToString().PadLeft(12) + "");
            }
        }

        public static void bessel_j1_spherical_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_J1_SPHERICAL_VALUES_TEST tests BESSEL_J1_SPHERICAL_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    06 January 2007
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
            Console.WriteLine("BESSEL_J1_SPHERICAL_VALUES_TEST:");
            Console.WriteLine("  BESSEL_J1_SPHERICAL_VALUES stores values of");
            Console.WriteLine("  the spherical Bessel j1 function.");
            Console.WriteLine("");
            Console.WriteLine("      X            j1(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Bessel.bessel_j1_spherical_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                     + x.ToString().PadLeft(12) + "  "
                     + fx.ToString().PadLeft(12) + "");
            }
        }

        public static void bessel_jn_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_JN_VALUES_TEST tests BESSEL_JN_VALUES.
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
            double fx = 0;
            int n = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("BESSEL_JN_VALUES_TEST:");
            Console.WriteLine("  BESSEL_JN_VALUES stores values of ");
            Console.WriteLine("  the Bessel Jn function.");
            Console.WriteLine("");
            Console.WriteLine("      N     X         JN(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Bessel.bessel_jn_values(ref n_data, ref n, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                     + n.ToString().PadLeft(6) + "  "
                     + x.ToString().PadLeft(12) + "  "
                     + fx.ToString().PadLeft(12) + "");
            }
        }

        public static void bessel_jx_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_JX_VALUES_TEST tests BESSEL_JX_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    01 April 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fx = 0;
            double nu = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("BESSEL_JX_VALUES_TEST:");
            Console.WriteLine("  BESSEL_JX_VALUES stores values of ");
            Console.WriteLine("  the Bessel Jn function for NONINTEGER order.");
            Console.WriteLine("");
            Console.WriteLine("      NU      X         JN(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Bessel.bessel_jx_values(ref n_data, ref nu, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                     + nu.ToString().PadLeft(12) + "  "
                     + x.ToString().PadLeft(12) + "  "
                     + fx.ToString().PadLeft(12) + "");
            }
        }

        public static void bessel_k0_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_K0_VALUES_TEST tests BESSEL_K0_VALUES.
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
            double fx = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("BESSEL_K0_VALUES_TEST:");
            Console.WriteLine("  BESSEL_K0_VALUES stores values of ");
            Console.WriteLine("  the Bessel K0 function.");
            Console.WriteLine("");
            Console.WriteLine("      X         K0(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Bessel.bessel_k0_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                     + x.ToString().PadLeft(12) + "  "
                     + fx.ToString().PadLeft(12) + "");
            }
        }

        public static void bessel_k0_int_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_K0_INT_VALUES_TEST tests BESSEL_K0_INT_VALUES.
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
            Console.WriteLine("BESSEL_K0_INT_VALUES_TEST:");
            Console.WriteLine("  BESSEL_K0_INT_VALUES stores values of ");
            Console.WriteLine("  the integral of the Bessel K0 function.");
            Console.WriteLine("");
            Console.WriteLine("                X                     FX");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Bessel.bessel_k0_int_values(ref n_data, ref x, ref bip);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString("0.################").PadLeft(24) + "  "
                                  + bip.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void bessel_k1_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_K1_VALUES_TEST tests BESSEL_K1_VALUES.
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
            double fx = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("BESSEL_K1_VALUES_TEST:");
            Console.WriteLine("  BESSEL_K1_VALUES stores values of ");
            Console.WriteLine("  the Bessel K1 function.");
            Console.WriteLine("");
            Console.WriteLine("      X         K1(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Bessel.bessel_k1_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                     + x.ToString().PadLeft(12) + "  "
                     + fx.ToString().PadLeft(12) + "");
            }
        }

        public static void bessel_kn_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_KN_VALUES_TEST tests BESSEL_KN_VALUES.
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
            double fx = 0;
            int n = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("BESSEL_KN_VALUES_TEST:");
            Console.WriteLine("  BESSEL_KN_VALUES stores values of ");
            Console.WriteLine("  the Bessel Kn function.");
            Console.WriteLine("");
            Console.WriteLine("      N      X         KN(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Bessel.bessel_kn_values(ref n_data, ref n, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                     + n.ToString().PadLeft(6) + "  "
                     + x.ToString().PadLeft(12) + "  "
                     + fx.ToString().PadLeft(12) + "");
            }
        }

        public static void bessel_kx_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_KX_VALUES_TEST tests BESSEL_KX_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    18 January 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fx = 0;
            double nu = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("BESSEL_KX_VALUES_TEST:");
            Console.WriteLine("  BESSEL_KX_VALUES stores values of ");
            Console.WriteLine("  the Bessel Kn function for NONINTEGER order.");
            Console.WriteLine("");
            Console.WriteLine("      NU     X         KN(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Bessel.bessel_kx_values(ref n_data, ref nu, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                     + nu.ToString().PadLeft(12) + "  "
                     + x.ToString().PadLeft(12) + "  "
                     + fx.ToString().PadLeft(12) + "");
            }
        }

        public static void bessel_y0_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_Y0_VALUES_TEST tests BESSEL_Y0_VALUES.
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
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("BESSEL_Y0_VALUES_TEST:");
            Console.WriteLine("  BESSEL_Y0_VALUES stores values of ");
            Console.WriteLine("  the Bessel Y0 function.");
            Console.WriteLine("");
            Console.WriteLine("      X         Y0(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Bessel.bessel_y0_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                     + x.ToString().PadLeft(12) + "  "
                     + fx.ToString().PadLeft(12) + "");
            }
        }

        public static void bessel_y0_int_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_Y0_INT_VALUES_TEST tests BESSEL_Y0_INT_VALUES.
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
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("BESSEL_Y0_INT_VALUES_TEST:");
            Console.WriteLine("  BESSEL_Y0_INT_VALUES stores values of ");
            Console.WriteLine("  the integral of the Bessel Y0 function.");
            Console.WriteLine("");
            Console.WriteLine("                X                     FX");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Bessel.bessel_y0_int_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString("0.################").PadLeft(24) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void bessel_y0_spherical_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_Y0_SPHERICAL_VALUES_TEST tests BESSEL_Y0_SPHERICAL_VALUES.
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
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("BESSEL_Y0_SPHERICAL_VALUES_TEST:");
            Console.WriteLine("  BESSEL_Y0_SPHERICAL_VALUES stores values of");
            Console.WriteLine("  the spherical Bessel y0 function.");
            Console.WriteLine("");
            Console.WriteLine("                X                      y0(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Bessel.bessel_y0_spherical_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString("0.################").PadLeft(24) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void bessel_y1_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_Y1_VALUES_TEST tests BESSEL_Y1_VALUES.
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
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("BESSEL_Y1_VALUES_TEST:");
            Console.WriteLine("  BESSEL_Y1_VALUES stores values of ");
            Console.WriteLine("  the Bessel Y1 function.");
            Console.WriteLine("");
            Console.WriteLine("                X                   Y1(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Bessel.bessel_y1_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString("0.################").PadLeft(24) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void bessel_y1_spherical_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_Y1_SPHERICAL_VALUES_TEST tests BESSEL_Y1_SPHERICAL_VALUES.
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
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("BESSEL_Y1_SPHERICAL_VALUES_TEST:");
            Console.WriteLine("  BESSEL_Y1_SPHERICAL_VALUES stores values of");
            Console.WriteLine("  the spherical Bessel y1 function.");
            Console.WriteLine("");
            Console.WriteLine("                X                      y1(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Bessel.bessel_y1_spherical_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString("0.################").PadLeft(24) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void bessel_yn_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_YN_VALUES_TEST tests BESSEL_YN_VALUES.
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
            Console.WriteLine("BESSEL_YN_VALUES_TEST:");
            Console.WriteLine("  BESSEL_YN_VALUES stores values of ");
            Console.WriteLine("  the Bessel Yn function.");
            Console.WriteLine("");
            Console.WriteLine("      N     X         YN(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Bessel.bessel_yn_values(ref n_data, ref n, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                     + n.ToString().PadLeft(6) + "  "
                     + x.ToString().PadLeft(12) + "  "
                     + fx.ToString().PadLeft(12) + "");
            }
        }

        public static void bessel_yx_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_YX_VALUES_TEST tests BESSEL_YX_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    01 April 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fx = 0;
            double nu = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("BESSEL_YX_VALUES_TEST:");
            Console.WriteLine("  BESSEL_YX_VALUES stores values of ");
            Console.WriteLine("  the Bessel Yn function for NONINTEGER order.");
            Console.WriteLine("");
            Console.WriteLine("      NU    X         YN(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Bessel.bessel_yx_values(ref n_data, ref nu, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                     + nu.ToString().PadLeft(12) + "  "
                     + x.ToString().PadLeft(12) + "  "
                     + fx.ToString().PadLeft(12) + "");
            }
        }

    }
}