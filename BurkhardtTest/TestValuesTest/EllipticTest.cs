using System;
using Burkardt.Values;

namespace TestValuesTest
{
    public static class EllipticTest
    {
        public static void elliptic_ea_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ELLIPTIC_EA_VALUES_TEST tests ELLIPTIC_EA_VALUES.
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
            Console.WriteLine("ELLIPTIC_EA_VALUES_TEST:");
            Console.WriteLine("  ELLIPTIC_EA_VALUES stores values of");
            Console.WriteLine("  the complete elliptic integral of the second");
            Console.WriteLine("  kind, with parameter angle A in degrees.");
            Console.WriteLine("");
            Console.WriteLine("    A            E(A)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Elliptic.elliptic_ea_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  " + x.ToString("0.########").PadLeft(12)
                                       + "  " + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void elliptic_ek_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ELLIPTIC_EK_VALUES_TEST tests ELLIPTIC_EK_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    29 May 2018
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
            Console.WriteLine("ELLIPTIC_EK_VALUES_TEST:");
            Console.WriteLine("  ELLIPTIC_EK_VALUES stores values of");
            Console.WriteLine("  the complete elliptic integral of the second");
            Console.WriteLine("  kind, with parameter K.");
            Console.WriteLine("");
            Console.WriteLine("      K            E(K)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Elliptic.elliptic_ek_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  " + x.ToString("0.########").PadLeft(12)
                                       + "  " + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void elliptic_em_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ELLIPTIC_EM_VALUES_TEST tests ELLIPTIC_EM_VALUES.
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
            Console.WriteLine("ELLIPTIC_EM_VALUES_TEST:");
            Console.WriteLine("  ELLIPTIC_EM_VALUES stores values of");
            Console.WriteLine("  the complete elliptic integral of the second");
            Console.WriteLine("  kind, with parameter modulus M.");
            Console.WriteLine("");
            Console.WriteLine("      M            E(M)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Elliptic.elliptic_em_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  " + x.ToString("0.########").PadLeft(12)
                                       + "  " + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void elliptic_fa_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ELLIPTIC_FA_VALUES_TEST tests ELLIPTIC_FA_VALUES.
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
            Console.WriteLine("ELLIPTIC_FA_VALUES_TEST:");
            Console.WriteLine("  ELLIPTIC_FA_VALUES stores values of");
            Console.WriteLine("  the complete elliptic integral of the first");
            Console.WriteLine("  kind, with parameter angle A in degrees.");
            Console.WriteLine("");
            Console.WriteLine("    A            F(A)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Elliptic.elliptic_fa_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  " + x.ToString("0.########").PadLeft(12)
                                       + "  " + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void elliptic_fk_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ELLIPTIC_FK_VALUES_TEST tests ELLIPTIC_FK_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    29 May 2018
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
            Console.WriteLine("ELLIPTIC_FK_VALUES_TEST:");
            Console.WriteLine("  ELLIPTIC_FK_VALUES stores values of");
            Console.WriteLine("  the complete elliptic integral of the first");
            Console.WriteLine("  kind, with parameter K.");
            Console.WriteLine("");
            Console.WriteLine("    K            F(K)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Elliptic.elliptic_fk_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  " + x.ToString("0.########").PadLeft(12)
                                       + "  " + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void elliptic_fm_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ELLIPTIC_FM_VALUES_TEST tests ELLIPTIC_FM_VALUES.
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
            Console.WriteLine("ELLIPTIC_FM_VALUES_TEST:");
            Console.WriteLine("  ELLIPTIC_FM_VALUES stores values of");
            Console.WriteLine("  the complete elliptic integral of the first");
            Console.WriteLine("  kind, with parameter modulus M.");
            Console.WriteLine("");
            Console.WriteLine("      M            F(M)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Elliptic.elliptic_fm_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  " + x.ToString("0.########").PadLeft(12)
                                  + "  " + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void elliptic_inc_ea_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ELLIPTIC_INC_EA_VALUES_TEST tests ELLIPTIC_INC_EA_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    24 June 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a = 0;
            double ea = 0;
            int n_data;
            double phi = 0;

            Console.WriteLine("");
            Console.WriteLine("ELLIPTIC_INC_EA_VALUES_TEST:");
            Console.WriteLine("  ELLIPTIC_INC_EA_VALUES stores values of");
            Console.WriteLine("  the incomplete elliptic integral of the second");
            Console.WriteLine("  kind, with parameters PHI, A.");
            Console.WriteLine("");
            Console.WriteLine("    PHI        A            E(PHI,A)");
            Console.WriteLine("");
            n_data = 0;
            while (true)
            {
                Elliptic.elliptic_inc_ea_values(ref n_data, ref phi, ref a, ref ea);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  " + phi.ToString("0.########").PadLeft(12)
                                       +"  " + a.ToString("0.########").PadLeft(12)
                                       + "  " + ea.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void elliptic_inc_ek_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ELLIPTIC_INC_EK_VALUES_TEST tests ELLIPTIC_INC_EK_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    24 June 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double ek = 0;
            double k = 0;
            int n_data;
            double phi = 0;

            Console.WriteLine("");
            Console.WriteLine("ELLIPTIC_INC_EK_VALUES_TEST:");
            Console.WriteLine("  ELLIPTIC_INC_EK_VALUES stores values of");
            Console.WriteLine("  the incomplete elliptic integral of the second");
            Console.WriteLine("  kind, with parameters PHI, K.");
            Console.WriteLine("");
            Console.WriteLine("    PHI        K            E(PHI,K)");
            Console.WriteLine("");
            n_data = 0;
            while (true)
            {
                Elliptic.elliptic_inc_ek_values(ref n_data, ref phi, ref k, ref ek);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  " + phi.ToString("0.########").PadLeft(12)
                                       +"  " + k.ToString("0.########").PadLeft(12)
                                       + "  " + (-ek).ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void elliptic_inc_em_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ELLIPTIC_INC_EM_VALUES_TEST tests ELLIPTIC_INC_EM_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    24 June 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double em = 0;
            double m = 0;
            int n_data;
            double phi = 0;

            Console.WriteLine("");
            Console.WriteLine("ELLIPTIC_INC_EM_VALUES_TEST:");
            Console.WriteLine("  ELLIPTIC_INC_EM_VALUES stores values of");
            Console.WriteLine("  the incomplete elliptic integral of the second");
            Console.WriteLine("  kind, with parameters PHI, M.");
            Console.WriteLine("");
            Console.WriteLine("    PHI        M            E(PHI,M)");
            Console.WriteLine("");
            n_data = 0;
            while (true)
            {
                Elliptic.elliptic_inc_em_values(ref n_data, ref phi, ref m, ref em);
                if (n_data == 0)
                {
                    break;
                }
                Console.WriteLine("  " + phi.ToString("0.########").PadLeft(12)
                                       +"  " + m.ToString("0.########").PadLeft(12)
                                       + "  " + em.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void elliptic_inc_fa_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ELLIPTIC_INC_FA_VALUES_TEST tests ELLIPTIC_INC_FA_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    24 June 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a = 0;
            double fa = 0;
            int n_data;
            double phi = 0;

            Console.WriteLine("");
            Console.WriteLine("ELLIPTIC_INC_FA_VALUES_TEST:");
            Console.WriteLine("  ELLIPTIC_INC_FA_VALUES stores values of");
            Console.WriteLine("  the incomplete elliptic integral of the first");
            Console.WriteLine("  kind, with parameters PHI, A.");
            Console.WriteLine("");
            Console.WriteLine("    PHI        A            F(PHI,A)");
            Console.WriteLine("");
            n_data = 0;
            while (true)
            {
                Elliptic.elliptic_inc_fa_values(ref n_data, ref phi, ref a, ref fa);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  " + phi.ToString("0.########").PadLeft(12)
                                       +"  " + a.ToString("0.########").PadLeft(12)
                                       + "  " + fa.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void elliptic_inc_fk_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ELLIPTIC_INC_FK_VALUES_TEST tests ELLIPTIC_INC_FK_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    24 June 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fk = 0;
            double k = 0;
            int n_data;
            double phi = 0;

            Console.WriteLine("");
            Console.WriteLine("ELLIPTIC_INC_FK_VALUES_TEST:");
            Console.WriteLine("  ELLIPTIC_INC_FK_VALUES stores values of");
            Console.WriteLine("  the incomplete elliptic integral of the first");
            Console.WriteLine("  kind, with parameters PHI, K.");
            Console.WriteLine("");
            Console.WriteLine("    PHI        K            F(PHI,K)");
            Console.WriteLine("");
            n_data = 0;
            while (true)
            {
                Elliptic.elliptic_inc_fk_values(ref n_data, ref phi, ref k, ref fk);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  " + phi.ToString("0.########").PadLeft(12)
                                       +"  " + k.ToString("0.########").PadLeft(12)
                                       + "  " + fk.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void elliptic_inc_fm_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ELLIPTIC_INC_FM_VALUES_TEST tests ELLIPTIC_INC_FM_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    24 June 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fm = 0;
            double m = 0;
            int n_data;
            double phi = 0;

            Console.WriteLine("");
            Console.WriteLine("ELLIPTIC_INC_FM_VALUES_TEST:");
            Console.WriteLine("  ELLIPTIC_INC_FM_VALUES stores values of");
            Console.WriteLine("  the incomplete elliptic integral of the first");
            Console.WriteLine("  kind, with parameters PHI, M.");
            Console.WriteLine("");
            Console.WriteLine("    PHI        M            F(PHI,M)");
            Console.WriteLine("");
            n_data = 0;
            while (true)
            {
                Elliptic.elliptic_inc_fm_values(ref n_data, ref phi, ref m, ref fm);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  " + phi.ToString("0.########").PadLeft(12)
                                       +"  " + m.ToString("0.########").PadLeft(12)
                                       + "  " + fm.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void elliptic_inc_pia_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ELLIPTIC_INC_PIA_VALUES_TEST tests ELLIPTIC_INC_PIA_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    24 June 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a = 0;
            double n = 0;
            int n_data;
            double phi = 0;
            double pia = 0;

            Console.WriteLine("");
            Console.WriteLine("ELLIPTIC_INC_PIA_VALUES_TEST:");
            Console.WriteLine("  ELLIPTIC_INC_PIA_VALUES stores values of");
            Console.WriteLine("  the incomplete elliptic integral of the third");
            Console.WriteLine("  kind, with parameters PHI, N and A.");
            Console.WriteLine("");
            Console.WriteLine("    PHI           N             A            Pi(PHI,N,A)");
            Console.WriteLine("");
            n_data = 0;
            while (true)
            {
                Elliptic.elliptic_inc_pia_values(ref n_data, ref phi, ref n, ref a, ref pia);
                if (n_data == 0)
                {
                    break;
                }
                Console.WriteLine("  " + phi.ToString("0.########").PadLeft(12)
                                       +"  " + n.ToString("0.########").PadLeft(12)
                                       +"  " + a.ToString("0.########").PadLeft(12)
                                       + "  " + pia.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void elliptic_inc_pik_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ELLIPTIC_INC_PIK_VALUES_TEST tests ELLIPTIC_INC_PIK_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    24 June 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double k = 0;
            double n = 0;
            int n_data;
            double phi = 0;
            double pik = 0;

            Console.WriteLine("");
            Console.WriteLine("ELLIPTIC_INC_PIK_VALUES_TEST:");
            Console.WriteLine("  ELLIPTIC_INC_PIK_VALUES stores values of");
            Console.WriteLine("  the incomplete elliptic integral of the third");
            Console.WriteLine("  kind, with parameters PHI, N and K.");
            Console.WriteLine("");
            Console.WriteLine("    PHI           N             K            Pi(PHI,N,K)");
            Console.WriteLine("");
            n_data = 0;
            while (true)
            {
                Elliptic.elliptic_inc_pik_values(ref n_data, ref phi, ref n, ref k, ref pik);
                if (n_data == 0)
                {
                    break;
                }
                Console.WriteLine("  " + phi.ToString("0.########").PadLeft(12)
                                       +"  " + n.ToString("0.########").PadLeft(12)
                                       +"  " + k.ToString("0.########").PadLeft(12)
                                       + "  " + pik.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void elliptic_inc_pim_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ELLIPTIC_INC_PIM_VALUES_TEST tests ELLIPTIC_INC_PIM_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    24 June 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double m = 0;
            double n = 0;
            int n_data;
            double phi = 0;
            double pim = 0;

            Console.WriteLine("");
            Console.WriteLine("ELLIPTIC_INC_PIM_VALUES_TEST:");
            Console.WriteLine("  ELLIPTIC_INC_PIM_VALUES stores values of");
            Console.WriteLine("  the incomplete elliptic integral of the third");
            Console.WriteLine("  kind, with parameters PHI, N and M.");
            Console.WriteLine("");
            Console.WriteLine("    PHI           N             M            Pi(PHI,N,M)");
            Console.WriteLine("");
            n_data = 0;
            while (true)
            {
                Elliptic.elliptic_inc_pim_values(ref n_data, ref phi, ref n, ref m, ref pim);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  " + phi.ToString("0.########").PadLeft(12)
                                       +"  " + n.ToString("0.########").PadLeft(12)
                                       +"  " + m.ToString("0.########").PadLeft(12)
                                       + "  " + pim.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void elliptic_pia_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ELLIPTIC_PIA_VALUES_TEST tests ELLIPTIC_PIA_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    31 May 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a = 0;
            double n = 0;
            int n_data;
            double pia = 0;
            Console.WriteLine("");
            Console.WriteLine("ELLIPTIC_PIA_VALUES_TEST:");
            Console.WriteLine("  ELLIPTIC_PIA_VALUES stores values of");
            Console.WriteLine("  the complete elliptic integral of the third");
            Console.WriteLine("  kind, with parameter angle A in degrees.");
            Console.WriteLine("");
            Console.WriteLine("    N             A            Pi(N,A)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Elliptic.elliptic_pia_values(ref n_data, ref n, ref a, ref pia);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  " + n.ToString("0.########").PadLeft(12)
                                       +"  " + a.ToString("0.########").PadLeft(12)
                                       + "  " + pia.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void elliptic_pik_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ELLIPTIC_PIK_VALUES_TEST tests ELLIPTIC_PIK_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    31 May 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double k = 0;
            double n = 0;
            int n_data;
            double pik = 0;
            Console.WriteLine("");
            Console.WriteLine("ELLIPTIC_PIK_VALUES_TEST:");
            Console.WriteLine("  ELLIPTIC_PIK_VALUES stores values of");
            Console.WriteLine("  the complete elliptic integral of the third");
            Console.WriteLine("  kind, with parameter K.");
            Console.WriteLine("");
            Console.WriteLine("    N             K            Pi(N,K)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Elliptic.elliptic_pik_values(ref n_data, ref n, ref k, ref pik);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  " + n.ToString("0.########").PadLeft(12)
                                       +"  " + k.ToString("0.########").PadLeft(12)
                                       + "  " + pik.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void elliptic_pim_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ELLIPTIC_PIM_VALUES_TEST tests ELLIPTIC_PIM_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    31 May 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double m = 0;
            double n = 0;
            int n_data;
            double pim = 0;
            Console.WriteLine("");
            Console.WriteLine("ELLIPTIC_PIM_VALUES_TEST:");
            Console.WriteLine("  ELLIPTIC_PIM_VALUES stores values of");
            Console.WriteLine("  the complete elliptic integral of the third");
            Console.WriteLine("  kind, with parameter M.");
            Console.WriteLine("");
            Console.WriteLine("    N             M            Pi(N,M)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Elliptic.elliptic_pim_values(ref n_data, ref n, ref m, ref pim);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  " + n.ToString("0.########").PadLeft(12)
                                       +"  " + m.ToString("0.########").PadLeft(12)
                                       + "  " + pim.ToString("0.################").PadLeft(24) + "");
            }
        }

    }
}