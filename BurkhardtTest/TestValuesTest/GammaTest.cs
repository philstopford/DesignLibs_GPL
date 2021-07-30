using System;
using Burkardt.AppliedStatistics;
using TestValues;

namespace TestValuesTest
{
    public static class GammaTest
    {
        public static void gamma_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GAMMA_VALUES_TEST tests GAMMA_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 January 2007
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
            Console.WriteLine("GAMMA_VALUES_TEST:");
            Console.WriteLine("  GAMMA_VALUES stores values of the Gamma function.");
            Console.WriteLine("");
            Console.WriteLine("      X            GAMMA(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Gamma.gamma_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  +x.ToString().PadLeft(12) + x + "  "
                    + fx.ToString("0.################").PadLeft(24) + fx + "");
            }
        }

        public static void gamma_01_pdf_values_test()
            //****************************************************************************80
            //
            //  Purpose: 
            //
            //    GAMMA_01_PDF_VALUES_TEST tests GAMMA_01_PDF_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    28 July 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double alpha = 0;
            ;
            double fx = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("GAMMA_01_PDF_VALUES_TEST:");
            Console.WriteLine("  GAMMA_01_PDF_VALUES stores values of");
            Console.WriteLine("  the standard Gamma Probability Density Function.");
            Console.WriteLine("");
            Console.WriteLine("       ALPHA        X                   PDF(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Gamma.gamma_01_pdf_values(ref n_data, ref alpha, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                    + alpha.ToString().PadLeft(12) + alpha + "  "
                    + x.ToString().PadLeft(12) + x + "  "
                    + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void gamma_cdf_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GAMMA_CDF_VALUES_TEST tests GAMMA_CDF_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 January 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fx = 0;
            double mu = 0;
            int n_data;
            double sigma = 0;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("GAMMA_CDF_VALUES_TEST:");
            Console.WriteLine("  GAMMA_CDF_VALUES stores values of");
            Console.WriteLine("  the Gamma CDF.");
            Console.WriteLine("");
            Console.WriteLine("      M    Sigma      X            CDF((X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Gamma.gamma_cdf_values(ref n_data, ref mu, ref sigma, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                    + mu.ToString().PadLeft(12) + mu + "  "
                    + sigma.ToString().PadLeft(12) + sigma + "  "
                    + x.ToString().PadLeft(12) + x + "  "
                    + fx.ToString("0.################").PadLeft(24) + fx + "");
            }
        }

        public static void gamma_inc_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GAMMA_INC_VALUES_TEST tests GAMMA_INC_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 January 2007
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
            Console.WriteLine("GAMMA_INC_VALUES_TEST:");
            Console.WriteLine("   GAMMA_INC_VALUES stores values of");
            Console.WriteLine("   the incomplete Gamma function.");
            Console.WriteLine("");
            Console.WriteLine("      A            X            GAMMA_INC(A)(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Gamma.gamma_inc_values(ref n_data, ref a, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  +a.ToString().PadLeft(12) + a + "  "
                    + x.ToString().PadLeft(12) + x + "  "
                    + fx.ToString("0.################").PadLeft(24) + fx + "");
            }
        }

        public static void gamma_inc_p_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GAMMA_INC_P_VALUES_TEST tests GAMMA_INC_P_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 January 2007
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
            Console.WriteLine("GAMMA_INC_P_VALUES_TEST:");
            Console.WriteLine("   GAMMA_INC_P_VALUES stores values of");
            Console.WriteLine("   the incomplete Gamma P function.");
            Console.WriteLine("");
            Console.WriteLine("      A            X            F(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Gamma.gamma_inc_p_values(ref n_data, ref a, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  +a.ToString().PadLeft(12) + "  "
                    + x.ToString().PadLeft(12) + "  "
                    + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void gamma_inc_q_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GAMMA_INC_Q_VALUES_TEST tests GAMMA_INC_Q_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 January 2007
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
            Console.WriteLine("GAMMA_INC_Q_VALUES_TEST:");
            Console.WriteLine("   GAMMA_INC_Q_VALUES stores values of");
            Console.WriteLine("   the incomplete Gamma Q function.");
            Console.WriteLine("");
            Console.WriteLine("      A            X            F(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Gamma.gamma_inc_q_values(ref n_data, ref a, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  +a.ToString().PadLeft(12) + "  "
                    + x.ToString().PadLeft(12) + "  "
                    + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void gamma_inc_tricomi_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GAMMA_INC_TRICOMI_VALUES_TEST tests GAMMA_INC_TRICOMI_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 January 2007
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
            Console.WriteLine("GAMMA_INC_TRICOMI_VALUES_TEST:");
            Console.WriteLine("   GAMMA_INC_TRICOMI_VALUES stores values of");
            Console.WriteLine("   the incomplete Tricomi Gamma function.");
            Console.WriteLine("");
            Console.WriteLine("      A            X            F(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Gamma.gamma_inc_tricomi_values(ref n_data, ref a, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  +a.ToString().PadLeft(12) + "  "
                    + x.ToString().PadLeft(12) + "  "
                    + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void gamma_log_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GAMMA_LOG_VALUES_TEST tests GAMMA_LOG_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 January 2007
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
            Console.WriteLine("GAMMA_LOG_VALUES_TEST:");
            Console.WriteLine("  GAMMA_LOG_VALUES stores values of");
            Console.WriteLine("  the logarithm of the Gamma function.");
            Console.WriteLine("");
            Console.WriteLine("      X            GAMMA_LOG(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Algorithms.gamma_log_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  +x.ToString().PadLeft(12) + x + "  "
                    + fx.ToString("0.################").PadLeft(24) + fx + "");
            }
        }

        public static void gamma_pdf_values_test()
            //****************************************************************************80
            //
            //  Purpose: 
            //
            //    GAMMA_PDF_VALUES_TEST tests GAMMA_PDF_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    29 July 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double alpha = 0;
            ;
            double beta = 0;
            double fx = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("GAMMA_PDF_VALUES_TEST:");
            Console.WriteLine("  GAMMA_PDF_VALUES stores values of");
            Console.WriteLine("  a Gamma Probability Density Function.");
            Console.WriteLine("");
            Console.WriteLine("       BETA          ALPHA        X                   PDF(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Gamma.gamma_pdf_values(ref n_data, ref beta, ref alpha, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                    + beta.ToString().PadLeft(12) + "  "
                    + alpha.ToString().PadLeft(12) + "  "
                    + x.ToString().PadLeft(12) + "  "
                    + fx.ToString("0.################").PadLeft(24) + "");
            }
        }
        
        public static void inverse_gamma_pdf_values_test ( )
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    INVERSE_GAMMA_PDF_VALUES_TEST tests INVERSE_GAMMA_PDF_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 August 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double alpha = 0;;
            double beta = 0;
            double fx = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("INVERSE_GAMMA_PDF_VALUES_TEST:");
            Console.WriteLine("  INVERSE_GAMMA_PDF_VALUES returns values of ");
            Console.WriteLine("  the inverse gamma Probability Density Function.");
            Console.WriteLine("");
            Console.WriteLine("     ALPHA        BETA      X    PDF");
            Console.WriteLine("");
            n_data = 0;
            for ( ; ; )
            {
                Gamma.inverse_gamma_pdf_values ( ref n_data, ref alpha, ref beta, ref x, ref fx );
                if ( n_data == 0 )
                {
                    break;
                }
                Console.WriteLine("  "
                                  + alpha.ToString().PadLeft(8) + "  "
                    + beta.ToString().PadLeft(8) + "  "
                    + x.ToString().PadLeft(8) + "  "
                    + fx.ToString().PadLeft(12) + "");
            }
        }
    }
}