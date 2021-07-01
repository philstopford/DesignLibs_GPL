using System;
using TestValues;

namespace TestValuesTest
{
    public class BetaTest
    {
        public static void beta_cdf_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BETA_CDF_VALUES_TEST tests BETA_CDF_VALUES.
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
            double a = 0;
            double b = 0;
            double fx = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("BETA_CDF_VALUES_TEST:");
            Console.WriteLine("  BETA_CDF_VALUES stores values of");
            Console.WriteLine("  the Beta CDF.");
            Console.WriteLine("");
            Console.WriteLine("      A            B            X            CDF(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Beta.beta_cdf_values(ref n_data, ref a, ref b, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + a.ToString().PadLeft(12) + "  "
                                  + b.ToString().PadLeft(12) + "  "
                                  + x.ToString().PadLeft(12) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void beta_inc_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BETA_INC_VALUES_TEST tests BETA_INC_VALUES.
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
            double a = 0;
            double b = 0;
            double fx = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("BETA_INC_VALUES_TEST:");
            Console.WriteLine("  BETA_INC_VALUES stores values of");
            Console.WriteLine("  the incomplete Beta function.");
            Console.WriteLine("");
            Console.WriteLine("      A            B            X            BETA_INC(A,B)(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Beta.beta_inc_values(ref n_data, ref a, ref b, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + a.ToString().PadLeft(12) + "  "
                                  + b.ToString().PadLeft(12) + "  "
                                  + x.ToString().PadLeft(12) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void beta_log_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BETA_LOG_VALUES_TEST tests BETA_LOG_VALUES.
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
            double fxy = 0;
            int n_data;
            double x = 0;
            double y = 0;
            Console.WriteLine("");
            Console.WriteLine("BETA_LOG_VALUES_TEST:");
            Console.WriteLine("  BETA_LOG_VALUES stores values of");
            Console.WriteLine("  the logarithm of the Beta function.");
            Console.WriteLine("");
            Console.WriteLine("      X              Y         BETA_LOG(X,Y)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Beta.beta_log_values(ref n_data, ref x, ref y, ref fxy);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString().PadLeft(12) + "  "
                                  + y.ToString().PadLeft(12) + "  "
                                  + fxy.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void beta_noncentral_cdf_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BETA_NONCENTRAL_CDF_VALUES_TEST tests BETA_NONCENTRAL_CDF_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    25 January 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a = 0;
            double b = 0;
            double fx = 0;
            double lambda = 0;
            ;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("BETA_NONCENTRAL_CDF_VALUES_TEST:");
            Console.WriteLine("  BETA_NONCENTRAL_CDF_VALUES stores values of");
            Console.WriteLine("  the noncentral Beta CDF.");
            Console.WriteLine("");
            Console.WriteLine("      A            B       LAMBDA             X            CDF(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Beta.beta_noncentral_cdf_values(ref n_data, ref a, ref b, ref lambda, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + a.ToString().PadLeft(12) + "  "
                                  + b.ToString().PadLeft(12) + "  "
                                  + lambda.ToString().PadLeft(12) + "  "
                                  + x.ToString().PadLeft(12) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void beta_pdf_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BETA_PDF_VALUES_TEST tests BETA_PDF_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 July 2015
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
            Console.WriteLine("BETA_PDF_VALUES_TEST:");
            Console.WriteLine("  BETA_PDF_VALUES stores values of");
            Console.WriteLine("  the Beta PDF.");
            Console.WriteLine("");
            Console.WriteLine("      ALPHA        BETA         X            PDF(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Beta.beta_pdf_values(ref n_data, ref alpha, ref beta, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + alpha.ToString().PadLeft(12) + "  "
                                  + beta.ToString().PadLeft(12) + "  "
                                  + x.ToString().PadLeft(12) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void beta_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BETA_VALUES_TEST tests BETA_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    04 March 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fxy = 0;
            int n_data;
            double x = 0;
            double y = 0;
            Console.WriteLine("");
            Console.WriteLine("BETA_VALUES_TEST:");
            Console.WriteLine("  BETA_VALUES stores values of");
            Console.WriteLine("  the Beta function.");
            Console.WriteLine("");
            Console.WriteLine("      X              Y         BETA(X,Y)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Beta.beta_values(ref n_data, ref x, ref y, ref fxy);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  " + x.ToString().PadLeft(12)
                                  + "  " + y.ToString().PadLeft(12)
                                  + "  " + fxy.ToString("0.################").PadLeft(24) + "");
            }
        }
    }
}