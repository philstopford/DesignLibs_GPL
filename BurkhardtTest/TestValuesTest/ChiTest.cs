using System;
using Burkardt.TestValues;

namespace TestValuesTest
{
    public class ChiTest
    {
        public static void chi_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CHI_VALUES_TEST tests CHI_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    11 June 2007
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
            Console.WriteLine("CHI_VALUES_TEST:");
            Console.WriteLine("  CHI_VALUES stores values of");
            Console.WriteLine("  the Hyperbolic Cosine Integral function CHI(X).");
            Console.WriteLine("");
            Console.WriteLine("      X            CHI(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Chi.chi_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString().PadLeft(12) + "  "
                                  + fx.ToString().PadLeft(12) + "");
            }
        }

        public static void chi_square_cdf_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CHI_SQUARE_CDF_VALUES_TEST tests CHI_SQUARE_CDF_VALUES.
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
            int a = 0;
            double fx = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("CHI_SQUARE_CDF_VALUES_TEST:");
            Console.WriteLine("  CHI_SQUARE_CDF_VALUES returns values of ");
            Console.WriteLine("  the Chi-Squared Cumulative Density Function.");
            Console.WriteLine("");
            Console.WriteLine("     N       X    CDF(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Chi.chi_square_cdf_values(ref n_data, ref a, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + a.ToString().PadLeft(6) + "  "
                                  + x.ToString().PadLeft(8) + "  "
                                  + fx.ToString().PadLeft(12) + "");
            }
        }

        public static void chi_square_pdf_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CHI_SQUARE_PDF_VALUES_TEST tests CHI_SQUARE_PDF_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    01 August 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double df = 0;
            double fx = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("CHI_SQUARE_PDF_VALUES_TEST:");
            Console.WriteLine("  CHI_SQUARE_PDF_VALUES returns values of ");
            Console.WriteLine("  the Chi-Squared Probability Density Function.");
            Console.WriteLine("");
            Console.WriteLine("     DF         X    PDF(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Chi.chi_square_pdf_values(ref n_data, ref df, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + df.ToString().PadLeft(8) + "  "
                                  + x.ToString().PadLeft(8) + "  "
                                  + fx.ToString().PadLeft(12) + "");
            }
        }

        public static void chi_square_noncentral_cdf_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CHI_SQUARE_NONCENTRAL_CDF_VALUES_TEST tests CHI_SQUARE_NONCENTRAL_CDF_VALUES.
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
            int df = 0;
            double fx = 0;
            double lambda = 0;
            ;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("CHI_SQUARE_NONCENTRAL_CDF_VALUES_TEST:");
            Console.WriteLine("  CHI_SQUARE_NONCENTRAL_CDF_VALUES returns values of");
            Console.WriteLine("  the noncentral Chi-Squared Cumulative Density Function.");
            Console.WriteLine("");
            Console.WriteLine("      X      LAMBDA     DF     CDF");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Chi.chi_square_noncentral_cdf_values(ref n_data, ref df, ref lambda, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString().PadLeft(10) + "  "
                                  + lambda.ToString().PadLeft(8) + "  "
                                  + df.ToString().PadLeft(4) + "  "
                                  + fx.ToString().PadLeft(12) + "");
            }
        }

        public static void inverse_chi_square_pdf_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    INVERSE_CHI_SQUARE_PDF_VALUES_TEST tests INVERSE_CHI_SQUARE_PDF_VALUES.
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
            double df = 0;
            double fx = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("INVERSE_CHI_SQUARE_PDF_VALUES_TEST:");
            Console.WriteLine("  INVERSE_CHI_SQUARE_PDF_VALUES returns values of ");
            Console.WriteLine("  the inverse Chi-Square Probability Density Function.");
            Console.WriteLine("");
            Console.WriteLine("     DF        X    PDF");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Chi.inverse_chi_square_pdf_values(ref n_data, ref df, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + df.ToString().PadLeft(8) + "  "
                                  + x.ToString().PadLeft(8) + "  "
                                  + fx.ToString().PadLeft(12) + "");
            }
        }

        public static void scaled_inverse_chi_square_pdf_values_test ( )
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SCALED_INVERSE_CHI_SQUARE_PDF_VALUES_TEST tests SCALED_INVERSE_CHI_SQUARE_PDF_VALUES.
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
            double df = 0;
            double fx = 0;
            int n_data;
            double x = 0;
            double xi = 0;
            Console.WriteLine("");
            Console.WriteLine("SCALED_INVERSE_CHI_SQUARE_PDF_VALUES_TEST:");
            Console.WriteLine("  SCALED_INVERSE_CHI_SQUARE_PDF_VALUES returns values of ");
            Console.WriteLine("  the scaled inverse Chi-Square Probability Density Function.");
            Console.WriteLine("");
            Console.WriteLine("     DF        XI        X    PDF");
            Console.WriteLine("");
            n_data = 0;
            for ( ; ; )
            {
                Chi.scaled_inverse_chi_square_pdf_values ( ref n_data, ref df, ref xi, ref x, ref fx );
                if ( n_data == 0 )
                {
                    break;
                }
                Console.WriteLine("  "
                                  + df.ToString().PadLeft(8) + "  "
                    + xi.ToString().PadLeft(8) + "  "
                    + x.ToString().PadLeft(8)  + "  "
                    + fx.ToString().PadLeft(12) + "");
            }
        }
    }
}