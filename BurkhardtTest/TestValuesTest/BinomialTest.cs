using System;
using TestValues;

namespace TestValuesTest
{
    public static class BinomialTest
    {
        public static void binomial_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BINOMIAL_VALUES_TEST tests BINOMIAL_VALUES.
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
            int a = 0;
            int b = 0;
            int c = 0;
            int n_data;

            Console.WriteLine("");
            Console.WriteLine("BINOMIAL_VALUES_TEST:");
            Console.WriteLine("  BINOMIAL_VALUES returns values of");
            Console.WriteLine("  the binomial numbers.");
            Console.WriteLine("");
            Console.WriteLine("     A     B        C(A,B)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Binomial.binomial_values(ref n_data, ref a, ref b, ref c);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  " + a.ToString().PadLeft(6)
                                       + "  " + b.ToString().PadLeft(6)
                                       + "  " + c.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void binomial_cdf_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BINOMIAL_CDF_VALUES_TEST tests BINOMIAL_CDF_VALUES.
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
            int a = 0;
            double b = 0;
            double fx = 0;
            int n_data;
            int x = 0;
            Console.WriteLine("");
            Console.WriteLine("BINOMIAL_CDF_VALUES_TEST:");
            Console.WriteLine("  BINOMIAL_CDF_VALUES returns values of ");
            Console.WriteLine("  the Binomial Cumulative Density Function.");
            Console.WriteLine("");
            Console.WriteLine("     A      B        X   CDF(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Binomial.binomial_cdf_values(ref n_data, ref a, ref b, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + a.ToString().PadLeft(6) + "  "
                                  + b.ToString().PadLeft(8) + "  "
                                  + x.ToString().PadLeft(4) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void binomial_pdf_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BINOMIAL_PDF_VALUES_TEST tests BINOMIAL_PDF_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    27 July 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int a = 0;
            double b = 0;
            double fx = 0;
            int n_data;
            int x = 0;
            Console.WriteLine("");
            Console.WriteLine("BINOMIAL_PDF_VALUES_TEST:");
            Console.WriteLine("  BINOMIAL_PDF_VALUES returns values of ");
            Console.WriteLine("  the Binomial Probability Density Function.");
            Console.WriteLine("");
            Console.WriteLine("       A          B                    X       PDF(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Binomial.binomial_pdf_values(ref n_data, ref a, ref b, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + a.ToString().PadLeft(6) + "  "
                                  + b.ToString("0.################").PadLeft(24) + "  "
                                  + x.ToString().PadLeft(4) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }

        public static void negative_binomial_cdf_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NEGATIVE_BINOMIAL_CDF_VALUES_TEST tests NEGATIVE_BINOMIAL_CDF_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    13 June 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double cdf = 0;
            int f = 0;
            int n_data;
            double p = 0;
            int s = 0;
            Console.WriteLine("");
            Console.WriteLine("NEGATIVE_BINOMIAL_CDF_VALUES_TEST:");
            Console.WriteLine("  NEGATIVE_BINOMIAL_CDF_VALUES stores values of");
            Console.WriteLine("  the Negative Binomial Cumulative Density Function.");
            Console.WriteLine("");
            Console.WriteLine("     F     S         P         CDF(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Binomial.negative_binomial_cdf_values(ref n_data, ref f, ref s, ref p, ref cdf);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + f.ToString().PadLeft(6) + "  "
                                  + s.ToString().PadLeft(6) + "  "
                                  + p.ToString().PadLeft(12) + "  "
                                  + cdf.ToString("0.################").PadLeft(24) + "");
            }
        }

    }
}