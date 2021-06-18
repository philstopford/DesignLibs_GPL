using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest
{
    partial class Program
    {
        static void log_series_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    LOG_SERIES_CDF_TEST tests LOG_SERIES_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
        {
            double a;
            double cdf;
            int i;
            double pdf;
            int seed = 123456789;
            int x;
            int x2;

            Console.WriteLine("");
            Console.WriteLine("LOG_SERIES_CDF_TEST");
            Console.WriteLine("  LOG_SERIES_CDF evaluates the Log Series CDF;");
            Console.WriteLine("  LOG_SERIES_CDF_INV inverts the Log Series CDF.");
            Console.WriteLine("  LOG_SERIES_PDF evaluates the Log Series PDF;");

            a = 0.25E+00;

            Console.WriteLine("");
            Console.WriteLine("  PDF parameter A =             " + a + "");

            if (!LogSeries.log_series_check(a))
            {
                Console.WriteLine("");
                Console.WriteLine("LOG_SERIES_CDF_TEST - Fatal error!");
                Console.WriteLine("  The parameters are not legal.");
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("       X            PDF           CDF            CDF_INV");
            Console.WriteLine("");

            for (i = 1; i <= 10; i++)
            {
                x = LogSeries.log_series_sample(a, ref seed);
                pdf = LogSeries.log_series_pdf(x, a);
                cdf = LogSeries.log_series_cdf(x, a);
                x2 = LogSeries.log_series_cdf_inv(cdf, a);

                Console.WriteLine("  "
                                  + x.ToString().PadLeft(12) + "  "
                                  + pdf.ToString().PadLeft(12) + "  "
                                  + cdf.ToString().PadLeft(12) + "  "
                                  + x2.ToString().PadLeft(12) + "");
            }

        }

        static void log_series_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    LOG_SERIES_SAMPLE_TEST tests LOG_SERIES_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
        {
            int SAMPLE_NUM = 1000;

            double a;
            int j;
            double mean;
            int seed = 123456789;
            double variance;
            int[] x = new int [SAMPLE_NUM];
            int xmax;
            int xmin;

            Console.WriteLine("");
            Console.WriteLine("LOG_SERIES_SAMPLE_TEST");
            Console.WriteLine("  LOG_SERIES_MEAN computes the Log Series mean;");
            Console.WriteLine("  LOG_SERIES_SAMPLE samples the Log Series distribution;");
            Console.WriteLine("  LOG_SERIES_VARIANCE computes the Log Series variance.");

            a = 0.25E+00;

            Console.WriteLine("");
            Console.WriteLine("  PDF parameter A =             " + a + "");

            if (!LogSeries.log_series_check(a))
            {
                Console.WriteLine("");
                Console.WriteLine("LOG_SERIES_SAMPLE_TEST - Fatal error!");
                Console.WriteLine("  The parameters are not legal.");
                return;
            }

            mean = LogSeries.log_series_mean(a);
            variance = LogSeries.log_series_variance(a);

            Console.WriteLine("  PDF mean =                    " + mean + "");
            Console.WriteLine("  PDF variance =                " + variance + "");

            for (j = 0; j < SAMPLE_NUM; j++)
            {
                x[j] = LogSeries.log_series_sample(a, ref seed);
            }

            mean = typeMethods.i4vec_mean(SAMPLE_NUM, x);
            variance = typeMethods.i4vec_variance(SAMPLE_NUM, x);
            xmax = typeMethods.i4vec_max(SAMPLE_NUM, x);
            xmin = typeMethods.i4vec_min(SAMPLE_NUM, x);

            Console.WriteLine("");
            Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
            Console.WriteLine("  Sample mean =     " + mean + "");
            Console.WriteLine("  Sample variance = " + variance + "");
            Console.WriteLine("  Sample maximum =  " + xmax + "");
            Console.WriteLine("  Sample minimum =  " + xmin + "");

        }

    }
}