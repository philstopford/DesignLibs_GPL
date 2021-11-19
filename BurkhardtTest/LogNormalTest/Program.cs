using System;
using Burkardt.CDFLib;
using Burkardt.PDFLib;
using Burkardt.Types;

namespace LogNormalTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for LOG_NORMAL_TEST.
        //
        //  Discussion:
        //
        //    LOG_NORMAL_TEST tests the LOG_NORMAL library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 March 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("LOG_NORMAL_TEST");
        Console.WriteLine("  Test the LOG_NORMAL library.");

        log_normal_cdf_test();
        log_normal_sample_test();

        Console.WriteLine("");
        Console.WriteLine("LOG_NORMAL_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void log_normal_cdf_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_NORMAL_CDF_TEST tests LOG_NORMAL_CDF.
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
        double cdf;
        int i;
        double mu;
        double pdf;
        int seed = 123456789;
        double sigma;
        double x;
        double x2;

        Console.WriteLine("");
        Console.WriteLine("LOG_NORMAL_CDF_TEST");
        Console.WriteLine("  LOG_NORMAL_CDF evaluates the Log Normal CDF;");
        Console.WriteLine("  LOG_NORMAL_CDF_INV inverts the Log Normal CDF.");
        Console.WriteLine("  LOG_NORMAL_PDF evaluates the Log Normal PDF;");

        mu = 10.0;
        sigma = 2.25;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter MU =      " + mu + "");
        Console.WriteLine("  PDF parameter SIGMA =   " + sigma + "");

        if (!PDF.log_normal_check(mu, sigma))
        {
            Console.WriteLine("");
            Console.WriteLine("LOG_NORMAL_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            x = PDF.log_normal_sample(mu, sigma, ref seed);
            pdf = PDF.log_normal_pdf(x, mu, sigma);
            cdf = CDF.log_normal_cdf(x, mu, sigma);
            x2 = CDF.log_normal_cdf_inv(cdf, mu, sigma);

            Console.WriteLine("  "
                              + x.ToString().PadLeft(12) + "  "
                              + pdf.ToString().PadLeft(12) + "  "
                              + cdf.ToString().PadLeft(12) + "  "
                              + x2.ToString().PadLeft(12) + "");
        }
    }

    private static void log_normal_sample_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_NORMAL_SAMPLE_TEST tests LOG_NORMAL_SAMPLE.
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

        int i;
        double mean;
        double mu;
        int seed = 123456789;
        double sigma;
        double variance;
        double[] x = new double[SAMPLE_NUM];
        double xmax;
        double xmin;

        Console.WriteLine("");
        Console.WriteLine("LOG_NORMAL_SAMPLE_TEST");
        Console.WriteLine("  LOG_NORMAL_MEAN computes the Log Normal mean;");
        Console.WriteLine("  LOG_NORMAL_SAMPLE samples the Log Normal distribution;");
        Console.WriteLine("  LOG_NORMAL_VARIANCE computes the Log Normal variance;");

        mu = 1.0;
        sigma = 2.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter MU =      " + mu + "");
        Console.WriteLine("  PDF parameter SIGMA =   " + sigma + "");

        if (!PDF.log_normal_check(mu, sigma))
        {
            Console.WriteLine("");
            Console.WriteLine("LOG_NORMAL_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        mean = PDF.log_normal_mean(mu, sigma);
        variance = PDF.log_normal_variance(mu, sigma);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = PDF.log_normal_sample(mu, sigma, ref seed);
        }

        mean = typeMethods.r8vec_mean(SAMPLE_NUM, x);
        variance = typeMethods.r8vec_variance(SAMPLE_NUM, x);
        xmax = typeMethods.r8vec_max(SAMPLE_NUM, x);
        xmin = typeMethods.r8vec_min(SAMPLE_NUM, x);

        Console.WriteLine("");
        Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
        Console.WriteLine("  Sample mean =     " + mean + "");
        Console.WriteLine("  Sample variance = " + variance + "");
        Console.WriteLine("  Sample maximum =  " + xmax + "");
        Console.WriteLine("  Sample minimum =  " + xmin + "");
    }
}