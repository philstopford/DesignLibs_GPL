using System;
using Burkardt.CDFLib;
using Burkardt.PDFLib;
using Burkardt.Types;

namespace LogNormalTruncatedABTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for LOG_NORMAL_TRUNCATED_AB_TEST.
        //
        //  Discussion:
        //
        //    LOG_NORMAL_TRUNCATED_AB_TEST tests the LOG_NORMAL_TRUNCATED_AB library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 March 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("LOG_NORMAL_TRUNCATED_AB_TEST");
        Console.WriteLine("  Test the LOG_NORMAL_TRUNCATED_AB library.");

        log_normal_truncated_ab_cdf_test();
        log_normal_truncated_ab_sample_test();
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("LOG_NORMAL_TRUNCATED_AB_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void log_normal_truncated_ab_cdf_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_NORMAL_TRUNCATED_AB_CDF_TEST tests LOG_NORMAL_TRUNCATED_AB_CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 March 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a;
        double b;
        double cdf;
        int i;
        double mu;
        double pdf;
        int seed = 123456789;
        double sigma;
        double x;
        double x2;

        Console.WriteLine("");
        Console.WriteLine("LOG_NORMAL_TRUNCATED_AB_CDF_TEST");
        Console.WriteLine("  LOG_NORMAL_TRUNCATED_AB_CDF evaluates the Log Normal Truncated AB CDF;");
        Console.WriteLine("  LOG_NORMAL_TRUNCATED_AB_CDF_INV inverts the Log Normal Truncated AB CDF.");
        Console.WriteLine("  LOG_NORMAL_TRUNCATED_AB_PDF evaluates the Log Normal Truncated AB PDF;");

        mu = 0.5;
        sigma = 3.0;
        a = Math.Exp(mu);
        b = Math.Exp(mu + 2.0 * sigma);

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter MU =     " + mu + "");
        Console.WriteLine("  PDF parameter SIGMA =  " + sigma + "");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!PDF.log_normal_truncated_ab_check(mu, sigma, a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("LOG_NORMAL_TRUNCATED_AB_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            x = PDF.log_normal_truncated_ab_sample(mu, sigma, a, b, ref seed);
            pdf = PDF.log_normal_truncated_ab_pdf(x, mu, sigma, a, b);
            cdf = CDF.log_normal_truncated_ab_cdf(x, mu, sigma, a, b);
            x2 = CDF.log_normal_truncated_ab_cdf_inv(cdf, mu, sigma, a, b);

            Console.WriteLine("  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    private static void log_normal_truncated_ab_sample_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_NORMAL_TRUNCATED_AB_SAMPLE_TEST tests LOG_NORMAL_TRUNCATED_AB_SAMPLE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 March 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int SAMPLE_NUM = 1000;

        double a;
        double b;
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
        Console.WriteLine("LOG_NORMAL_TRUNCATED_AB_SAMPLE_TEST");
        Console.WriteLine("  LOG_NORMAL_TRUNCATED_AB_MEAN computes the Log Normal Truncated AB mean;");
        Console.WriteLine("  LOG_NORMAL_TRUNCATED_AB_SAMPLE samples the Log Normal Truncated AB distribution;");
        Console.WriteLine("  LOG_NORMAL_TRUNCATED_AB_VARIANCE computes the Log Normal Truncated AB variance;");

        mu = 0.5;
        sigma = 3.0;
        a = Math.Exp(mu);
        b = Math.Exp(mu + 2.0 * sigma);

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter MU =     " + mu + "");
        Console.WriteLine("  PDF parameter SIGMA =  " + sigma + "");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!PDF.log_normal_truncated_ab_check(mu, sigma, a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("LOG_NORMAL_TRUNCATED_AB_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        mean = PDF.log_normal_truncated_ab_mean(mu, sigma, a, b);
        variance = PDF.log_normal_truncated_ab_variance(mu, sigma, a, b);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = PDF.log_normal_truncated_ab_sample(mu, sigma, a, b, ref seed);
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