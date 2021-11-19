using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal partial class Program
{
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
        double a;
        double b;
        double cdf;
        int i;
        double pdf;
        int seed = 123456789;
        double x;
        double x2;

        Console.WriteLine("");
        Console.WriteLine("LOG_NORMAL_CDF_TEST");
        Console.WriteLine("  LOG_NORMAL_CDF evaluates the Log Normal CDF;");
        Console.WriteLine("  LOG_NORMAL_CDF_INV inverts the Log Normal CDF.");
        Console.WriteLine("  LOG_NORMAL_PDF evaluates the Log Normal PDF;");

        a = 10.0;
        b = 2.25;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!LogNormal.log_normal_check(a, b))
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
            x = LogNormal.log_normal_sample(a, b, ref seed);
            pdf = LogNormal.log_normal_pdf(x, a, b);
            cdf = LogNormal.log_normal_cdf(x, a, b);
            x2 = LogNormal.log_normal_cdf_inv(cdf, a, b);

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

        double a;
        double b;
        int i;
        double mean;
        int seed = 123456789;
        double variance;
        double[] x = new double [SAMPLE_NUM];
        double xmax;
        double xmin;

        Console.WriteLine("");
        Console.WriteLine("LOG_NORMAL_SAMPLE_TEST");
        Console.WriteLine("  LOG_NORMAL_MEAN computes the Log Normal mean;");
        Console.WriteLine("  LOG_NORMAL_SAMPLE samples the Log Normal distribution;");
        Console.WriteLine("  LOG_NORMAL_VARIANCE computes the Log Normal variance;");

        a = 1.0;
        b = 2.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!LogNormal.log_normal_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("LOG_NORMAL_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        mean = LogNormal.log_normal_mean(a, b);
        variance = LogNormal.log_normal_variance(a, b);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = LogNormal.log_normal_sample(a, b, ref seed);
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