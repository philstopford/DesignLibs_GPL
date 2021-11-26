using System;
using System.Globalization;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal static partial class Program
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
        int i;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("LOG_NORMAL_CDF_TEST");
        Console.WriteLine("  LOG_NORMAL_CDF evaluates the Log Normal CDF;");
        Console.WriteLine("  LOG_NORMAL_CDF_INV inverts the Log Normal CDF.");
        Console.WriteLine("  LOG_NORMAL_PDF evaluates the Log Normal PDF;");

        const double a = 10.0;
        const double b = 2.25;

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
            double x = LogNormal.log_normal_sample(a, b, ref seed);
            double pdf = LogNormal.log_normal_pdf(x, a, b);
            double cdf = LogNormal.log_normal_cdf(x, a, b);
            double x2 = LogNormal.log_normal_cdf_inv(cdf, a, b);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
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
        const int SAMPLE_NUM = 1000;

        int i;
        int seed = 123456789;
        double[] x = new double [SAMPLE_NUM];

        Console.WriteLine("");
        Console.WriteLine("LOG_NORMAL_SAMPLE_TEST");
        Console.WriteLine("  LOG_NORMAL_MEAN computes the Log Normal mean;");
        Console.WriteLine("  LOG_NORMAL_SAMPLE samples the Log Normal distribution;");
        Console.WriteLine("  LOG_NORMAL_VARIANCE computes the Log Normal variance;");

        const double a = 1.0;
        const double b = 2.0;

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

        double mean = LogNormal.log_normal_mean(a, b);
        double variance = LogNormal.log_normal_variance(a, b);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = LogNormal.log_normal_sample(a, b, ref seed);
        }

        mean = typeMethods.r8vec_mean(SAMPLE_NUM, x);
        variance = typeMethods.r8vec_variance(SAMPLE_NUM, x);
        double xmax = typeMethods.r8vec_max(SAMPLE_NUM, x);
        double xmin = typeMethods.r8vec_min(SAMPLE_NUM, x);

        Console.WriteLine("");
        Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
        Console.WriteLine("  Sample mean =     " + mean + "");
        Console.WriteLine("  Sample variance = " + variance + "");
        Console.WriteLine("  Sample maximum =  " + xmax + "");
        Console.WriteLine("  Sample minimum =  " + xmin + "");
    }

}