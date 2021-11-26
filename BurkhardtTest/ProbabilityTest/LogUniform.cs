using System;
using System.Globalization;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void log_uniform_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    LOG_UNIFORM_CDF_TEST tests LOG_UNIFORM_CDF.
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
        Console.WriteLine("LOG_UNIFORM_CDF_TEST");
        Console.WriteLine("  LOG_UNIFORM_CDF evaluates the Log Uniform CDF;");
        Console.WriteLine("  LOG_UNIFORM_CDF_INV inverts the Log Uniform CDF.");
        Console.WriteLine("  LOG_UNIFORM_PDF evaluates the Log Uniform PDF;");

        const double a = 2.0;
        const double b = 20.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!LogUniform.log_uniform_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("LOG_UNIFORM_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            double x = LogUniform.log_uniform_sample(a, b, ref seed);
            double pdf = LogUniform.log_uniform_pdf(x, a, b);
            double cdf = LogUniform.log_uniform_cdf(x, a, b);
            double x2 = LogUniform.log_uniform_cdf_inv(cdf, a, b);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

    }

    private static void log_uniform_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    LOG_UNIFORM_SAMPLE_TEST tests LOG_UNIFORM_SAMPLE;
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
        Console.WriteLine("LOG_UNIFORM_SAMPLE_TEST");
        Console.WriteLine("  LOG_UNIFORM_MEAN computes the Log Uniform mean;");
        Console.WriteLine("  LOG_UNIFORM_SAMPLE samples the Log Uniform distribution;");

        const double a = 2.0;
        const double b = 20.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!LogUniform.log_uniform_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("LOG_UNIFORM_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        double mean = LogUniform.log_uniform_mean(a, b);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = LogUniform.log_uniform_sample(a, b, ref seed);
        }

        mean = typeMethods.r8vec_mean(SAMPLE_NUM, x);
        double variance = typeMethods.r8vec_variance(SAMPLE_NUM, x);
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