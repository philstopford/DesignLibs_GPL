using System;
using System.Globalization;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void logistic_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    LOGISTIC_CDF_TEST tests LOGISTIC_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 April 2016
//
//  Author:
//
//    John Burkardt
//
    {
        int i;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("LOGISTIC_CDF_TEST");
        Console.WriteLine("  LOGISTIC_CDF evaluates the Logistic CDF;");
        Console.WriteLine("  LOGISTIC_CDF_INV inverts the Logistic CDF.");
        Console.WriteLine("  LOGISTIC_PDF evaluates the Logistic PDF;");

        const double a = 1.0;
        const double b = 2.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!Logistic.logistic_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("LOGISTIC_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            double x = Logistic.logistic_sample(a, b, ref seed);
            double pdf = Logistic.logistic_pdf(x, a, b);
            double cdf = Logistic.logistic_cdf(x, a, b);
            double x2 = Logistic.logistic_cdf_inv(cdf, a, b);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    private static void logistic_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    LOGISTIC_SAMPLE_TEST tests LOGISTIC_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 April 2016
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
        Console.WriteLine("LOGISTIC_SAMPLE_TEST");
        Console.WriteLine("  LOGISTIC_MEAN computes the Logistic mean;");
        Console.WriteLine("  LOGISTIC_SAMPLE samples the Logistic distribution;");
        Console.WriteLine("  LOGISTIC_VARIANCE computes the Logistic variance;");

        const double a = 2.0;
        const double b = 3.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!Logistic.logistic_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("LOGISTIC_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        double mean = Logistic.logistic_mean(a, b);
        double variance = Logistic.logistic_variance(a, b);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Logistic.logistic_sample(a, b, ref seed);
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