using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal partial class Program
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
        double a;
        double b;
        double cdf;
        int i;
        double pdf;
        int seed = 123456789;
        double x;
        double x2;

        Console.WriteLine("");
        Console.WriteLine("LOGISTIC_CDF_TEST");
        Console.WriteLine("  LOGISTIC_CDF evaluates the Logistic CDF;");
        Console.WriteLine("  LOGISTIC_CDF_INV inverts the Logistic CDF.");
        Console.WriteLine("  LOGISTIC_PDF evaluates the Logistic PDF;");

        a = 1.0;
        b = 2.0;

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
            x = Logistic.logistic_sample(a, b, ref seed);
            pdf = Logistic.logistic_pdf(x, a, b);
            cdf = Logistic.logistic_cdf(x, a, b);
            x2 = Logistic.logistic_cdf_inv(cdf, a, b);

            Console.WriteLine("  "
                              + x.ToString().PadLeft(12) + "  "
                              + pdf.ToString().PadLeft(12) + "  "
                              + cdf.ToString().PadLeft(12) + "  "
                              + x2.ToString().PadLeft(12) + "");
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
        Console.WriteLine("LOGISTIC_SAMPLE_TEST");
        Console.WriteLine("  LOGISTIC_MEAN computes the Logistic mean;");
        Console.WriteLine("  LOGISTIC_SAMPLE samples the Logistic distribution;");
        Console.WriteLine("  LOGISTIC_VARIANCE computes the Logistic variance;");

        a = 2.0;
        b = 3.0;

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

        mean = Logistic.logistic_mean(a, b);
        variance = Logistic.logistic_variance(a, b);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Logistic.logistic_sample(a, b, ref seed);
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