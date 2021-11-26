using System;
using System.Globalization;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void cosine_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    COSINE_CDF_TEST tests COSINE_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
    {
        int i;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("COSINE_CDF_TEST");
        Console.WriteLine("  COSINE_CDF evaluates the Cosine CDF;");
        Console.WriteLine("  COSINE_CDF_INV inverts the Cosine CDF.");
        Console.WriteLine("  COSINE_PDF evaluates theCosine  PDF;");

        const double a = 2.0;
        const double b = 1.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!Cosine.cosine_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("COSINE_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            double x = Cosine.cosine_sample(a, b, ref seed);
            double pdf = Cosine.cosine_pdf(x, a, b);
            double cdf = Cosine.cosine_cdf(x, a, b);
            double x2 = Cosine.cosine_cdf_inv(cdf, a, b);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    private static void cosine_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    COSINE_SAMPLE_TEST tests COSINE_SAMPLE.
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
        const int SAMPLE_NUM = 1000;

        int i;
        int seed = 123456789;
        double[] x = new double [SAMPLE_NUM];

        Console.WriteLine("");
        Console.WriteLine("COSINE_SAMPLE_TEST");
        Console.WriteLine("  COSINE_MEAN computes the Cosine mean;");
        Console.WriteLine("  COSINE_SAMPLE samples the Cosine distribution;");
        Console.WriteLine("  COSINE_VARIANCE computes the Cosine variance;");

        const double a = 2.0;
        const double b = 1.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!Cosine.cosine_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("COSINE_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        double mean = Cosine.cosine_mean(a, b);
        double variance = Cosine.cosine_variance(a, b);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Cosine.cosine_sample(a, b, ref seed);
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