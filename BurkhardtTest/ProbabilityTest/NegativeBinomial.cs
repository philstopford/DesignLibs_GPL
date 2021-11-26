using System;
using System.Globalization;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void negative_binomial_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    NEGATIVE_BINOMIAL_CDF_TEST tests NEGATIVE_BINOMIAL_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 April 2016
//
//  Author:
//
//    John Burkardt
//
    {
        int i;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("NEGATIVE_BINOMIAL_CDF_TEST");
        Console.WriteLine("  NEGATIVE_BINOMIAL_CDF evaluates the Negative Binomial CDF;");
        Console.WriteLine("  NEGATIVE_BINOMIAL_CDF_INV inverts the Negative Binomial CDF.");
        Console.WriteLine("  NEGATIVE_BINOMIAL_PDF evaluates the Negative Binomial PDF;");

        const int a = 2;
        const double b = 0.25;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!NegativeBinomial.negative_binomial_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("NEGATIVE_BINOMIAL_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            int x = NegativeBinomial.negative_binomial_sample(a, b, ref seed);
            double pdf = NegativeBinomial.negative_binomial_pdf(x, a, b);
            double cdf = NegativeBinomial.negative_binomial_cdf(x, a, b);
            int x2 = NegativeBinomial.negative_binomial_cdf_inv(cdf, a, b);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    private static void negative_binomial_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    NEGATIVE_BINOMIAL_SAMPLE_TEST tests NEGATIVE_BINOMIAL_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 April 2016
//
//  Author:
//
//    John Burkardt
//
    {
        const int SAMPLE_NUM = 1000;

        int i;
        int seed = 123456789;
        int[] x = new int[SAMPLE_NUM];

        Console.WriteLine("");
        Console.WriteLine("NEGATIVE_BINOMIAL_SAMPLE_TEST");
        Console.WriteLine("  NEGATIVE_BINOMIAL_MEAN computes the Negative Binomial mean;");
        Console.WriteLine("  NEGATIVE_BINOMIAL_SAMPLE samples the Negative Binomial distribution;");
        Console.WriteLine("  NEGATIVE_BINOMIAL_VARIANCE computes the Negative Binomial variance;");

        const int a = 2;
        const double b = 0.75;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!NegativeBinomial.negative_binomial_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("NEGATIVE_BINOMIAL_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        double mean = NegativeBinomial.negative_binomial_mean(a, b);
        double variance = NegativeBinomial.negative_binomial_variance(a, b);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = NegativeBinomial.negative_binomial_sample(a, b, ref seed);
        }

        mean = typeMethods.i4vec_mean(SAMPLE_NUM, x);
        variance = typeMethods.i4vec_variance(SAMPLE_NUM, x);
        int xmax = typeMethods.i4vec_max(SAMPLE_NUM, x);
        int xmin = typeMethods.i4vec_min(SAMPLE_NUM, x);

        Console.WriteLine("");
        Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
        Console.WriteLine("  Sample mean =     " + mean + "");
        Console.WriteLine("  Sample variance = " + variance + "");
        Console.WriteLine("  Sample maximum =  " + xmax + "");
        Console.WriteLine("  Sample minimum =  " + xmin + "");

    }

}