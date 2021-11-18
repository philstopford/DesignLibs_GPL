using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal partial class Program
{
    private static void binomial_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    BINOMIAL_CDF_TEST tests BINOMIAL_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
    {
        int a;
        double b;
        double cdf;
        int i;
        double pdf;
        int seed = 123456789;
        int x;
        int x2;

        Console.WriteLine("");
        Console.WriteLine("BINOMIAL_CDF_TEST");
        Console.WriteLine("  BINOMIAL_CDF evaluates the Binomial CDF;");
        Console.WriteLine("  BINOMIAL_CDF_INV inverts the Binomial CDF.");
        Console.WriteLine("  BINOMIAL_PDF evaluates the Binomial PDF;");

        a = 5;
        b = 0.65;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!Binomial.binomial_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("BINOMIAL_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            x = Binomial.binomial_sample(a, b, ref seed);
            pdf = Binomial.binomial_pdf(x, a, b);
            cdf = Binomial.binomial_cdf(x, a, b);
            x2 = Binomial.binomial_cdf_inv(cdf, a, b);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

    }

    private static void binomial_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    BINOMIAL_SAMPLE_TEST tests BINOMIAL_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
    {
        int SAMPLE_NUM = 1000;

        int a;
        double b;
        int i;
        double mean;
        int seed = 123456789;
        double variance;
        int[] x = new int [SAMPLE_NUM];
        int xmax;
        int xmin;

        Console.WriteLine("");
        Console.WriteLine("BINOMIAL_SAMPLE_TEST");
        Console.WriteLine("  BINOMIAL_MEAN computes the Binomial mean;");
        Console.WriteLine("  BINOMIAL_SAMPLE samples the Binomial distribution;");
        Console.WriteLine("  BINOMIAL_VARIANCE computes the Binomial variance;");

        a = 5;
        b = 0.30;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!Binomial.binomial_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("BINOMIAL_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        mean = Binomial.binomial_mean(a, b);
        variance = Binomial.binomial_variance(a, b);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Binomial.binomial_sample(a, b, ref seed);
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