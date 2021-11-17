using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal partial class Program
{
    private static void frechet_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    FRECHET_CDF_TEST tests FRECHET_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 April 2016
//
//  Author:
//
//    John Burkardt
//
    {
        double alpha;
        double cdf;
        int i;
        double pdf;
        int seed = 123456789;
        double x;
        double x2;

        Console.WriteLine("");
        Console.WriteLine("FRECHET_CDF_TEST");
        Console.WriteLine("  FRECHET_CDF evaluates the Frechet CDF;");
        Console.WriteLine("  FRECHET_CDF_INV inverts the Frechet CDF.");
        Console.WriteLine("  FRECHET_PDF evaluates the Frechet PDF;");

        alpha = 3.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter ALPHA =  " + alpha + "");

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            x = Frechet.frechet_sample(alpha, ref seed);
            pdf = Frechet.frechet_pdf(x, alpha);
            cdf = Frechet.frechet_cdf(x, alpha);
            x2 = Frechet.frechet_cdf_inv(cdf, alpha);

            Console.WriteLine("  "
                              + x.ToString().PadLeft(12) + "  "
                              + pdf.ToString().PadLeft(12) + "  "
                              + cdf.ToString().PadLeft(12) + "  "
                              + x2.ToString().PadLeft(12) + "");
        }

    }

    private static void frechet_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    FRECHET_SAMPLE_TEST tests FRECHET_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 April 2016
//
//  Author:
//
//    John Burkardt
//
    {
        int SAMPLE_NUM = 1000;

        double alpha;
        int i;
        double mean;
        int seed = 123456789;
        double variance;
        double[] x = new double [SAMPLE_NUM];
        double xmax;
        double xmin;

        Console.WriteLine("");
        Console.WriteLine("FRECHET_SAMPLE_TEST");
        Console.WriteLine("  FRECHET_MEAN computes the Frechet mean;");
        Console.WriteLine("  FRECHET_SAMPLE samples the Frechet distribution;");
        Console.WriteLine("  FRECHET_VARIANCE computes the Frechet variance;");

        alpha = 3.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter ALPHA =  " + alpha + "");

        mean = Frechet.frechet_mean(alpha);
        variance = Frechet.frechet_variance(alpha);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Frechet.frechet_sample(alpha, ref seed);
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