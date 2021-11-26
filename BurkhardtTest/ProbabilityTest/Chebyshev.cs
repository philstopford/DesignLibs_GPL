using System;
using System.Globalization;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void chebyshev1_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    CHEBYSHEV1_CDF_TEST tests CHEBYSHEV1_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 August 2016
//
//  Author:
//
//    John Burkardt
//
    {
        int i;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("CHEBYSHEV1_CDF_TEST");

        Console.WriteLine("  CHEBYSHEV1_CDF evaluates the Chebyshev1 CDF;");
        Console.WriteLine("  CHEBYSHEV1_CDF_INV inverts the Chebyshev1 CDF.");
        Console.WriteLine("  CHEBYSHEV1_PDF evaluates the Chebyshev1 PDF;");

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            double x = Chebyshevi.chebyshev1_sample(ref seed);
            double pdf = Chebyshevi.chebyshev1_pdf(x);
            double cdf = Chebyshevi.chebyshev1_cdf(x);
            double x2 = Chebyshevi.chebyshev1_cdf_inv(cdf);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

    }

    private static void chebyshev1_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    CHEBYSHEV1_SAMPLE_TEST tests CHEBYSHEV1_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 August 2016
//
//  Author:
//
//    John Burkardt
//
    {
        const int SAMPLE_NUM = 1000;

        int i;
        int seed = 123456789;
        double[] x = new double[SAMPLE_NUM];

        Console.WriteLine("");
        Console.WriteLine("CHEBYSHEV1_SAMPLE_TEST");
        Console.WriteLine("  CHEBYSHEV1_MEAN computes the Chebyshev1 mean;");
        Console.WriteLine("  CHEBYSHEV1_SAMPLE samples the Chebyshev1 distribution;");
        Console.WriteLine("  CHEBYSHEV1_VARIANCE computes the Chebyshev1 variance.");

        double mean = Chebyshevi.chebyshev1_mean();
        double variance = Chebyshevi.chebyshev1_variance();

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Chebyshevi.chebyshev1_sample(ref seed);
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