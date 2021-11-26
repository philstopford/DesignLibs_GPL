using System;
using System.Globalization;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void gumbel_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    GUMBEL_CDF_TEST tests GUMBEL_CDF.
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
        Console.WriteLine("GUMBEL_CDF_TEST");
        Console.WriteLine("  GUMBEL_CDF evaluates the Gumbel CDF;");
        Console.WriteLine("  GUMBEL_CDF_INV inverts the Gumbel CDF.");
        Console.WriteLine("  GUMBEL_PDF evaluates the Gumbel PDF;");

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            double x = Gumbel.gumbel_sample(ref seed);
            double pdf = Gumbel.gumbel_pdf(x);
            double cdf = Gumbel.gumbel_cdf(x);
            double x2 = Gumbel.gumbel_cdf_inv(cdf);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

    }

    private static void gumbel_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    GUMBEL_SAMPLE_TEST tests GUMBEL_SAMPLE.
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
        double[] x = new double [SAMPLE_NUM];

        Console.WriteLine("");
        Console.WriteLine("GUMBEL_SAMPLE_TEST");
        Console.WriteLine("  GUMBEL_MEAN computes the Gumbel mean;");
        Console.WriteLine("  GUMBEL_SAMPLE samples the Gumbel distribution;");
        Console.WriteLine("  GUMBEL_VARIANCE computes the Gumbel variance.");

        double mean = Gumbel.gumbel_mean();
        double variance = Gumbel.gumbel_variance();

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Gumbel.gumbel_sample(ref seed);
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