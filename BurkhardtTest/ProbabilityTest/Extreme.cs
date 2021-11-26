using System;
using System.Globalization;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void extreme_values_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    EXTREME_VALUES_CDF_TEST tests EXTREME_VALUES_CDF.
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
        Console.WriteLine("EXTREME_VALUES_CDF_TEST");
        Console.WriteLine("  EXTREME_VALUES_CDF evaluates the Extreme Values CDF;");
        Console.WriteLine("  EXTREME_VALUES_CDF_INV inverts the Extreme Values CDF.");
        Console.WriteLine("  EXTREME_VALUES_PDF evaluates the Extreme Values PDF;");

        const double a = 2.0;
        const double b = 3.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!Extreme.extreme_values_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("EXTREME_VALUES_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            double x = Extreme.extreme_values_sample(a, b, ref seed);
            double pdf = Extreme.extreme_values_pdf(x, a, b);
            double cdf = Extreme.extreme_values_cdf(x, a, b);
            double x2 = Extreme.extreme_values_cdf_inv(cdf, a, b);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    private static void extreme_values_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    EXTREME_VALUES_SAMPLE_TEST tests EXTREME_VALUES_SAMPLE.
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
        const int SAMPLE_NUM = 1000;

        int i;
        int seed = 123456789;
        double[] x = new double [SAMPLE_NUM];

        Console.WriteLine("");
        Console.WriteLine("EXTREME_VALUES_SAMPLE_TEST");
        Console.WriteLine("  EXTREME_VALUES_MEAN computes the Extreme Values mean;");
        Console.WriteLine("  EXTREME_VALUES_SAMPLE samples the Extreme Values distribution;");
        Console.WriteLine("  EXTREME_VALUES_VARIANCE computes the Extreme Values variance;");

        const double a = 2.0;
        const double b = 3.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!Extreme.extreme_values_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("EXTREME_VALUES_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        double mean = Extreme.extreme_values_mean(a, b);
        double variance = Extreme.extreme_values_variance(a, b);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Extreme.extreme_values_sample(a, b, ref seed);
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