using System;
using System.Globalization;
using Burkardt.Types;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void power_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    POWER_CDF_TEST tests POWER_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 March 2016
//
//  Author:
//
//    John Burkardt
//
    {
        int i;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("POWER_CDF_TEST");
        Console.WriteLine("  POWER_CDF evaluates the Power CDF;");
        Console.WriteLine("  POWER_CDF_INV inverts the Power CDF.");
        Console.WriteLine("  POWER_PDF evaluates the Power PDF;");

        double a = 2.0;
        double b = 3.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!Burkardt.Probability.Power.power_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("POWER_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            double x = Burkardt.Probability.Power.power_sample(a, b, ref seed);
            double pdf = Burkardt.Probability.Power.power_pdf(x, a, b);
            double cdf = Burkardt.Probability.Power.power_cdf(x, a, b);
            double x2 = Burkardt.Probability.Power.power_cdf_inv(cdf, a, b);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    private static void power_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    POWER_SAMPLE_TEST tests POWER_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 March 2016
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
        Console.WriteLine("POWER_SAMPLE_TEST");
        Console.WriteLine("  POWER_MEAN computes the Power mean;");
        Console.WriteLine("  POWER_SAMPLE samples the Power distribution;");
        Console.WriteLine("  POWER_VARIANCE computes the Power variance;");

        const double a = 2.0;
        const double b = 3.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!Burkardt.Probability.Power.power_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("POWER_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        double mean = Burkardt.Probability.Power.power_mean(a, b);
        double variance = Burkardt.Probability.Power.power_variance(a, b);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Burkardt.Probability.Power.power_sample(a, b, ref seed);
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

internal class Probability
{
}