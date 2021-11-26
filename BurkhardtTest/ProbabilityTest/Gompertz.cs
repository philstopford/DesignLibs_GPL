using System;
using System.Globalization;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void gompertz_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    GOMPERTZ_CDF_TEST tests GOMPERTZ_CDF
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
        Console.WriteLine("GOMPERTZ_CDF_TEST");
        Console.WriteLine("  GOMPERTZ_CDF evaluates the Gompertz CDF;");
        Console.WriteLine("  GOMPERTZ_CDF_INV inverts the Gompertz CDF.");
        Console.WriteLine("  GOMPERTZ_PDF evaluates the Gompertz PDF;");

        const double a = 2.0;
        const double b = 3.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!Gompertz.gompertz_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("GOMPERTZ_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            double x = Gompertz.gompertz_sample(a, b, ref seed);
            double pdf = Gompertz.gompertz_pdf(x, a, b);
            double cdf = Gompertz.gompertz_cdf(x, a, b);
            double x2 = Gompertz.gompertz_cdf_inv(cdf, a, b);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

    }

    private static void gompertz_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    GOMPERTZ_SAMPLE_TEST tests GOMPERTZ_SAMPLE.
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
        Console.WriteLine("GOMPERTZ_SAMPLE_TEST");
        Console.WriteLine("  GOMPERTZ_MEAN computes the Gompertz mean;");
        Console.WriteLine("  GOMPERTZ_SAMPLE samples the Gompertz distribution;");
        Console.WriteLine("  GOMPERTZ_VARIANCE computes the Gompertz variance;");

        const double a = 2.0;
        const double b = 3.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!Gompertz.gompertz_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("GOMPERTZ_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Gompertz.gompertz_sample(a, b, ref seed);
        }

        double mean = typeMethods.r8vec_mean(SAMPLE_NUM, x);
        double variance = typeMethods.r8vec_variance(SAMPLE_NUM, x);
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