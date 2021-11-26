using System;
using System.Globalization;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void von_mises_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    VON_MISES_CDF_TEST tests VON_MISES_CDF.
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
        int i;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("VON_MISES_CDF_TEST");
        Console.WriteLine("  VON_MISES_CDF evaluates the Von Mises CDF;");
        Console.WriteLine("  VON_MISES_CDF_INV inverts the Von Mises CDF.");
        Console.WriteLine("  VON_MISES_PDF evaluates the Von Mises PDF;");

        const double a = 1.0;
        const double b = 2.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!VonMises.von_mises_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("VON_MISES_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            double x = VonMises.von_mises_sample(a, b, ref seed);
            double pdf = VonMises.von_mises_pdf(x, a, b);
            double cdf = VonMises.von_mises_cdf(x, a, b);
            double x2 = VonMises.von_mises_cdf_inv(cdf, a, b);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    private static void von_mises_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    VON_MISES_SAMPLE_TEST tests VON_MISES_SAMPLE.
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
        const int SAMPLE_NUM = 1000;

        int i;
        int seed = 123456789;
        double[] x = new double [SAMPLE_NUM];

        Console.WriteLine("");
        Console.WriteLine("VON_MISES_SAMPLE_TEST");
        Console.WriteLine("  VON_MISES_MEAN computes the Von Mises mean;");
        Console.WriteLine("  VON_MISES_SAMPLE samples the Von Mises distribution;");
        Console.WriteLine("  VON_MISES_CIRCULAR_VARIANCE computes the Von Mises circular variance;");

        const double a = 1.0;
        const double b = 2.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!VonMises.von_mises_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("VON_MISES_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        double mean = VonMises.von_mises_mean(a, b);
        double variance = VonMises.von_mises_circular_variance(a, b);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =              " + mean + "");
        Console.WriteLine("  PDF circular variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = VonMises.von_mises_sample(a, b, ref seed);
        }

        mean = typeMethods.r8vec_mean(SAMPLE_NUM, x);
        variance = typeMethods.r8vec_circular_variance(SAMPLE_NUM, x);
        double xmax = typeMethods.r8vec_max(SAMPLE_NUM, x);
        double xmin = typeMethods.r8vec_min(SAMPLE_NUM, x);

        Console.WriteLine("");
        Console.WriteLine("  Sample size =              " + SAMPLE_NUM + "");
        Console.WriteLine("  Sample mean =              " + mean + "");
        Console.WriteLine("  Sample circular variance = " + variance + "");
        Console.WriteLine("  Sample maximum =           " + xmax + "");
        Console.WriteLine("  Sample minimum =           " + xmin + "");

    }

}