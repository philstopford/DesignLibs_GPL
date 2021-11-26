using System;
using System.Globalization;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void burr_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    BURR_CDF_TEST tests BURR_CDF, BURR_CDF_INV, BURR_PDF;
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
        int i;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("BURR_CDF_TEST");
        Console.WriteLine("  BURR_CDF evaluates the Burr CDF;");
        Console.WriteLine("  BURR_CDF_INV inverts the Burr CDF.");
        Console.WriteLine("  BURR_PDF evaluates the Burr PDF;");

        const double a = 1.0;
        const double b = 2.0;
        const double c = 3.0;
        const double d = 2.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");
        Console.WriteLine("  PDF parameter C =      " + c + "");
        Console.WriteLine("  PDF parameter D =      " + d + "");

        if (!Burr.burr_check(a, b, c, d))
        {
            Console.WriteLine("");
            Console.WriteLine("BURR_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            double x = Burr.burr_sample(a, b, c, d, ref seed);
            double pdf = Burr.burr_pdf(x, a, b, c, d);
            double cdf = Burr.burr_cdf(x, a, b, c, d);
            double x2 = Burr.burr_cdf_inv(cdf, a, b, c, d);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

    }

    private static void burr_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    BURR_SAMPLE_TEST tests BURR_MEAN, BURR_SAMPLE, BURR_VARIANCE;
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
        const int SAMPLE_NUM = 1000;

        int i;
        int seed = 123456789;
        double[] x = new double [SAMPLE_NUM];

        Console.WriteLine("");
        Console.WriteLine("BURR_SAMPLE_TEST");
        Console.WriteLine("  BURR_MEAN computes the Burr mean;");
        Console.WriteLine("  BURR_SAMPLE samples the Burr distribution;");
        Console.WriteLine("  BURR_VARIANCE computes the Burr variance;");

        const double a = 1.0;
        const double b = 2.0;
        const double c = 3.0;
        const double d = 2.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");
        Console.WriteLine("  PDF parameter C =      " + c + "");
        Console.WriteLine("  PDF parameter D =      " + d + "");

        if (!Burr.burr_check(a, b, c, d))
        {
            Console.WriteLine("");
            Console.WriteLine("BURR_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        double mean = Burr.burr_mean(a, b, c, d);
        double variance = Burr.burr_variance(a, b, c, d);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Burr.burr_sample(a, b, c, d, ref seed);
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