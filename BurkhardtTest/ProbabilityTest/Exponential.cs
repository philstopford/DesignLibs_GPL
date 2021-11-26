using System;
using System.Globalization;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void exponential_01_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    EXPONENTIAL_01_CDF_TEST tests EXPONENTIAL_01_CDF.
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
        Console.WriteLine("EXPONENTIAL_01_CDF_TEST");
        Console.WriteLine("  EXPONENTIAL_01_CDF evaluates the Exponential 01 CDF;");
        Console.WriteLine("  EXPONENTIAL_01_CDF_INV inverts the Exponential 01 CDF.");
        Console.WriteLine("  EXPONENTIAL_01_PDF evaluates the Exponential 01 PDF;");

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            double x = Exponential.exponential_01_sample(ref seed);
            double pdf = Exponential.exponential_01_pdf(x);
            double cdf = Exponential.exponential_01_cdf(x);
            double x2 = Exponential.exponential_01_cdf_inv(cdf);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    private static void exponential_01_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    EXPONENTIAL_01_SAMPLE_TEST tests EXPONENTIAL_01_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 March 2016
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
        Console.WriteLine("EXPONENTIAL_01_SAMPLE_TEST");
        Console.WriteLine("  EXPONENTIAL_01_MEAN computes the Exponential 01 mean;");
        Console.WriteLine("  EXPONENTIAL_01_SAMPLE samples the Exponential 01 distribution;");
        Console.WriteLine("  EXPONENTIAL_01_VARIANCE computes the Exponential 01 variance.");

        double mean = Exponential.exponential_01_mean();
        double variance = Exponential.exponential_01_variance();

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Exponential.exponential_01_sample(ref seed);
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

    private static void exponential_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    EXPONENTIAL_CDF_TEST tests EXPONENTIAL_CDF.
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
        Console.WriteLine("EXPONENTIAL_CDF_TEST");
        Console.WriteLine("  EXPONENTIAL_CDF evaluates the Exponential CDF;");
        Console.WriteLine("  EXPONENTIAL_CDF_INV inverts the Exponential CDF.");
        Console.WriteLine("  EXPONENTIAL_PDF evaluates the Exponential PDF;");

        const double a = 1.0;
        const double b = 2.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!Exponential.exponential_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("EXPONENTIAL_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            double x = Exponential.exponential_sample(a, b, ref seed);
            double pdf = Exponential.exponential_pdf(x, a, b);
            double cdf = Exponential.exponential_cdf(x, a, b);
            double x2 = Exponential.exponential_cdf_inv(cdf, a, b);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    private static void exponential_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    EXPONENTIAL_SAMPLE_TEST tests EXPONENTIAL_SAMPLE.
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
        const int SAMPLE_NUM = 1000;

        int i;
        int seed = 123456789;
        double[] x = new double [SAMPLE_NUM];

        Console.WriteLine("");
        Console.WriteLine("EXPONENTIAL_SAMPLE_TEST");
        Console.WriteLine("  EXPONENTIAL_MEAN computes the Exponential mean;");
        Console.WriteLine("  EXPONENTIAL_SAMPLE samples the Exponential distribution;");
        Console.WriteLine("  EXPONENTIAL_VARIANCE computes the Exponential variance;");

        const double a = 1.0;
        const double b = 10.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!Exponential.exponential_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("EXPONENTIAL_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        double mean = Exponential.exponential_mean(a, b);
        double variance = Exponential.exponential_variance(a, b);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Exponential.exponential_sample(a, b, ref seed);
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