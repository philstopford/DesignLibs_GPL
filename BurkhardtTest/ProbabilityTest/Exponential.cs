using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal partial class Program
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
        double cdf;
        int i;
        double pdf;
        int seed = 123456789;
        double x;
        double x2;

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
            x = Exponential.exponential_01_sample(ref seed);
            pdf = Exponential.exponential_01_pdf(x);
            cdf = Exponential.exponential_01_cdf(x);
            x2 = Exponential.exponential_01_cdf_inv(cdf);

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
        int SAMPLE_NUM = 1000;

        int i;
        double mean;
        int seed = 123456789;
        double variance;
        double[] x = new double [SAMPLE_NUM];
        double xmax;
        double xmin;

        Console.WriteLine("");
        Console.WriteLine("EXPONENTIAL_01_SAMPLE_TEST");
        Console.WriteLine("  EXPONENTIAL_01_MEAN computes the Exponential 01 mean;");
        Console.WriteLine("  EXPONENTIAL_01_SAMPLE samples the Exponential 01 distribution;");
        Console.WriteLine("  EXPONENTIAL_01_VARIANCE computes the Exponential 01 variance.");

        mean = Exponential.exponential_01_mean();
        variance = Exponential.exponential_01_variance();

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Exponential.exponential_01_sample(ref seed);
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
        double a;
        double b;
        double cdf;
        int i;
        double pdf;
        int seed = 123456789;
        double x;
        double x2;

        Console.WriteLine("");
        Console.WriteLine("EXPONENTIAL_CDF_TEST");
        Console.WriteLine("  EXPONENTIAL_CDF evaluates the Exponential CDF;");
        Console.WriteLine("  EXPONENTIAL_CDF_INV inverts the Exponential CDF.");
        Console.WriteLine("  EXPONENTIAL_PDF evaluates the Exponential PDF;");

        a = 1.0;
        b = 2.0;

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
            x = Exponential.exponential_sample(a, b, ref seed);
            pdf = Exponential.exponential_pdf(x, a, b);
            cdf = Exponential.exponential_cdf(x, a, b);
            x2 = Exponential.exponential_cdf_inv(cdf, a, b);

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
        int SAMPLE_NUM = 1000;

        double a;
        double b;
        int i;
        double mean;
        int seed = 123456789;
        double variance;
        double[] x = new double [SAMPLE_NUM];
        double xmax;
        double xmin;

        Console.WriteLine("");
        Console.WriteLine("EXPONENTIAL_SAMPLE_TEST");
        Console.WriteLine("  EXPONENTIAL_MEAN computes the Exponential mean;");
        Console.WriteLine("  EXPONENTIAL_SAMPLE samples the Exponential distribution;");
        Console.WriteLine("  EXPONENTIAL_VARIANCE computes the Exponential variance;");

        a = 1.0;
        b = 10.0;

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

        mean = Exponential.exponential_mean(a, b);
        variance = Exponential.exponential_variance(a, b);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Exponential.exponential_sample(a, b, ref seed);
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