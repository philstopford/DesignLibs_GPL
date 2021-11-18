using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal partial class Program
{
    private static void reciprocal_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    RECIPROCAL_CDF_TEST tests RECIPROCAL_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 March 2016
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
        Console.WriteLine("RECIPROCAL_CDF_TEST");
        Console.WriteLine("  RECIPROCAL_CDF evaluates the Reciprocal CDF;");
        Console.WriteLine("  RECIPROCAL_CDF_INV inverts the Reciprocal CDF.");
        Console.WriteLine("  RECIPROCAL_PDF evaluates the Reciprocal PDF;");

        a = 1.0;
        b = 3.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!Reciprocal.reciprocal_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("RECIPROCAL_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            x = Reciprocal.reciprocal_sample(a, b, ref seed);
            pdf = Reciprocal.reciprocal_pdf(x, a, b);
            cdf = Reciprocal.reciprocal_cdf(x, a, b);
            x2 = Reciprocal.reciprocal_cdf_inv(cdf, a, b);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    private static void reciprocal_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    RECIPROCAL_SAMPLE_TEST tests RECIPROCAL_SAMPLE.
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
        Console.WriteLine("RECIPROCAL_SAMPLE_TEST");
        Console.WriteLine("  RECIPROCAL_MEAN computes the Reciprocal mean;");
        Console.WriteLine("  RECIPROCAL_SAMPLE samples the Reciprocal distribution;");
        Console.WriteLine("  RECIPROCAL_VARIANCE computes the Reciprocal variance;");

        a = 1.0;
        b = 3.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!Reciprocal.reciprocal_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("RECIPROCAL_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        mean = Reciprocal.reciprocal_mean(a, b);
        variance = Reciprocal.reciprocal_variance(a, b);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Reciprocal.reciprocal_sample(a, b, ref seed);
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