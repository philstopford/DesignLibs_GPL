using System;
using System.Globalization;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void triangle_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_CDF_TEST tests TRIANGLE_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
    {
        int i;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("TRIANGLE_CDF_TEST");
        Console.WriteLine("  TRIANGLE_CDF evaluates the Triangle CDF;");
        Console.WriteLine("  TRIANGLE_CDF_INV inverts the Triangle CDF.");
        Console.WriteLine("  TRIANGLE_PDF evaluates the Triangle PDF;");

        const double a = 1.0;
        const double b = 3.0;
        const double c = 10.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");
        Console.WriteLine("  PDF parameter C =      " + c + "");

        if (!Triangle.triangle_check(a, b, c))
        {
            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            double x = Triangle.triangle_sample(a, b, c, ref seed);
            double pdf = Triangle.triangle_pdf(x, a, b, c);
            double cdf = Triangle.triangle_cdf(x, a, b, c);
            double x2 = Triangle.triangle_cdf_inv(cdf, a, b, c);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

    }

    private static void triangle_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_SAMPLE_TEST tests TRIANGLE_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 March 2016
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
        Console.WriteLine("TRIANGLE_SAMPLE_TEST");
        Console.WriteLine("  TRIANGLE_MEAN computes the Triangle mean;");
        Console.WriteLine("  TRIANGLE_SAMPLE samples the Triangle distribution;");
        Console.WriteLine("  TRIANGLE_VARIANCE computes the Triangle variance;");

        const double a = 1.0;
        const double b = 3.0;
        const double c = 10.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =        " + a + "");
        Console.WriteLine("  PDF parameter B =        " + b + "");
        Console.WriteLine("  PDF parameter C =        " + c + "");

        if (!Triangle.triangle_check(a, b, c))
        {
            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        double mean = Triangle.triangle_mean(a, b, c);
        double variance = Triangle.triangle_variance(a, b, c);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Triangle.triangle_sample(a, b, c, ref seed);
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