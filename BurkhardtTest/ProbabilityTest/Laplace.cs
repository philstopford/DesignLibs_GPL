using System;
using System.Globalization;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void laplace_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    LAPLACE_CDF_TEST tests LAPLACE_CDF.
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
        Console.WriteLine("LAPLACE_CDF_TEST");
        Console.WriteLine("  LAPLACE_CDF evaluates the Laplace CDF;");
        Console.WriteLine("  LAPLACE_CDF_INV inverts the Laplace CDF.");
        Console.WriteLine("  LAPLACE_PDF evaluates the Laplace PDF;");

        const double a = 1.0;
        const double b = 2.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!Laplace.laplace_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("LAPLACE_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            double x = Laplace.laplace_sample(a, b, ref seed);
            double pdf = Laplace.laplace_pdf(x, a, b);
            double cdf = Laplace.laplace_cdf(x, a, b);
            double x2 = Laplace.laplace_cdf_inv(cdf, a, b);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    private static void laplace_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    LAPLACE_SAMPLE_TEST tests LAPLACE_SAMPLE.
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
        Console.WriteLine("LAPLACE_SAMPLE_TEST");
        Console.WriteLine("  LAPLACE_MEAN computes the Laplace mean;");
        Console.WriteLine("  LAPLACE_SAMPLE samples the Laplace distribution;");
        Console.WriteLine("  LAPLACE_VARIANCE computes the Laplace variance;");

        double a = 1.0;
        double b = 2.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!Laplace.laplace_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("LAPLACE_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        double mean = Laplace.laplace_mean(a, b);
        double variance = Laplace.laplace_variance(a, b);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Laplace.laplace_sample(a, b, ref seed);
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