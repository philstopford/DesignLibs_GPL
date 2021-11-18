using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal partial class Program
{
    private static void semicircular_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    SEMICIRCULAR_CDF_TEST tests SEMICIRCULAR_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 March 2016
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
        Console.WriteLine("SEMICIRCULAR_CDF_TEST");
        Console.WriteLine("  SEMICIRCULAR_CDF evaluates the Semicircular CDF;");
        Console.WriteLine("  SEMICIRCULAR_CDF_INV inverts the Semicircular CDF.");
        Console.WriteLine("  SEMICIRCULAR_PDF evaluates the Semicircular PDF;");

        a = 3.0;
        b = 2.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!Semicircular.semicircular_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("SEMICIRCULAR_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            x = Semicircular.semicircular_sample(a, b, ref seed);
            pdf = Semicircular.semicircular_pdf(x, a, b);
            cdf = Semicircular.semicircular_cdf(x, a, b);
            x2 = Semicircular.semicircular_cdf_inv(cdf, a, b);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    private static void semicircular_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    SEMICIRCULAR_SAMPLE_TEST tests SEMICIRCULAR_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 March 2016
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
        Console.WriteLine("SEMICIRCULAR_SAMPLE_TEST");
        Console.WriteLine("  SEMICIRCULAR_MEAN computes the Semicircular mean;");
        Console.WriteLine("  SEMICIRCULAR_SAMPLE samples the Semicircular distribution;");
        Console.WriteLine("  SEMICIRCULAR_VARIANCE computes the Semicircular variance;");

        a = 3.0;
        b = 2.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!Semicircular.semicircular_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("SEMICIRCULAR_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        mean = Semicircular.semicircular_mean(a, b);
        variance = Semicircular.semicircular_variance(a, b);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Semicircular.semicircular_sample(a, b, ref seed);
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