using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal partial class Program
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
        double a;
        double b;
        double cdf;
        int i;
        double pdf;
        int seed = 123456789;
        double x;
        double x2;

        Console.WriteLine("");
        Console.WriteLine("LAPLACE_CDF_TEST");
        Console.WriteLine("  LAPLACE_CDF evaluates the Laplace CDF;");
        Console.WriteLine("  LAPLACE_CDF_INV inverts the Laplace CDF.");
        Console.WriteLine("  LAPLACE_PDF evaluates the Laplace PDF;");

        a = 1.0;
        b = 2.0;

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
            x = Laplace.laplace_sample(a, b, ref seed);
            pdf = Laplace.laplace_pdf(x, a, b);
            cdf = Laplace.laplace_cdf(x, a, b);
            x2 = Laplace.laplace_cdf_inv(cdf, a, b);

            Console.WriteLine("  "
                              + x.ToString().PadLeft(12) + "  "
                              + pdf.ToString().PadLeft(12) + "  "
                              + cdf.ToString().PadLeft(12) + "  "
                              + x2.ToString().PadLeft(12) + "");
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
        Console.WriteLine("LAPLACE_SAMPLE_TEST");
        Console.WriteLine("  LAPLACE_MEAN computes the Laplace mean;");
        Console.WriteLine("  LAPLACE_SAMPLE samples the Laplace distribution;");
        Console.WriteLine("  LAPLACE_VARIANCE computes the Laplace variance;");

        a = 1.0;
        b = 2.0;

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

        mean = Laplace.laplace_mean(a, b);
        variance = Laplace.laplace_variance(a, b);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Laplace.laplace_sample(a, b, ref seed);
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