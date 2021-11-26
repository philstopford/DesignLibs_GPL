using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void genlogistic_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    GENLOGISTIC_CDF_TEST tests GENLOGISTIC_CDF.
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
        double a;
        double b;
        double c;
        double cdf;
        int i;
        double pdf;
        int seed = 123456789;
        double x;
        double x2;

        Console.WriteLine("");
        Console.WriteLine("GENLOGISTIC_CDF_TEST");
        Console.WriteLine("  GENLOGISTIC_CDF evaluates the Genlogistic CDF;");
        Console.WriteLine("  GENLOGISTIC_CDF_INV inverts the Genlogistic CDF.");
        Console.WriteLine("  GENLOGISTIC_PDF evaluates the Genlogistic PDF;");

        a = 1.0;
        b = 2.0;
        c = 3.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");
        Console.WriteLine("  PDF parameter C =      " + c + "");

        if (!Genlogistic.genlogistic_check(a, b, c))
        {
            Console.WriteLine("");
            Console.WriteLine("GENLOGISTIC_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            x = Genlogistic.genlogistic_sample(a, b, c, ref seed);
            pdf = Genlogistic.genlogistic_pdf(x, a, b, c);
            cdf = Genlogistic.genlogistic_cdf(x, a, b, c);
            x2 = Genlogistic.genlogistic_cdf_inv(cdf, a, b, c);

            Console.WriteLine("  "
                              + x.ToString().PadLeft(12) + "  "
                              + pdf.ToString().PadLeft(12) + "  "
                              + cdf.ToString().PadLeft(12) + "  "
                              + x2.ToString().PadLeft(12) + "");
        }
    }

    private static void genlogistic_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    GENLOGISTIC_SAMPLE_TEST tests GENLOGISTIC_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2007
//
//  Author:
//
//    John Burkardt
//
    {
        int SAMPLE_NUM = 1000;

        double a;
        double b;
        double c;
        int i;
        double mean;
        int seed = 123456789;
        double variance;
        double[] x = new double [SAMPLE_NUM];
        double xmax;
        double xmin;

        Console.WriteLine("");
        Console.WriteLine("GENLOGISTIC_SAMPLE_TEST");
        Console.WriteLine("  GENLOGISTIC_MEAN computes the Genlogistic mean;");
        Console.WriteLine("  GENLOGISTIC_SAMPLE samples the Genlogistic distribution;");
        Console.WriteLine("  GENLOGISTIC_VARIANCE computes the Genlogistic variance;");

        a = 1.0;
        b = 2.0;
        c = 3.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");
        Console.WriteLine("  PDF parameter C =      " + c + "");

        if (!Genlogistic.genlogistic_check(a, b, c))
        {
            Console.WriteLine("");
            Console.WriteLine("GENLOGISTIC_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        mean = Genlogistic.genlogistic_mean(a, b, c);
        variance = Genlogistic.genlogistic_variance(a, b, c);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Genlogistic.genlogistic_sample(a, b, c, ref seed);
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