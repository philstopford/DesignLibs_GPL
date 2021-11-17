using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal partial class Program
{
    private static void gamma_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_CDF_TEST tests GAMMA_CDF.
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
        double a;
        double b;
        double c;
        double cdf;
        int i;
        double pdf;
        int seed = 123456789;
        double x;

        Console.WriteLine("");
        Console.WriteLine("GAMMA_CDF_TEST");
        Console.WriteLine("  GAMMA_CDF evaluates the Gamma CDF;");
        Console.WriteLine("  GAMMA_PDF evaluates the Gamma PDF;");
        Console.WriteLine("  GAMMA_SAMPLE samples the Gamma PDF;");

        a = 1.0;
        b = 1.5;
        c = 3.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A = " + a + "");
        Console.WriteLine("  PDF parameter B = " + b + "");
        Console.WriteLine("  PDF parameter B = " + c + "");

        if (!Gamma.gamma_check(a, b, c))
        {
            Console.WriteLine("");
            Console.WriteLine("GAMMA_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameter values are illegal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            x = Gamma.gamma_sample(a, b, c, ref seed);
            pdf = Gamma.gamma_pdf(x, a, b, c);
            cdf = Gamma.gamma_cdf(x, a, b, c);

            Console.WriteLine("  "
                              + x.ToString().PadLeft(12) + "  "
                              + pdf.ToString().PadLeft(12) + "  "
                              + cdf.ToString().PadLeft(12) + "");
        }
    }

    private static void gamma_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_SAMPLE_TEST tests GAMMA_SAMPLE.
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
        Console.WriteLine("GAMMA_SAMPLE_TEST");
        Console.WriteLine("  GAMMA_MEAN computes the Gamma mean;");
        Console.WriteLine("  GAMMA_SAMPLE samples the Gamma distribution;");
        Console.WriteLine("  GAMMA_VARIANCE computes the Gamma variance;");

        a = 1.0;
        b = 3.0;
        c = 2.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");
        Console.WriteLine("  PDF parameter C =      " + c + "");

        if (!Gamma.gamma_check(a, b, c))
        {
            Console.WriteLine("");
            Console.WriteLine("GAMMA_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        mean = Gamma.gamma_mean(a, b, c);
        variance = Gamma.gamma_variance(a, b, c);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Gamma.gamma_sample(a, b, c, ref seed);
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