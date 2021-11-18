using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal partial class Program
{
    private static void f_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    F_CDF_TEST tests F_CDF.
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
        int m;
        int n;
        double pdf;
        int seed = 123456789;
        double x;

        Console.WriteLine("");
        Console.WriteLine("F_CDF_TEST");
        Console.WriteLine("  F_CDF evaluates the F CDF;");
        Console.WriteLine("  F_PDF evaluates the F PDF;");
        Console.WriteLine("  F_SAMPLE samples the F PDF;");

        m = 1;
        n = 1;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter M = " + m + "");
        Console.WriteLine("  PDF parameter N = " + n + "");

        if (!FProb.f_check(m, n))
        {
            Console.WriteLine("");
            Console.WriteLine("F_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameter values are illegal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            x = FProb.f_sample(m, n, ref seed);
            pdf = FProb.f_pdf(x, m, n);
            cdf = FProb.f_cdf(x, m, n);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    private static void f_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    F_SAMPLE_TEST tests F_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 April 2016
//
//  Author:
//
//    John Burkardt
//
    {
        int SAMPLE_NUM = 1000;

        int i;
        int m;
        double mean;
        int n;
        int seed = 123456789;
        double variance;
        double[] x = new double [SAMPLE_NUM];
        double xmax;
        double xmin;

        Console.WriteLine("");
        Console.WriteLine("F_SAMPLE_TEST");
        Console.WriteLine("  F_MEAN computes the F mean;");
        Console.WriteLine("  F_SAMPLE samples the F distribution;");
        Console.WriteLine("  F_VARIANCE computes the F variance;");

        m = 8;
        n = 6;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter M = " + m + "");
        Console.WriteLine("  PDF parameter N = " + n + "");

        if (!FProb.f_check(m, n))
        {
            Console.WriteLine("");
            Console.WriteLine("F_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        mean = FProb.f_mean(m, n);
        variance = FProb.f_variance(m, n);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = FProb.f_sample(m, n, ref seed);
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