using System;
using System.Globalization;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void erlang_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    ERLANG_CDF_TEST tests ERLANG_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 March 2016
//
//  Author:
//
//    John Burkardt
//
    {
        int i;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("ERLANG_CDF_TEST");
        Console.WriteLine("  ERLANG_CDF evaluates the Erlang CDF;");
        Console.WriteLine("  ERLANG_CDF_INV inverts the Erlang CDF.");
        Console.WriteLine("  ERLANG_PDF evaluates the Erlang PDF;");

        const double a = 1.0;
        const double b = 2.0;
        const int c = 3;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");
        Console.WriteLine("  PDF parameter C =      " + c + "");

        if (!Erlang.erlang_check(a, b, c))
        {
            Console.WriteLine("");
            Console.WriteLine("ERLANG_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            double x = Erlang.erlang_sample(a, b, c, ref seed);
            double pdf = Erlang.erlang_pdf(x, a, b, c);
            double cdf = Erlang.erlang_cdf(x, a, b, c);
            double x2 = Erlang.erlang_cdf_inv(cdf, a, b, c);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

    }

    private static void erlang_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    ERLANG_SAMPLE_TEST tests ERLANG_SAMPLE.
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
        Console.WriteLine("ERLANG_SAMPLE_TEST");
        Console.WriteLine("  ERLANG_MEAN computes the Erlang mean;");
        Console.WriteLine("  ERLANG_SAMPLE samples the Erlang distribution;");
        Console.WriteLine("  ERLANG_VARIANCE computes the Erlang variance;");

        const double a = 1.0;
        const double b = 2.0;
        const int c = 3;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");
        Console.WriteLine("  PDF parameter C =      " + c + "");

        if (!Erlang.erlang_check(a, b, c))
        {
            Console.WriteLine("");
            Console.WriteLine("ERLANG_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        double mean = Erlang.erlang_mean(a, b, c);
        double variance = Erlang.erlang_variance(a, b, c);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Erlang.erlang_sample(a, b, c, ref seed);
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