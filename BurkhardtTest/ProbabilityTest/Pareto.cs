using System;
using System.Globalization;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void pareto_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    PARETO_CDF_TEST tests PARETO_CDF.
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
        Console.WriteLine("PARETO_CDF_TEST");
        Console.WriteLine("  PARETO_CDF evaluates the Pareto CDF;");
        Console.WriteLine("  PARETO_CDF_INV inverts the Pareto CDF.");
        Console.WriteLine("  PARETO_PDF evaluates the Pareto PDF;");

        const double a = 0.5;
        const double b = 5.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!Pareto.pareto_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("PARETO_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            double x = Pareto.pareto_sample(a, b, ref seed);
            double pdf = Pareto.pareto_pdf(x, a, b);
            double cdf = Pareto.pareto_cdf(x, a, b);
            double x2 = Pareto.pareto_cdf_inv(cdf, a, b);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    private static void pareto_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    PARETO_SAMPLE_TEST tests PARETO_SAMPLE.
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
        Console.WriteLine("PARETO_SAMPLE_TEST");
        Console.WriteLine("  PARETO_MEAN computes the Pareto mean;");
        Console.WriteLine("  PARETO_SAMPLE samples the Pareto distribution;");
        Console.WriteLine("  PARETO_VARIANCE computes the Pareto variance;");

        const double a = 0.5;
        const double b = 5.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!Pareto.pareto_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("PARETO_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        double mean = Pareto.pareto_mean(a, b);
        double variance = Pareto.pareto_variance(a, b);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Pareto.pareto_sample(a, b, ref seed);
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