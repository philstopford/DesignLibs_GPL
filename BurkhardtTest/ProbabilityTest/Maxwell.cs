using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal partial class Program
{
    private static void maxwell_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    MAXWELL_CDF_TEST tests MAXWELL_CDF.
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
        double cdf;
        int i;
        double pdf;
        int seed = 123456789;
        double x;
        double x2;

        Console.WriteLine("");
        Console.WriteLine("MAXWELL_CDF_TEST8");
        Console.WriteLine("  MAXWELL_CDF evaluates the Maxwell CDF;");
        Console.WriteLine("  MAXWELL_CDF_INV inverts the Maxwell CDF.");
        Console.WriteLine("  MAXWELL_PDF evaluates the Maxwell PDF;");

        a = 2.0E+00;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =             " + a + "");

        if (!Maxwell.maxwell_check(a))
        {
            Console.WriteLine("");
            Console.WriteLine("MAXWELL_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            x = Maxwell.maxwell_sample(a, ref seed);
            pdf = Maxwell.maxwell_pdf(x, a);
            cdf = Maxwell.maxwell_cdf(x, a);
            x2 = Maxwell.maxwell_cdf_inv(cdf, a);

            Console.WriteLine("  "
                              + x.ToString().PadLeft(12) + "  "
                              + pdf.ToString().PadLeft(12) + "  "
                              + cdf.ToString().PadLeft(12) + "  "
                              + x2.ToString().PadLeft(12) + "");
        }
    }

    private static void maxwell_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    MAXWELL_SAMPLE_TEST tests MAXWELL_SAMPLE.
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
        int j;
        double mean;
        int seed = 123456789;
        double variance;
        double[] x = new double [SAMPLE_NUM];
        double xmax;
        double xmin;

        Console.WriteLine("");
        Console.WriteLine("MAXWELL_SAMPLE_TEST");
        Console.WriteLine("  MAXWELL_MEAN computes the Maxwell mean;");
        Console.WriteLine("  MAXWELL_SAMPLE samples the Maxwell distribution;");
        Console.WriteLine("  MAXWELL_VARIANCE computes the Maxwell variance.");

        a = 2.0E+00;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =             " + a + "");

        if (!Maxwell.maxwell_check(a))
        {
            Console.WriteLine("");
            Console.WriteLine("MAXWELL_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        mean = Maxwell.maxwell_mean(a);
        variance = Maxwell.maxwell_variance(a);

        Console.WriteLine("  PDF mean =                    " + mean + "");
        Console.WriteLine("  PDF variance =                " + variance + "");

        for (j = 0; j < SAMPLE_NUM; j++)
        {
            x[j] = Maxwell.maxwell_sample(a, ref seed);
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