using System;
using System.Globalization;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void planck_pdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    PLANCK_PDF_TEST tests PLANCK_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 April 2016
//
//  Author:
//
//    John Burkardt
//
    {
        int i;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("PLANCK_PDF_TEST");
        Console.WriteLine("  PLANCK_PDF evaluates the Planck PDF.");
        Console.WriteLine("  PLANCK_SAMPLE samples the Planck PDF.");

        const double a = 2.0E+00;
        const double b = 3.0E+00;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A = " + a + "");
        Console.WriteLine("  PDF parameter B = " + b + "");

        if (!Planck.planck_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("PLANCK_PDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            double x = Planck.planck_sample(a, b, ref seed);

            double pdf = Planck.planck_pdf(x, a, b);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    private static void planck_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    PLANCK_SAMPLE_TEST tests PLANCK_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 April 2016
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
        Console.WriteLine("PLANCK_SAMPLE_TEST");
        Console.WriteLine("  PLANCK_MEAN computes the Planck mean;");
        Console.WriteLine("  PLANCK_SAMPLE samples the Planck distribution;");
        Console.WriteLine("  PLANCK_VARIANCE computes the Planck variance;");

        const double a = 2.0;
        const double b = 3.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!Planck.planck_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("PLANCK_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        double mean = Planck.planck_mean(a, b);
        double variance = Planck.planck_variance(a, b);

        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Planck.planck_sample(a, b, ref seed);
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