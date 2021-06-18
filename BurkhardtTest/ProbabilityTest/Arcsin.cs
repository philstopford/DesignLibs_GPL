using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest
{
    partial class Program
    {
    static void arcsin_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    ARCSIN_CDF_TEST tests ARCSIN_CDF, ARCSIN_CDF_INV, ARCSIN_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2006
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
        Console.WriteLine("ARCSIN_CDF_TEST");
        Console.WriteLine("  ARCSIN_CDF evaluates the Arcsin CDF;");
        Console.WriteLine("  ARCSIN_CDF_INV inverts the Arcsin CDF.");
        Console.WriteLine("  ARCSIN_PDF evaluates the Arcsin PDF;");

        a = 1.0E+00;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =             " + a + "");

        if (!Arcsin.arcsin_check(a))
        {
            Console.WriteLine("");
            Console.WriteLine("TEST006 - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            x = Arcsin.arcsin_sample(a, ref seed);
            pdf = Arcsin.arcsin_pdf(x, a);
            cdf = Arcsin.arcsin_cdf(x, a);
            x2 = Arcsin.arcsin_cdf_inv(cdf, a);

            Console.WriteLine("  "
                              + x.ToString().PadLeft(12) + "  "
                              + pdf.ToString().PadLeft(12) + "  "
                              + cdf.ToString().PadLeft(12) + "  "
                              + x2.ToString().PadLeft(12) + "");
        }

        return;
    }

    static void arcsin_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    ARCSIN_SAMPLE_TEST tests ARCSIN_MEAN, ARCSIN_SAMPLE, ARCSIN_VARIANCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
    {
        int SAMPLE_NUM = 1000;

        double a = 0;
        int i;
        int j;
        double mean;
        int seed = 123456789;
        double variance;
        double[] x = new double[SAMPLE_NUM];
        double xmax;
        double xmin;

        Console.WriteLine("");
        Console.WriteLine("ARCSIN_SAMPLE_TEST");
        Console.WriteLine("  ARCSIN_MEAN computes the Arcsin mean;");
        Console.WriteLine("  ARCSIN_SAMPLE samples the Arcsin distribution;");
        Console.WriteLine("  ARCSIN_VARIANCE computes the Arcsin variance.");

        for (i = 1; i <= 2; i++)
        {
            if (i == 1)
            {
                a = 1.0E+00;
            }
            else if (i == 2)
            {
                a = 16.0E+00;
            }

            Console.WriteLine("");
            Console.WriteLine("  PDF parameter A =             " + a + "");

            if (!Arcsin.arcsin_check(a))
            {
                Console.WriteLine("");
                Console.WriteLine("ARCSIN_SAMPLE_TEST - Fatal error!");
                Console.WriteLine("  The parameters are not legal.");
                return;
            }

            mean = Arcsin.arcsin_mean(a);
            variance = Arcsin.arcsin_variance(a);

            Console.WriteLine("  PDF mean =                    " + mean + "");
            Console.WriteLine("  PDF variance =                " + variance + "");

            for (j = 0; j < SAMPLE_NUM; j++)
            {
                x[j] = Arcsin.arcsin_sample(a, ref seed);
            }

            mean = typeMethods.r8vec_mean(j, x);
            variance = typeMethods.r8vec_variance(j, x);
            xmax = typeMethods.r8vec_max(j, x);
            xmin = typeMethods.r8vec_min(j, x);

            Console.WriteLine("");
            Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
            Console.WriteLine("  Sample mean =     " + mean + "");
            Console.WriteLine("  Sample variance = " + variance + "");
            Console.WriteLine("  Sample maximum =  " + xmax + "");
            Console.WriteLine("  Sample minimum =  " + xmin + "");
        }

    }
        
    }
}