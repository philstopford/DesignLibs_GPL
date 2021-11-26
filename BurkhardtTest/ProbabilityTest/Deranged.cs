using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void deranged_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    DERANGED_CDF_TEST tests DERANGED_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
    {
        int a;
        double cdf;
        double pdf;
        int x;
        int x2;

        Console.WriteLine("");
        Console.WriteLine("DERANGED_CDF_TEST");
        Console.WriteLine("  DERANGED_CDF evaluates the Deranged CDF;");
        Console.WriteLine("  DERANGED_CDF_INV inverts the Deranged CDF.");
        Console.WriteLine("  DERANGED_PDF evaluates the Deranged PDF;");

        a = 7;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =             " + a + "");

        if (!Deranged.deranged_check(a))
        {
            Console.WriteLine("");
            Console.WriteLine("DERANGED_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (x = 0; x <= a; x++)
        {
            pdf = Deranged.deranged_pdf(x, a);
            cdf = Deranged.deranged_cdf(x, a);
            x2 = Deranged.deranged_cdf_inv(cdf, a);

            Console.WriteLine("  "
                              + x.ToString().PadLeft(12) + "  "
                              + pdf.ToString().PadLeft(12) + "  "
                              + cdf.ToString().PadLeft(12) + "  "
                              + x2.ToString().PadLeft(12) + "");
        }
    }

    private static void deranged_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    DERANGED_SAMPLE_TEST tests DERANGED_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
    {
        int SAMPLE_NUM = 1000;

        int a;
        int j;
        double mean;
        int seed = 123456789;
        double variance;
        int[] x = new int[SAMPLE_NUM];
        int xmax;
        int xmin;

        Console.WriteLine("");
        Console.WriteLine("DERANGED_SAMPLE_TEST");
        Console.WriteLine("  DERANGED_MEAN computes the Deranged mean;");
        Console.WriteLine("  DERANGED_SAMPLE samples the Deranged distribution;");
        Console.WriteLine("  DERANGED_VARIANCE computes the Deranged variance.");

        a = 7;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =             " + a + "");

        if (!Deranged.deranged_check(a))
        {
            Console.WriteLine("");
            Console.WriteLine("DERANGED_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        mean = Deranged.deranged_mean(a);
        variance = Deranged.deranged_variance(a);

        Console.WriteLine("  PDF mean =                    " + mean + "");
        Console.WriteLine("  PDF variance =                " + variance + "");

        for (j = 0; j < SAMPLE_NUM; j++)
        {
            x[j] = Deranged.deranged_sample(a, ref seed);
        }

        mean = typeMethods.i4vec_mean(SAMPLE_NUM, x);
        variance = typeMethods.i4vec_variance(SAMPLE_NUM, x);
        xmax = typeMethods.i4vec_max(SAMPLE_NUM, x);
        xmin = typeMethods.i4vec_min(SAMPLE_NUM, x);

        Console.WriteLine("");
        Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
        Console.WriteLine("  Sample mean =     " + mean + "");
        Console.WriteLine("  Sample variance = " + variance + "");
        Console.WriteLine("  Sample maximum =  " + xmax + "");
        Console.WriteLine("  Sample minimum =  " + xmin + "");
    }

}