using System;
using System.Globalization;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void zipf_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    ZIPF_CDF_TEST tests ZIPF_CDF, ZIPF_CDF_INV, ZIPF_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 March 2016
//
//  Author:
//
//    John Burkardt
//
    {
        int x;

        Console.WriteLine("");
        Console.WriteLine("ZIPF_CDF_TEST");
        Console.WriteLine("  ZIPF_CDF evaluates the Zipf CDF;");
        Console.WriteLine("  ZIPF_CDF_INV inverts the Zipf CDF;");
        Console.WriteLine("  ZIPF_PDF evaluates the Zipf PDF;");

        double a = 2.0E+00;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =             " + a + "");

        if (!Zipf.zipf_check(a))
        {
            Console.WriteLine("");
            Console.WriteLine("ZIPF_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF    CDF_INV()");
        Console.WriteLine("");

        for (x = 1; x <= 20; x++)
        {
            double pdf = Zipf.zipf_pdf(x, a);
            double cdf = Zipf.zipf_cdf(x, a);
            int x2 = Zipf.zipf_cdf_inv(a, cdf);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    private static void zipf_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    ZIPF_SAMPLE_TEST tests ZIPF_MEAN, ZIPF_SAMPLE, ZIPF_VARIANCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 March 2016
//
//  Author:
//
//    John Burkardt
//
    {
        const int SAMPLE_NUM = 1000;

        int j;
        int seed = 123456789;
        int[] x = new int[SAMPLE_NUM];

        Console.WriteLine("");
        Console.WriteLine("ZIPF_SAMPLE_TEST");
        Console.WriteLine("  ZIPF_MEAN computes the Zipf mean;");
        Console.WriteLine("  ZIPF_SAMPLE samples the Zipf distribution;");
        Console.WriteLine("  ZIPF_VARIANCE computes the Zipf variance.");

        double a = 4.0E+00;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =             " + a + "");

        if (!Zipf.zipf_check(a))
        {
            Console.WriteLine("");
            Console.WriteLine("ZIPF_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        double mean = Zipf.zipf_mean(a);
        double variance = Zipf.zipf_variance(a);

        Console.WriteLine("  PDF mean =                    " + mean + "");
        Console.WriteLine("  PDF variance =                " + variance + "");

        for (j = 0; j < SAMPLE_NUM; j++)
        {
            x[j] = Zipf.zipf_sample(a, ref seed);
        }

        mean = typeMethods.i4vec_mean(SAMPLE_NUM, x);
        variance = typeMethods.i4vec_variance(SAMPLE_NUM, x);
        int xmax = typeMethods.i4vec_max(SAMPLE_NUM, x);
        int xmin = typeMethods.i4vec_min(SAMPLE_NUM, x);

        Console.WriteLine("");
        Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
        Console.WriteLine("  Sample mean =     " + mean + "");
        Console.WriteLine("  Sample variance = " + variance + "");
        Console.WriteLine("  Sample maximum =  " + xmax + "");
        Console.WriteLine("  Sample minimum =  " + xmin + "");
    }        
}