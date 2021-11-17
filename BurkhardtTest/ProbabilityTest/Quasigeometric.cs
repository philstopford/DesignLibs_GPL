using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal partial class Program
{
    private static void quasigeometric_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    QUASIGEOMETRIC_CDF_TEST tests QUASIGEOMETRIC_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 April 2016
//
//  Author:
//
//    John Burkardt
//
    {
        double a;
        double b;
        double cdf;
        int i;
        double pdf;
        int seed = 123456789;
        int x;
        int x2;

        Console.WriteLine("");
        Console.WriteLine("QUASIGEOMETRIC_CDF_TEST");
        Console.WriteLine("  QUASIGEOMETRIC_CDF evaluates the Quasigeometric CDF;");
        Console.WriteLine("  QUASIGEOMETRIC_CDF_INV inverts the Quasigeometric CDF.");
        Console.WriteLine("  QUASIGEOMETRIC_PDF evaluates the Quasigeometric PDF;");

        a = 0.4825;
        b = 0.5893;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A = " + a + "");
        Console.WriteLine("  PDF parameter B = " + b + "");

        if (!Quasigeometric.quasigeometric_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("QUASIGEOMETRIC_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            x = Quasigeometric.quasigeometric_sample(a, b, ref seed);
            pdf = Quasigeometric.quasigeometric_pdf(x, a, b);
            cdf = Quasigeometric.quasigeometric_cdf(x, a, b);
            x2 = Quasigeometric.quasigeometric_cdf_inv(cdf, a, b);

            Console.WriteLine("  "
                              + x.ToString().PadLeft(12) + "  "
                              + pdf.ToString().PadLeft(12) + "  "
                              + cdf.ToString().PadLeft(12) + "  "
                              + x2.ToString().PadLeft(12) + "");
        }
    }

    private static void quasigeometric_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    QUASIGEOMETRIC_SAMPLE_TEST tests QUASIGEOMETRIC_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 April 2016
//
//  Author:
//
//    John Burkardt
//
    {
        double a;
        double b;
        int j;
        double mean;
        int sample_num = 1000;
        int seed = 123456789;
        double variance;
        int[] x;
        int xmax;
        int xmin;

        Console.WriteLine("");
        Console.WriteLine("QUASIGEOMETRIC_SAMPLE_TEST");
        Console.WriteLine("  QUASIGEOMETRIC_MEAN computes the Quasigeometric mean;");
        Console.WriteLine("  QUASIGEOMETRIC_SAMPLE samples the Quasigeometric distribution;");
        Console.WriteLine("  QUASIGEOMETRIC_VARIANCE computes the Quasigeometric variance.");

        a = 0.4825;
        b = 0.5893;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A = " + a + "");
        Console.WriteLine("  PDF parameter B = " + b + "");

        if (!Quasigeometric.quasigeometric_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("QUASIGEOMETRIC_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        mean = Quasigeometric.quasigeometric_mean(a, b);
        variance = Quasigeometric.quasigeometric_variance(a, b);

        Console.WriteLine("  PDF mean =                    " + mean + "");
        Console.WriteLine("  PDF variance =                " + variance + "");

        x = new int[sample_num];

        for (j = 0; j < sample_num; j++)
        {
            x[j] = Quasigeometric.quasigeometric_sample(a, b, ref seed);
        }

        mean = typeMethods.i4vec_mean(sample_num, x);
        variance = typeMethods.i4vec_variance(sample_num, x);
        xmax = typeMethods.i4vec_max(sample_num, x);
        xmin = typeMethods.i4vec_min(sample_num, x);

        Console.WriteLine("");
        Console.WriteLine("  Sample size =     " + sample_num + "");
        Console.WriteLine("  Sample mean =     " + mean + "");
        Console.WriteLine("  Sample variance = " + variance + "");
        Console.WriteLine("  Sample maximum =  " + xmax + "");
        Console.WriteLine("  Sample minimum =  " + xmin + "");
    }

}