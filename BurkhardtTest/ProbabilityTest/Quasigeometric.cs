using System;
using System.Globalization;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal static partial class Program
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
        int i;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("QUASIGEOMETRIC_CDF_TEST");
        Console.WriteLine("  QUASIGEOMETRIC_CDF evaluates the Quasigeometric CDF;");
        Console.WriteLine("  QUASIGEOMETRIC_CDF_INV inverts the Quasigeometric CDF.");
        Console.WriteLine("  QUASIGEOMETRIC_PDF evaluates the Quasigeometric PDF;");

        const double a = 0.4825;
        const double b = 0.5893;

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
            int x = Quasigeometric.quasigeometric_sample(a, b, ref seed);
            double pdf = Quasigeometric.quasigeometric_pdf(x, a, b);
            double cdf = Quasigeometric.quasigeometric_cdf(x, a, b);
            int x2 = Quasigeometric.quasigeometric_cdf_inv(cdf, a, b);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
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
        int j;
        const int sample_num = 1000;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("QUASIGEOMETRIC_SAMPLE_TEST");
        Console.WriteLine("  QUASIGEOMETRIC_MEAN computes the Quasigeometric mean;");
        Console.WriteLine("  QUASIGEOMETRIC_SAMPLE samples the Quasigeometric distribution;");
        Console.WriteLine("  QUASIGEOMETRIC_VARIANCE computes the Quasigeometric variance.");

        const double a = 0.4825;
        const double b = 0.5893;

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

        double mean = Quasigeometric.quasigeometric_mean(a, b);
        double variance = Quasigeometric.quasigeometric_variance(a, b);

        Console.WriteLine("  PDF mean =                    " + mean + "");
        Console.WriteLine("  PDF variance =                " + variance + "");

        int[] x = new int[sample_num];

        for (j = 0; j < sample_num; j++)
        {
            x[j] = Quasigeometric.quasigeometric_sample(a, b, ref seed);
        }

        mean = typeMethods.i4vec_mean(sample_num, x);
        variance = typeMethods.i4vec_variance(sample_num, x);
        int xmax = typeMethods.i4vec_max(sample_num, x);
        int xmin = typeMethods.i4vec_min(sample_num, x);

        Console.WriteLine("");
        Console.WriteLine("  Sample size =     " + sample_num + "");
        Console.WriteLine("  Sample mean =     " + mean + "");
        Console.WriteLine("  Sample variance = " + variance + "");
        Console.WriteLine("  Sample maximum =  " + xmax + "");
        Console.WriteLine("  Sample minimum =  " + xmin + "");
    }

}