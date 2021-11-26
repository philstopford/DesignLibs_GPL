using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void hypergeometric_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    HYPERGEOMETRIC_CDF_TEST tests HYPERGEOMETRIC_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 April 2016
//
//  Author:
//
//    John Burkardt
//
    {
        Console.WriteLine("");
        Console.WriteLine("HYPERGEOMETRIC_CDF_TEST");
        Console.WriteLine("  HYPERGEOMETRIC_CDF evaluates the Hypergeometric CDF.");
        Console.WriteLine("  HYPERGEOMETRIC_PDF evaluates the Hypergeometric PDF.");

        const int n = 10;
        const int m = 7;
        const int l = 100;

        Console.WriteLine("");
        Console.WriteLine("  Total number of balls L =         " + l + "");
        Console.WriteLine("  Number of white balls M =         " + m + "");
        Console.WriteLine("  Number of balls taken N =         " + n + "");

        if (!Hypergeometric.hypergeometric_check(n, m, l))
        {
            Console.WriteLine("");
            Console.WriteLine("HYPERGEOMETRIC_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        int x = 7;

        double pdf = Hypergeometric.hypergeometric_pdf(x, n, m, l);

        double cdf = Hypergeometric.hypergeometric_cdf(x, n, m, l);

        Console.WriteLine("  PDF argument X =                " + x + "");
        Console.WriteLine("  PDF value =                   = " + pdf + "");
        Console.WriteLine("  CDF value =                   = " + cdf + "");
    }

    private static void hypergeometric_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    HYPERGEOMETRIC_SAMPLE_TEST tests HYPERGEOMETRIC_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 April 2016
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
        Console.WriteLine("HYPERGEOMETRIC_SAMPLE_TEST");
        Console.WriteLine("  HYPERGEOMETRIC_MEAN computes the Hypergeometric mean;");
        Console.WriteLine("  HYPERGEOMETRIC_SAMPLE samples the Hypergeometric distribution;");
        Console.WriteLine("  HYPERGEOMETRIC_VARIANCE computes the Hypergeometric variance.");

        const int n = 10;
        const int m = 7;
        const int l = 100;

        Console.WriteLine("");
        Console.WriteLine("  Total number of balls L =         " + l + "");
        Console.WriteLine("  Number of white balls M =         " + m + "");
        Console.WriteLine("  Number of balls taken N =         " + n + "");

        if (!Hypergeometric.hypergeometric_check(n, m, l))
        {
            Console.WriteLine("");
            Console.WriteLine("HYPERGEOMETRIC_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        double mean = Hypergeometric.hypergeometric_mean(n, m, l);
        double variance = Hypergeometric.hypergeometric_variance(n, m, l);

        Console.WriteLine("  PDF mean =                    " + mean + "");
        Console.WriteLine("  PDF variance =                " + variance + "");

        for (j = 0; j < SAMPLE_NUM; j++)
        {
            x[j] = Hypergeometric.hypergeometric_sample(n, m, l, ref seed );
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