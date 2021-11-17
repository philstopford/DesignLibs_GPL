using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal partial class Program
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
        double cdf;
        int l;
        int m;
        int n;
        double pdf;
        int x;

        Console.WriteLine("");
        Console.WriteLine("HYPERGEOMETRIC_CDF_TEST");
        Console.WriteLine("  HYPERGEOMETRIC_CDF evaluates the Hypergeometric CDF.");
        Console.WriteLine("  HYPERGEOMETRIC_PDF evaluates the Hypergeometric PDF.");

        n = 10;
        m = 7;
        l = 100;

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

        x = 7;

        pdf = Hypergeometric.hypergeometric_pdf(x, n, m, l);

        cdf = Hypergeometric.hypergeometric_cdf(x, n, m, l);

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
        int SAMPLE_NUM = 1000;

        int j;
        int l;
        int m;
        double mean;
        int n;
        int seed = 123456789;
        double variance;
        int[] x = new int[SAMPLE_NUM];
        int xmax;
        int xmin;

        Console.WriteLine("");
        Console.WriteLine("HYPERGEOMETRIC_SAMPLE_TEST");
        Console.WriteLine("  HYPERGEOMETRIC_MEAN computes the Hypergeometric mean;");
        Console.WriteLine("  HYPERGEOMETRIC_SAMPLE samples the Hypergeometric distribution;");
        Console.WriteLine("  HYPERGEOMETRIC_VARIANCE computes the Hypergeometric variance.");

        n = 10;
        m = 7;
        l = 100;

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

        mean = Hypergeometric.hypergeometric_mean(n, m, l);
        variance = Hypergeometric.hypergeometric_variance(n, m, l);

        Console.WriteLine("  PDF mean =                    " + mean + "");
        Console.WriteLine("  PDF variance =                " + variance + "");

        for (j = 0; j < SAMPLE_NUM; j++)
        {
            x[j] = Hypergeometric.hypergeometric_sample(n, m, l, ref seed );
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