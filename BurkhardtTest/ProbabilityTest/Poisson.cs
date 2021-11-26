using System;
using System.Globalization;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void poisson_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    POISSON_CDF_TEST tests POISSON_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 March 2016
//
//  Author:
//
//    John Burkardt
//
    {
        int i;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("POISSON_CDF_TEST");
        Console.WriteLine("  POISSON_CDF evaluates the Poisson CDF;");
        Console.WriteLine("  POISSON_CDF_INV inverts the Poisson CDF.");
        Console.WriteLine("  POISSON_PDF evaluates the Poisson PDF;");

        double a = 10.0E+00;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =             " + a + "");

        if (!Poisson.poisson_check(a))
        {
            Console.WriteLine("");
            Console.WriteLine("POISSON_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            int x = Poisson.poisson_sample(a, ref seed);
            double pdf = Poisson.poisson_pdf(x, a);
            double cdf = Poisson.poisson_cdf(x, a);
            int x2 = Poisson.poisson_cdf_inv(cdf, a);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    private static void poisson_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    POISSON_SAMPLE_TEST tests POISSON_SAMPLE.
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
        int[] x = new int[SAMPLE_NUM];

        Console.WriteLine("");
        Console.WriteLine("POISSON_SAMPLE_TEST");
        Console.WriteLine("  POISSON_SAMPLE samples the Poisson PDF.");
        Console.WriteLine("  POISSON_SAMPLE samples the Poisson PDF.");
        Console.WriteLine("  POISSON_SAMPLE samples the Poisson PDF.");

        double a = 10.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");

        if (!Poisson.poisson_check(a))
        {
            Console.WriteLine("");
            Console.WriteLine("POISSON_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        double mean = Poisson.poisson_mean(a);
        double variance = Poisson.poisson_variance(a);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Poisson.poisson_sample(a, ref seed);
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