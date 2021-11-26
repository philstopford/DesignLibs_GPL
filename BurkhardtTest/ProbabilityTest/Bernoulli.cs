﻿using System;
using System.Globalization;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void bernoulli_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    BERNOULLI_CDF_TEST tests BERNOULLI_CDF, BERNOULLI_CDF_INV, BERNOULLI_PDF.
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
        int i;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("BERNOULLI_CDF_TEST");
        Console.WriteLine("  BERNOULLI_CDF evaluates the Bernoulli CDF;");
        Console.WriteLine("  BERNOULLI_CDF_INV inverts the Bernoulli CDF;");
        Console.WriteLine("  BERNOULLI_PDF evaluates the Bernoulli PDF.");

        double a = 0.75E+00;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =             " + a + "");

        if (!Bernoulli.bernoulli_check(a))
        {
            Console.WriteLine("");
            Console.WriteLine("BERNOULLI_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            int x = Bernoulli.bernoulli_sample(a, ref seed);
            double pdf = Bernoulli.bernoulli_pdf(x, a);
            double cdf = Bernoulli.bernoulli_cdf(x, a);
            int x2 = Bernoulli.bernoulli_cdf_inv(cdf, a);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    private static void bernoulli_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    BERNOULLI_SAMPLE_TEST tests BERNOULLI_SAMPLE.
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
        const int SAMPLE_NUM = 1000;

        int i;
        int seed = 123456789;
        int[] x = new int [SAMPLE_NUM];

        Console.WriteLine("");
        Console.WriteLine("BERNOULLI_SAMPLE_TEST");
        Console.WriteLine("  BERNOULLI_MEAN computes the Bernoulli mean;");
        Console.WriteLine("  BERNOULLI_SAMPLE samples the Bernoulli distribution;");
        Console.WriteLine("  BERNOULLI_VARIANCE computes the Bernoulli variance.");

        double a = 0.75E+00;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A = " + a + "");

        if (!Bernoulli.bernoulli_check(a))
        {
            Console.WriteLine("");
            Console.WriteLine("BERNOULLI_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        double mean = Bernoulli.bernoulli_mean(a);
        double variance = Bernoulli.bernoulli_variance(a);

        Console.WriteLine("  PDF mean = " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Bernoulli.bernoulli_sample(a, ref seed);
        }

        mean = typeMethods.i4vec_mean(i, x);
        variance = typeMethods.i4vec_variance(i, x);
        int xmax = typeMethods.i4vec_max(i, x);
        int xmin = typeMethods.i4vec_min(i, x);

        Console.WriteLine("");
        Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
        Console.WriteLine("  Sample mean =     " + mean + "");
        Console.WriteLine("  Sample variance = " + variance + "");
        Console.WriteLine("  Sample maximum =  " + xmax + "");
        Console.WriteLine("  Sample minimum =  " + xmin + "");
    }
}