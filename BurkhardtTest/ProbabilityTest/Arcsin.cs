﻿using System;
using System.Globalization;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void arcsin_cdf_test()

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
        int i;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("ARCSIN_CDF_TEST");
        Console.WriteLine("  ARCSIN_CDF evaluates the Arcsin CDF;");
        Console.WriteLine("  ARCSIN_CDF_INV inverts the Arcsin CDF.");
        Console.WriteLine("  ARCSIN_PDF evaluates the Arcsin PDF;");

        double a = 1.0E+00;

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
            double x = Arcsin.arcsin_sample(a, ref seed);
            double pdf = Arcsin.arcsin_pdf(x, a);
            double cdf = Arcsin.arcsin_cdf(x, a);
            double x2 = Arcsin.arcsin_cdf_inv(cdf, a);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    private static void arcsin_sample_test()

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
        const int SAMPLE_NUM = 1000;

        double a = 0;
        int i;
        int seed = 123456789;
        double[] x = new double[SAMPLE_NUM];

        Console.WriteLine("");
        Console.WriteLine("ARCSIN_SAMPLE_TEST");
        Console.WriteLine("  ARCSIN_MEAN computes the Arcsin mean;");
        Console.WriteLine("  ARCSIN_SAMPLE samples the Arcsin distribution;");
        Console.WriteLine("  ARCSIN_VARIANCE computes the Arcsin variance.");

        for (i = 1; i <= 2; i++)
        {
            a = i switch
            {
                1 => 1.0E+00,
                2 => 16.0E+00,
                _ => a
            };

            Console.WriteLine("");
            Console.WriteLine("  PDF parameter A =             " + a + "");

            if (!Arcsin.arcsin_check(a))
            {
                Console.WriteLine("");
                Console.WriteLine("ARCSIN_SAMPLE_TEST - Fatal error!");
                Console.WriteLine("  The parameters are not legal.");
                return;
            }

            double mean = Arcsin.arcsin_mean(a);
            double variance = Arcsin.arcsin_variance(a);

            Console.WriteLine("  PDF mean =                    " + mean + "");
            Console.WriteLine("  PDF variance =                " + variance + "");

            int j;
            for (j = 0; j < SAMPLE_NUM; j++)
            {
                x[j] = Arcsin.arcsin_sample(a, ref seed);
            }

            mean = typeMethods.r8vec_mean(j, x);
            variance = typeMethods.r8vec_variance(j, x);
            double xmax = typeMethods.r8vec_max(j, x);
            double xmin = typeMethods.r8vec_min(j, x);

            Console.WriteLine("");
            Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
            Console.WriteLine("  Sample mean =     " + mean + "");
            Console.WriteLine("  Sample variance = " + variance + "");
            Console.WriteLine("  Sample maximum =  " + xmax + "");
            Console.WriteLine("  Sample minimum =  " + xmin + "");
        }

    }
        
}