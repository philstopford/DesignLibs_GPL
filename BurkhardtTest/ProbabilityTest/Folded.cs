﻿using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace Burkardt.ProbabilityTest
{
    partial class Program
    {
        static void folded_normal_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    FOLDED_NORMAL_CDF_TEST tests FOLDED_NORMAL_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 April 2016
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
            double x;
            double x2;

            Console.WriteLine("");
            Console.WriteLine("FOLDED_NORMAL_CDF_TEST");
            Console.WriteLine("  FOLDED_NORMAL_CDF evaluates the Folded Normal CDF;");
            Console.WriteLine("  FOLDED_NORMAL_CDF_INV inverts the Folded Normal CDF.");
            Console.WriteLine("  FOLDED_NORMAL_PDF evaluates the Folded Normal PDF;");

            a = 2.0;
            b = 3.0;

            Console.WriteLine("");
            Console.WriteLine("  PDF parameter A =      " + a + "");
            Console.WriteLine("  PDF parameter B =      " + b + "");

            if (!Folded.folded_normal_check(a, b))
            {
                Console.WriteLine("");
                Console.WriteLine("FOLDED_NORMAL_CDF_TEST - Fatal error!");
                Console.WriteLine("  The parameters are not legal.");
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("       X            PDF           CDF            CDF_INV");
            Console.WriteLine("");

            for (i = 1; i <= 10; i++)
            {
                x = Folded.folded_normal_sample(a, b, ref seed);
                pdf = Folded.folded_normal_pdf(x, a, b);
                cdf = Folded.folded_normal_cdf(x, a, b);
                x2 = Folded.folded_normal_cdf_inv(cdf, a, b);

                Console.WriteLine("  "
                                  + x.ToString().PadLeft(12) + "  "
                                  + pdf.ToString().PadLeft(12) + "  "
                                  + cdf.ToString().PadLeft(12) + "  "
                                  + x2.ToString().PadLeft(12) + "");
            }

        }

        static void folded_normal_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    FOLDED_NORMAL_SAMPLE_TEST tests FOLDED_NORMAL_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 April 2016
//
//  Author:
//
//    John Burkardt
//
        {
            int SAMPLE_NUM = 1000;

            double a;
            double b;
            int i;
            double mean;
            int seed = 123456789;
            double variance;
            double[] x = new double [SAMPLE_NUM];
            double xmax;
            double xmin;

            Console.WriteLine("");
            Console.WriteLine("FOLDED_NORMAL_SAMPLE_TEST");
            Console.WriteLine("  FOLDED_NORMAL_MEAN computes the Folded Normal mean;");
            Console.WriteLine("  FOLDED_NORMAL_SAMPLE samples the Folded Normal distribution;");
            Console.WriteLine("  FOLDED_NORMAL_VARIANCE computes the Folded Normal variance;");

            a = 2.0;
            b = 3.0;

            Console.WriteLine("");
            Console.WriteLine("  PDF parameter A =      " + a + "");
            Console.WriteLine("  PDF parameter B =      " + b + "");

            if (!Folded.folded_normal_check(a, b))
            {
                Console.WriteLine("");
                Console.WriteLine("FOLDED_NORMAL_SAMPLE_TEST - Fatal error!");
                Console.WriteLine("  The parameters are not legal.");
                return;
            }

            mean = Folded.folded_normal_mean(a, b);
            variance = Folded.folded_normal_variance(a, b);

            Console.WriteLine("");
            Console.WriteLine("  PDF mean =     " + mean + "");
            Console.WriteLine("  PDF variance = " + variance + "");

            for (i = 0; i < SAMPLE_NUM; i++)
            {
                x[i] = Folded.folded_normal_sample(a, b, ref seed);
            }

            mean = typeMethods.r8vec_mean(SAMPLE_NUM, x);
            variance = typeMethods.r8vec_variance(SAMPLE_NUM, x);
            xmax = typeMethods.r8vec_max(SAMPLE_NUM, x);
            xmin = typeMethods.r8vec_min(SAMPLE_NUM, x);

            Console.WriteLine("");
            Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
            Console.WriteLine("  Sample mean =     " + mean + "");
            Console.WriteLine("  Sample variance = " + variance + "");
            Console.WriteLine("  Sample maximum =  " + xmax + "");
            Console.WriteLine("  Sample minimum =  " + xmin + "");

        }

    }
}