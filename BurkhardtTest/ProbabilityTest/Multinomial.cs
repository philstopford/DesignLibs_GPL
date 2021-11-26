﻿using System;
using System.Globalization;
using Burkardt.Probability;
using Burkardt.Types;
using Burkardt.Uniform;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void multinomial_coef_test()

//****************************************************************************80
//
//  Purpose:
//
//    MULTINOMIAL_TEST tests MULTINOMIAL_COEF1, MULTINOMIAL_COEF2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 April 2016
//
//  Author:
//
//    John Burkardt
//
    {
        const int MAXFACTOR = 5;

        int[] factor = new int[MAXFACTOR];
        int i;
        int ncomb1;
        int ncomb2;

        Console.WriteLine("");
        Console.WriteLine("MULTINOMIAL_TEST");
        Console.WriteLine("  MULTINOMIAL_COEF1 computes multinomial");
        Console.WriteLine("  coefficients using the Gamma function;");
        Console.WriteLine("  MULTINOMIAL_COEF2 computes multinomial");
        Console.WriteLine("  coefficients directly.");

        Console.WriteLine("");
        Console.WriteLine("  Line 10 of the BINOMIAL table:");
        Console.WriteLine("");

        int n = 10;
        int nfactor = 2;

        for (i = 0; i <= n; i++)
        {
            factor[0] = i;
            factor[1] = n - i;

            ncomb1 = Burkardt.PolynomialNS.Multinomial.multinomial_coef1(nfactor, factor);

            ncomb2 = Burkardt.PolynomialNS.Multinomial.multinomial_coef2(nfactor, factor);

            Console.WriteLine("  "
                              + factor[0].ToString(CultureInfo.InvariantCulture).PadLeft(2) + "  "
                              + factor[1].ToString(CultureInfo.InvariantCulture).PadLeft(2) + "  "
                              + ncomb1.ToString(CultureInfo.InvariantCulture).PadLeft(5) + "  "
                              + ncomb2.ToString(CultureInfo.InvariantCulture).PadLeft(5) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  Level 5 of the TRINOMIAL coefficients:");

        n = 5;
        nfactor = 3;

        for (i = 0; i <= n; i++)
        {
            factor[0] = i;

            Console.WriteLine("");

            int j;
            for (j = 0; j <= n - factor[0]; j++)
            {
                factor[1] = j;
                factor[2] = n - factor[0] - factor[1];

                ncomb1 = Burkardt.PolynomialNS.Multinomial.multinomial_coef1(nfactor, factor);

                ncomb2 = Burkardt.PolynomialNS.Multinomial.multinomial_coef2(nfactor, factor);

                Console.WriteLine("  "
                                  + factor[0].ToString(CultureInfo.InvariantCulture).PadLeft(2) + "  "
                                  + factor[1].ToString(CultureInfo.InvariantCulture).PadLeft(2) + "  "
                                  + factor[2].ToString(CultureInfo.InvariantCulture).PadLeft(2) + "  "
                                  + ncomb1.ToString(CultureInfo.InvariantCulture).PadLeft(5) + "  "
                                  + ncomb2.ToString(CultureInfo.InvariantCulture).PadLeft(5) + "");
            }
        }

    }

    private static void multinomial_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    MULTINOMIAL_SAMPLE_TEST tests MULTINOMIAL_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 April 2016
//
//  Author:
//
//    John Burkardt
//
    {
        const int B = 3;
        const int SAMPLE_NUM = 1000;

        double[] c =  {
                0.125, 0.500, 0.375
            }
            ;
        int i;
        int j;
        int seed = 123456789;
        int[] x = new int [B * SAMPLE_NUM];

        Console.WriteLine("");
        Console.WriteLine("MULTINOMIAL_SAMPLE_TEST");
        Console.WriteLine("  MULTINOMIAL_MEAN computes the Multinomial mean;");
        Console.WriteLine("  MULTINOMIAL_SAMPLE samples the Multinomial distribution;");
        Console.WriteLine("  MULTINOMIAL_VARIANCE computes the Multinomial variance;");

        int a = 5;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + B + "");
        typeMethods.r8vec_print(B, c, "  PDF parameter C:");

        if (!Multinomial.multinomial_check(a, B, c))
        {
            Console.WriteLine("");
            Console.WriteLine("MULTINOMIAL_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        double[] mean = Multinomial.multinomial_mean(a, B, c);
        double[] variance = Multinomial.multinomial_variance(a, B, c);
        typeMethods.r8vec_print(B, mean, "  PDF mean:");
        typeMethods.r8vec_print(B, variance, "  PDF variance:");

        for (j = 0; j < SAMPLE_NUM; j++)
        {
            int[] y = Multinomial.multinomial_sample(a, B, c, ref seed);
            for (i = 0; i < B; i++)
            {
                x[i + j * B] = y[i];
            }
        }

        mean = typeMethods.i4row_mean(B, SAMPLE_NUM, x);
        variance = typeMethods.i4row_variance(B, SAMPLE_NUM, x);
        int[] xmax = typeMethods.i4row_max(B, SAMPLE_NUM, x);
        int[] xmin = typeMethods.i4row_min(B, SAMPLE_NUM, x);

        Console.WriteLine("");
        Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
        Console.WriteLine("");
        Console.WriteLine("  Component Mean, Variance, Min, Max:");
        Console.WriteLine("");

        for (i = 0; i < B; i++)
        {
            Console.WriteLine("  "
                              + (i + 1).ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + mean[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + variance[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + xmin[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + xmax[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

    }

    private static void multinomial_pdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    MULTINOMIAL_PDF_TEST tests MULTINOMIAL_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 April 2016
//
//  Author:
//
//    John Burkardt
//
    {
        const int B = 3;

        double[] c =  {
                0.1, 0.5, 0.4
            }
            ;
        int[] x =  {
                0, 2, 3
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("MULTINOMIAL_PDF_TEST");
        Console.WriteLine("  MULTINOMIAL_PDF evaluates the Multinomial PDF;");

        const int a = 5;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + B + "");
        typeMethods.r8vec_print(B, c, "  PDF parameter C:");

        if (!Multinomial.multinomial_check(a, B, c))
        {
            Console.WriteLine("");
            Console.WriteLine("MULTINOMIAL_PDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        typeMethods.i4vec_print(B, x, "  PDF argument X:");

        double pdf = Multinomial.multinomial_pdf(x, a, B, c);

        Console.WriteLine("");
        Console.WriteLine("  PDF value = " + pdf + "");

    }

    private static void multinoulli_pdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    MULTINOULLLI_PDF_TEST tests MULTINOULLI_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2018
//
//  Author:
//
//    John Burkardt
//
    {
        const int n = 5;
        int x;

        Console.WriteLine("");
        Console.WriteLine("MULTINOULLI_PDF_TEST");
        Console.WriteLine("  MULTINOULLI_PDF evaluates the Multinoulli PDF.");

        int seed = 123456789;
        double[] theta = UniformRNG.r8vec_uniform_01_new(n, ref seed);
        double theta_sum = 0.0;
        for (x = 0; x < n; x++)
        {
            theta_sum += theta[x];
        }

        for (x = 0; x < n; x++)
        {
            theta[x] /= theta_sum;
        }

        Console.WriteLine("");
        Console.WriteLine("   X     pdf(X)");
        Console.WriteLine("");
        for (x = -1; x <= n; x++)
        {
            double pdf = Multinoulli.multinoulli_pdf(x, n, theta);
            Console.WriteLine("  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

    }

}