using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void beta_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    BETA_CDF_TEST tests BETA_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 April 2013
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
        Console.WriteLine("BETA_CDF_TEST");
        Console.WriteLine("  BETA_CDF evaluates the Beta CDF;");
        Console.WriteLine("  BETA_CDF_INV inverts the Beta CDF.");
        Console.WriteLine("  BETA_PDF evaluates the Beta PDF;");

        a = 12.0;
        b = 12.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!Beta.beta_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("BETA_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("             A             B        X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            x = Beta.beta_sample(a, b, ref seed);
            pdf = Beta.beta_pdf(x, a, b);
            cdf = Beta.beta_cdf(x, a, b);
            x2 = Beta.beta_cdf_inv(cdf, a, b);

            Console.WriteLine("  " + a.ToString().PadLeft(12)
                                   + "  " + b.ToString().PadLeft(12)
                                   + "  " + x.ToString().PadLeft(12)
                                   + "  " + pdf.ToString().PadLeft(12)
                                   + "  " + cdf.ToString().PadLeft(12)
                                   + "  " + x2.ToString().PadLeft(12) + "");
        }

    }

    private static void beta_inc_test()

//****************************************************************************80
//
//  Purpose:
//
//    BETA_INC_TEST tests BETA_INC.
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
        double a = 0;
        double b = 0;
        double fx = 0;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("BETA_INC_TEST:");
        Console.WriteLine("  BETA_INC evaluates the normalized incomplete Beta");
        Console.WriteLine("  function BETA_INC(A,B,X).");
        Console.WriteLine("");
        Console.WriteLine("         A         B         X       Exact F       BETA_INC(A,B,X)");
        Console.WriteLine("");

        int n_data = 0;

        for (;;)
        {
            Beta.beta_inc_values(ref n_data, ref a, ref b, ref x, ref fx);

            if (n_data == 0)
            {
                break;
            }

            double fx2 = Beta.beta_inc(a, b, x);

            Console.WriteLine("  "
                              + a.ToString().PadLeft(8) + "  "
                              + b.ToString().PadLeft(8) + "  "
                              + x.ToString().PadLeft(8) + "  "
                              + fx.ToString().PadLeft(16) + "  "
                              + fx2.ToString().PadLeft(16) + "");
        }

    }

    private static void beta_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    BETA_SAMPLE_TEST tests BETA_SAMPLE.
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
        int SAMPLE_NUM = 1000;

        double a;
        double b;
        int i;
        double mean;
        int seed = 123456789;
        double variance;
        double[] x = new double[SAMPLE_NUM];
        double xmax;
        double xmin;

        Console.WriteLine("");
        Console.WriteLine("BETA_SAMPLE_TEST");
        Console.WriteLine("  BETA_MEAN computes the Beta mean;");
        Console.WriteLine("  BETA_SAMPLE samples the Beta distribution;");
        Console.WriteLine("  BETA_VARIANCE computes the Beta variance;");

        a = 2.0;
        b = 3.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!Beta.beta_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("BETA_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        mean = Beta.beta_mean(a, b);
        variance = Beta.beta_variance(a, b);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Beta.beta_sample(a, b, ref seed);
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

    private static void beta_binomial_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    BETA_BINOMIAL_CDF_TEST tests BETA_BINOMIAL_CDF.
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
        double a;
        double b;
        int c;
        double cdf;
        int i;
        double pdf;
        int seed = 123456789;
        int x;
        int x2;

        Console.WriteLine("");
        Console.WriteLine("BETA_BINOMIAL_CDF_TEST");
        Console.WriteLine("  BETA_BINOMIAL_CDF evaluates the Beta Binomial CDF;");
        Console.WriteLine("  BETA_BINOMIAL_CDF_INV inverts the Beta Binomial CDF.");
        Console.WriteLine("  BETA_BINOMIAL_PDF evaluates the Beta Binomial PDF;");

        a = 2.0;
        b = 3.0;
        c = 4;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");
        Console.WriteLine("  PDF parameter C =      " + c + "");

        if (!Beta.beta_binomial_check(a, b, c))
        {
            Console.WriteLine("");
            Console.WriteLine("BETA_BINOMIAL_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            x = Beta.beta_binomial_sample(a, b, c, ref seed);
            pdf = Beta.beta_binomial_pdf(x, a, b, c);
            cdf = Beta.beta_binomial_cdf(x, a, b, c);
            x2 = Beta.beta_binomial_cdf_inv(cdf, a, b, c);

            Console.WriteLine("  "
                              + x.ToString().PadLeft(12) + "  "
                              + pdf.ToString().PadLeft(12) + "  "
                              + cdf.ToString().PadLeft(12) + "  "
                              + x2.ToString().PadLeft(12) + "");
        }
    }

    private static void beta_binomial_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    BETA_BINOMIAL_SAMPLE_TEST tests BETA_BINOMIAL_SAMPLE.
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
        int SAMPLE_NUM = 1000;

        double a;
        double b;
        int c;
        int i;
        double mean;
        int seed = 123456789;
        double variance;
        int[] x = new int[SAMPLE_NUM];
        int xmax;
        int xmin;

        Console.WriteLine("");
        Console.WriteLine("BETA_BINOMIAL_SAMPLE_TEST");
        Console.WriteLine("  BETA_BINOMIAL_MEAN computes the Beta Binomial mean;");
        Console.WriteLine("  BETA_BINOMIAL_SAMPLE samples the Beta Binomial distribution;");
        Console.WriteLine("  BETA_BINOMIAL_VARIANCE computes the Beta Binomial variance;");

        a = 2.0;
        b = 3.0;
        c = 4;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");
        Console.WriteLine("  PDF parameter C =      " + c + "");

        if (!Beta.beta_binomial_check(a, b, c))
        {
            Console.WriteLine("");
            Console.WriteLine("BETA_BINOMIAL_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        mean = Beta.beta_binomial_mean(a, b, c);
        variance = Beta.beta_binomial_variance(a, b, c);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Beta.beta_binomial_sample(a, b, c, ref seed);
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