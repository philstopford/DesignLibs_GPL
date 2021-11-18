using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal partial class Program
{
    private static void chi_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    CHI_CDF_TEST tests CHI_CDF, CHI_CDF_INV, CHI_PDF;
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
        double c;
        double cdf;
        int i;
        double pdf;
        int seed = 123456789;
        double x;
        double x2;

        Console.WriteLine("");
        Console.WriteLine("CHI_CDF_TEST");
        Console.WriteLine("  CHI_CDF evaluates the Chi CDF;");
        Console.WriteLine("  CHI_CDF_INV inverts the Chi CDF.");
        Console.WriteLine("  CHI_PDF evaluates the Chi PDF;");

        a = 1.0;
        b = 2.0;
        c = 3.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");
        Console.WriteLine("  PDF parameter C =      " + c + "");

        if (!Chi.chi_check(a, b, c))
        {
            Console.WriteLine("");
            Console.WriteLine("CHI_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            x = Chi.chi_sample(a, b, c, ref seed);
            pdf = Chi.chi_pdf(x, a, b, c);
            cdf = Chi.chi_cdf(x, a, b, c);
            x2 = Chi.chi_cdf_inv(cdf, a, b, c);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    private static void chi_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    CHI_SAMPLE_TEST tests CHI_MEAN, CHI_SAMPLE, CHI_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
    {
        int SAMPLE_NUM = 1000;

        double a;
        double b;
        double c;
        int i;
        double mean;
        int seed = 123456789;
        double variance;
        double[] x = new double[SAMPLE_NUM];
        double xmax;
        double xmin;

        Console.WriteLine("");
        Console.WriteLine("CHI_SAMPLE_TEST");
        Console.WriteLine("  CHI_MEAN computes the Chi mean;");
        Console.WriteLine("  CHI_SAMPLE samples the Chi distribution;");
        Console.WriteLine("  CHI_VARIANCE computes the Chi variance;");

        a = 1.0;
        b = 2.0;
        c = 3.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");
        Console.WriteLine("  PDF parameter C =      " + c + "");

        if (!Chi.chi_check(a, b, c))
        {
            Console.WriteLine("");
            Console.WriteLine("CHI_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        mean = Chi.chi_mean(a, b, c);
        variance = Chi.chi_variance(a, b, c);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Chi.chi_sample(a, b, c, ref seed);
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

    private static void chi_square_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    CHI_SQUARE_CDF_TEST tests CHI_SQUARE_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
    {
        double a;
        double cdf;
        int i;
        double pdf;
        int seed = 123456789;
        double x;
        double x2;

        Console.WriteLine("");
        Console.WriteLine("CHI_SQUARE_CDF_TEST");
        Console.WriteLine("  CHI_SQUARE_CDF evaluates the Chi Square CDF;");
        Console.WriteLine("  CHI_SQUARE_CDF_INV inverts the Chi Square CDF.");
        Console.WriteLine("  CHI_SQUARE_PDF evaluates the Chi Square PDF;");

        a = 4.0E+00;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =             " + a + "");

        if (!Chi.chi_square_check(a))
        {
            Console.WriteLine("");
            Console.WriteLine("CHI_SQUARE_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            x = Chi.chi_square_sample(a, ref seed);
            pdf = Chi.chi_square_pdf(x, a);
            cdf = Chi.chi_square_cdf(x, a);
            x2 = Chi.chi_square_cdf_inv(cdf, a);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    private static void chi_square_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    CHI_SQUARE_SAMPLE_TEST tests CHI_SQUARE_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
    {
        int SAMPLE_NUM = 1000;

        double a;
        int j;
        double mean;
        int seed = 123456789;
        double variance;
        double[] x = new double[SAMPLE_NUM];
        double xmax;
        double xmin;

        Console.WriteLine("");
        Console.WriteLine("CHI_SQUARE_SAMPLE_TEST");
        Console.WriteLine("  CHI_SQUARE_MEAN computes the Chi Square mean;");
        Console.WriteLine("  CHI_SQUARE_SAMPLE samples the Chi Square distribution;");
        Console.WriteLine("  CHI_SQUARE_VARIANCE computes the Chi Square variance.");

        a = 10.0E+00;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =             " + a + "");

        if (!Chi.chi_square_check(a))
        {
            Console.WriteLine("");
            Console.WriteLine("CHI_SQUARE_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        mean = Chi.chi_square_mean(a);
        variance = Chi.chi_square_variance(a);

        Console.WriteLine("  PDF mean =                    " + mean + "");
        Console.WriteLine("  PDF variance =                " + variance + "");

        for (j = 0; j < SAMPLE_NUM; j++)
        {
            x[j] = Chi.chi_square_sample(a, ref seed);
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

    private static void chi_square_noncentral_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    CHI_SQUARE_NONCENTRAL_SAMPLE_TEST tests CHI_SQUARE_NONCENTRAL_SAMPLE.
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
        int SAMPLE_NUM = 1000;

        double a;
        double b;
        int i;
        double mean;
        int seed;
        int seed_init = 123456789;
        double variance;
        double[] x = new double [SAMPLE_NUM];
        double xmax;
        double xmin;

        Console.WriteLine("");
        Console.WriteLine("CHI_SQUARE_NONCENTRAL_SAMPLE_TEST");
        Console.WriteLine("  CHI_SQUARE_NONCENTRAL_MEAN computes the Chi Square Noncentral mean;");
        Console.WriteLine("  CHI_SQUARE_NONCENTRAL_SAMPLE samples the Chi Square Noncentral distribution;");
        Console.WriteLine("  CHI_SQUARE_NONCENTRAL_VARIANCE computes the Chi Square Noncentral variance;");

        a = 3.0;
        b = 2.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!Chi.chi_square_noncentral_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("CHI_SQUARE_NONCENTRAL_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        mean = Chi.chi_square_noncentral_mean(a, b);
        variance = Chi.chi_square_noncentral_variance(a, b);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        seed = seed_init;

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Chi.chi_square_noncentral_sample(a, b, ref seed);
        }

        mean = typeMethods.r8vec_mean(SAMPLE_NUM, x);
        variance = typeMethods.r8vec_variance(SAMPLE_NUM, x);
        xmax = typeMethods.r8vec_max(SAMPLE_NUM, x);
        xmin = typeMethods.r8vec_min(SAMPLE_NUM, x);

        Console.WriteLine("");
        Console.WriteLine("  Initial seed =     " + seed_init + "");
        Console.WriteLine("  Final seed =       " + seed + "");

        Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
        Console.WriteLine("  Sample mean =     " + mean + "");
        Console.WriteLine("  Sample variance = " + variance + "");
        Console.WriteLine("  Sample maximum =  " + xmax + "");
        Console.WriteLine("  Sample minimum =  " + xmin + "");

    }

}