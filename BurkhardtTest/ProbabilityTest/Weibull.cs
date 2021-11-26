using System;
using System.Globalization;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void weibull_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    WEIBULL_CDF_TEST tests WEIBULL_CDF.
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
        Console.WriteLine("WEIBULL_CDF_TEST");
        Console.WriteLine("  WEIBULL_CDF evaluates the Weibull CDF;");
        Console.WriteLine("  WEIBULL_CDF_INV inverts the Weibull CDF.");
        Console.WriteLine("  WEIBULL_PDF evaluates the Weibull PDF;");

        double a = 2.0;
        double b = 3.0;
        double c = 4.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");
        Console.WriteLine("  PDF parameter C =      " + c + "");

        if (!Weibull.weibull_check(a, b, c))
        {
            Console.WriteLine("");
            Console.WriteLine("WEIBULL_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            double x = Weibull.weibull_sample(a, b, c, ref seed);
            double pdf = Weibull.weibull_pdf(x, a, b, c);
            double cdf = Weibull.weibull_cdf(x, a, b, c);
            double x2 = Weibull.weibull_cdf_inv(cdf, a, b, c);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    private static void weibull_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    WEIBULL_SAMPLE_TEST tests WEIBULL_SAMPLE.
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
        const int SAMPLE_NUM = 1000;

        int i;
        int seed = 123456789;
        double[] x = new double [SAMPLE_NUM];

        Console.WriteLine("");
        Console.WriteLine("WEIBULL_SAMPLE_TEST");
        Console.WriteLine("  WEIBULL_MEAN computes the Weibull mean;");
        Console.WriteLine("  WEIBULL_SAMPLE samples the Weibull distribution;");
        Console.WriteLine("  WEIBULL_VARIANCE computes the Weibull variance.");

        double a = 2.0;
        double b = 3.0;
        double c = 4.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");
        Console.WriteLine("  PDF parameter C =      " + c + "");

        if (!Weibull.weibull_check(a, b, c))
        {
            Console.WriteLine("");
            Console.WriteLine("WEIBULL_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        double mean = Weibull.weibull_mean(a, b, c);
        double variance = Weibull.weibull_variance(a, b, c);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Weibull.weibull_sample(a, b, c, ref seed);
        }

        mean = typeMethods.r8vec_mean(SAMPLE_NUM, x);
        variance = typeMethods.r8vec_variance(SAMPLE_NUM, x);
        double xmax = typeMethods.r8vec_max(SAMPLE_NUM, x);
        double xmin = typeMethods.r8vec_min(SAMPLE_NUM, x);

        Console.WriteLine("");
        Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
        Console.WriteLine("  Sample mean =     " + mean + "");
        Console.WriteLine("  Sample variance = " + variance + "");
        Console.WriteLine("  Sample maximum =  " + xmax + "");
        Console.WriteLine("  Sample minimum =  " + xmin + "");
    }

    private static void weibull_discrete_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    WEIBULL_DISCRETE_CDF_TEST tests WEIBULL_DISCRETE_CDF.
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
        int i;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("WEIBULL_DISCRETE_CDF_TEST");
        Console.WriteLine("  WEIBULL_DISCRETE_CDF evaluates the Weibull Discrete CDF;");
        Console.WriteLine("  WEIBULL_DISCRETE_CDF_INV inverts the Weibull Discrete CDF.");
        Console.WriteLine("  WEIBULL_DISCRETE_PDF evaluates the Weibull Discrete PDF;");

        const double a = 0.5;
        const double b = 1.5;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!Weibull.weibull_discrete_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("WEIBULL_DISCRETE_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            int x = Weibull.weibull_discrete_sample(a, b, ref seed);
            double pdf = Weibull.weibull_discrete_pdf(x, a, b);
            double cdf = Weibull.weibull_discrete_cdf(x, a, b);
            int x2 = Weibull.weibull_discrete_cdf_inv(cdf, a, b);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    private static void weibull_discrete_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    WEIBULL_DISCRETE_SAMPLE_TEST tests WEIBULL_DISCRETE_SAMPLE.
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

        int i;
        int seed = 123456789;
        double[] x = new double [SAMPLE_NUM];

        Console.WriteLine("");
        Console.WriteLine("WEIBULL_DISCRETE_SAMPLE_TEST");
        Console.WriteLine("  WEIBULL_DISCRETE_SAMPLE samples the Weibull Discrete distribution;");

        const double a = 0.5;
        const double b = 1.5;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!Weibull.weibull_discrete_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("WEIBULL_DISCRETE_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Weibull.weibull_discrete_sample(a, b, ref seed);
        }

        double mean = typeMethods.r8vec_mean(SAMPLE_NUM, x);
        double variance = typeMethods.r8vec_variance(SAMPLE_NUM, x);
        double xmax = typeMethods.r8vec_max(SAMPLE_NUM, x);
        double xmin = typeMethods.r8vec_min(SAMPLE_NUM, x);

        Console.WriteLine("");
        Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
        Console.WriteLine("  Sample mean =     " + mean + "");
        Console.WriteLine("  Sample variance = " + variance + "");
        Console.WriteLine("  Sample maximum =  " + xmax + "");
        Console.WriteLine("  Sample minimum =  " + xmin + "");

    }

}