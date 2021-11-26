using System;
using System.Globalization;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void empirical_discrete_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    EMPIRICAL_DISCRETE_CDF_TEST tests EMPIRICAL_DISCRETE_CDF.
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
        const int A = 6;

        double[] b =
            {
                1.0, 1.0, 3.0, 2.0, 1.0, 2.0
            }
            ;
        double[] c =
            {
                0.0, 1.0, 2.0, 4.5, 6.0, 10.0
            }
            ;
        int i;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("EMPIRICAL_DISCRETE_CDF_TEST");
        Console.WriteLine("  EMPIRICAL_DISCRETE_CDF evaluates the Empirical Discrete CDF;");
        Console.WriteLine("  EMPIRICAL_DISCRETE_CDF_INV inverts the Empirical Discrete CDF.");
        Console.WriteLine("  EMPIRICAL_DISCRETE_PDF evaluates the Empirical Discrete PDF;");

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A = " + A + "");
        typeMethods.r8vec_print(A, b, "  PDF parameter B = ");
        typeMethods.r8vec_print(A, c, "  PDF parameter C = ");

        if (!EmpiricalDiscrete.empirical_discrete_check(A, b, c))
        {
            Console.WriteLine("");
            Console.WriteLine("EMPIRICAL_DISCRETE_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            double x = EmpiricalDiscrete.empirical_discrete_sample(A, b, c, ref seed);
            double pdf = EmpiricalDiscrete.empirical_discrete_pdf(x, A, b, c);
            double cdf = EmpiricalDiscrete.empirical_discrete_cdf(x, A, b, c);
            double x2 = EmpiricalDiscrete.empirical_discrete_cdf_inv(cdf, A, b, c);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

    }

    private static void empirical_discrete_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    EMPIRICAL_DISCRETE_SAMPLE_TEST tests EMPIRICAL_DISCRETE_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2008
//
//  Author:
//
//    John Burkardt
//
    {
        const int A = 6;
        const int SAMPLE_NUM = 1000;

        double[] b =
            {
                1.0, 1.0, 3.0, 2.0, 1.0, 2.0
            }
            ;
        double[] c =
            {
                0.0, 1.0, 2.0, 4.5, 6.0, 10.0
            }
            ;
        int i;
        int seed = 123456789;
        double[] x = new double [SAMPLE_NUM];

        Console.WriteLine("");
        Console.WriteLine("EMPIRICAL_DISCRETE_SAMPLE_TEST");
        Console.WriteLine("  EMPIRICAL_DISCRETE_MEAN computes the Empirical Discrete mean;");
        Console.WriteLine("  EMPIRICAL_DISCRETE_SAMPLE samples the Empirical Discrete distribution;");
        Console.WriteLine("  EMPIRICAL_DISCRETE_VARIANCE computes the Empirical Discrete variance.");

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A = " + A + "");
        typeMethods.r8vec_print(A, b, "  PDF parameter B = ");
        typeMethods.r8vec_print(A, c, "  PDF parameter C = ");

        if (!EmpiricalDiscrete.empirical_discrete_check(A, b, c))
        {
            Console.WriteLine("");
            Console.WriteLine("EMPIRICAL_DISCRETE_SAMPLE_TEST- Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        double mean = EmpiricalDiscrete.empirical_discrete_mean(A, b, c);
        double variance = EmpiricalDiscrete.empirical_discrete_variance(A, b, c);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = EmpiricalDiscrete.empirical_discrete_sample(A, b, c, ref seed);
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

}