using System;
using System.Globalization;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void discrete_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    DISCRETE_CDF_TEST tests DISCRETE_CDF.
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
                1.0, 2.0, 6.0, 2.0, 4.0, 1.0
            }
            ;
        int i;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("DISCRETE_CDF_TEST");
        Console.WriteLine("  DISCRETE_CDF evaluates the Discrete CDF;");
        Console.WriteLine("  DISCRETE_CDF_INV inverts the Discrete CDF.");
        Console.WriteLine("  DISCRETE_PDF evaluates the Discrete PDF;");

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + A + "");
        typeMethods.r8vec_print(A, b, "  PDF parameter B:");

        if (!Discrete.discrete_check(A, b))
        {
            Console.WriteLine("");
            Console.WriteLine("DISCRETE_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }


        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            int x = Discrete.discrete_sample(A, b, ref seed);
            double pdf = Discrete.discrete_pdf(x, A, b);
            double cdf = Discrete.discrete_cdf(x, A, b);
            int x2 = Discrete.discrete_cdf_inv(cdf, A, b);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

    }

    private static void discrete_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    DISCRETE_SAMPLE_TEST tests DISCRETE_SAMPLE.
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
        const int SAMPLE_NUM = 1000;

        double[] b =
            {
                1.0, 2.0, 6.0, 2.0, 4.0, 1.0
            }
            ;
        int i;
        int seed = 123456789;
        int[] x = new int[SAMPLE_NUM];

        Console.WriteLine("");
        Console.WriteLine("DISCRETE_SAMPLE_TEST");
        Console.WriteLine("  DISCRETE_MEAN computes the Discrete mean;");
        Console.WriteLine("  DISCRETE_SAMPLE samples the Discrete distribution;");
        Console.WriteLine("  DISCRETE_VARIANCE computes the Discrete variance;");

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + A + "");
        typeMethods.r8vec_print(A, b, "  PDF parameter B:");

        if (!Discrete.discrete_check(A, b))
        {
            Console.WriteLine("");
            Console.WriteLine("DISCRETE_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        double mean = Discrete.discrete_mean(A, b);
        double variance = Discrete.discrete_variance(A, b);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Discrete.discrete_sample(A, b, ref seed);
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