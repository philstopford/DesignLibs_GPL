using System;
using System.Globalization;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal static partial class Program
{
    public static void triangular_cdf_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULAR_CDF_TEST tests TRIANGULAR_CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 March 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("TRIANGULAR_CDF_TEST");
        Console.WriteLine("  TRIANGULAR_CDF evaluates the Triangular CDF;");
        Console.WriteLine("  TRIANGULAR_CDF_INV inverts the Triangular CDF.");
        Console.WriteLine("  TRIANGULAR_PDF evaluates the Triangular PDF;");

        double a = 1.0;
        double b = 10.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!Triangular.triangular_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("TRIANGULAR_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            double x = Triangular.triangular_sample(a, b, ref seed);
            double pdf = Triangular.triangular_pdf(x, a, b);
            double cdf = Triangular.triangular_cdf(x, a, b);
            double x2 = Triangular.triangular_cdf_inv(cdf, a, b);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

    }

    private static void triangular_sample_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULAR_SAMPLE_TEST tests TRIANGULAR_SAMPLE.
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
        const int SAMPLE_NUM = 1000;

        int i;
        int seed = 123456789;
        double[] x = new double [SAMPLE_NUM];

        Console.WriteLine("");
        Console.WriteLine("TRIANGULAR_SAMPLE_TEST");
        Console.WriteLine("  TRIANGULAR_MEAN computes the Triangular mean;");
        Console.WriteLine("  TRIANGULAR_SAMPLE samples the Triangular distribution;");
        Console.WriteLine("  TRIANGULAR_VARIANCE computes the Triangular variance;");

        double a = 1.0;
        double b = 10.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        if (!Triangular.triangular_check(a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("TRIANGULAR_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        double mean = Triangular.triangular_mean(a, b);
        double variance = Triangular.triangular_variance(a, b);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =     " + mean + "");
        Console.WriteLine("  PDF variance = " + variance + "");

        for (i = 0; i < SAMPLE_NUM; i++)
        {
            x[i] = Triangular.triangular_sample(a, b, ref seed);
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