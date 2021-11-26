using System;
using System.Globalization;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void geometric_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    GEOMETRIC_CDF_TEST tests GEOMETRIC_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2007
//
//  Author:
//
//    John Burkardt
//
    {
        int i;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("GEOMETRIC_CDF_TEST");
        Console.WriteLine("  GEOMETRIC_CDF evaluates the Geometric CDF;");
        Console.WriteLine("  GEOMETRIC_CDF_INV inverts the Geometric CDF.");
        Console.WriteLine("  GEOMETRIC_PDF evaluates the Geometric PDF;");

        const double a = 0.25E+00;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =             " + a + "");

        if (!Geometric.geometric_check(a))
        {
            Console.WriteLine("");
            Console.WriteLine("GEOMETRIC_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            int x = Geometric.geometric_sample(a, ref seed);
            double pdf = Geometric.geometric_pdf(x, a);
            double cdf = Geometric.geometric_cdf(x, a);
            int x2 = Geometric.geometric_cdf_inv(cdf, a);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    private static void geometric_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    GEOMETRIC_SAMPLE_TEST tests GEOMETRIC_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2007
//
//  Author:
//
//    John Burkardt
//
    {
        const int SAMPLE_NUM = 1000;

        int j;
        int seed = 123456789;
        int[] x = new int[SAMPLE_NUM];

        Console.WriteLine("");
        Console.WriteLine("GEOMETRIC_SAMPLE_TEST");
        Console.WriteLine("  GEOMETRIC_MEAN computes the Geometric mean;");
        Console.WriteLine("  GEOMETRIC_SAMPLE samples the Geometric distribution;");
        Console.WriteLine("  GEOMETRIC_VARIANCE computes the Geometric variance.");

        const double a = 0.25E+00;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =             " + a + "");

        if (!Geometric.geometric_check(a))
        {
            Console.WriteLine("");
            Console.WriteLine("GEOMETRIC_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        double mean = Geometric.geometric_mean(a);
        double variance = Geometric.geometric_variance(a);

        Console.WriteLine("  PDF mean =                    " + mean + "");
        Console.WriteLine("  PDF variance =                " + variance + "");

        for (j = 0; j < SAMPLE_NUM; j++)
        {
            x[j] = Geometric.geometric_sample(a, ref seed);
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