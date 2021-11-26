using System;
using System.Globalization;
using Burkardt.Probability;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void nakagami_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    NAKAGAMI_CDF_TEST tests NAKAGAMI_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 August 2016
//
//  Author:
//
//    John Burkardt
//
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("NAKAGAMI_CDF_TEST");
        Console.WriteLine("  NAKAGAMI_CDF evaluates the Nakagami CDF;");
        Console.WriteLine("  NAKAGAMI_CDF_INV inverts the Nakagami CDF;");
        Console.WriteLine("  NAKAGAMI_PDF evaluates the Nakagami PDF;");

        const double a = 1.0;
        const double b = 2.0;
        const double c = 3.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");
        Console.WriteLine("  PDF parameter C =      " + c + "");

        if (!Nakagami.nakagami_check(a, b, c))
        {
            Console.WriteLine("");
            Console.WriteLine("NAKAGAMI_CDF_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF         CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            double x = a + b * Math.Sqrt(i / c / 10.0);
            double pdf = Nakagami.nakagami_pdf(x, a, b, c);
            double cdf = Nakagami.nakagami_cdf(x, a, b, c);
            double x2 = Nakagami.nakagami_cdf_inv(cdf, a, b, c);

            Console.WriteLine("  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    private static void nakagami_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    NAKAGAMI_SAMPLE_TEST tests NAKAGAMI_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 April 2016
//
//  Author:
//
//    John Burkardt
//
    {
        Console.WriteLine("");
        Console.WriteLine("NAKAGAMI_SAMPLE_TEST");
        Console.WriteLine("  NAKAGAMI_MEAN evaluates the Nakagami mean;");
        Console.WriteLine("  NAKAGAMI_VARIANCE evaluates the Nakagami variance;");

        const double a = 1.0;
        const double b = 2.0;
        const double c = 3.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");
        Console.WriteLine("  PDF parameter C =      " + c + "");

        if (!Nakagami.nakagami_check(a, b, c))
        {
            Console.WriteLine("");
            Console.WriteLine("NAKAGAMI_SAMPLE_TEST - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        double mean = Nakagami.nakagami_mean(a, b, c);
        double variance = Nakagami.nakagami_variance(a, b, c);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =      " + mean + "");
        Console.WriteLine("  PDF variance =  " + variance + "");
    }

}