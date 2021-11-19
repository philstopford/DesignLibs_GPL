using System;
using Burkardt.Probability;

namespace ProbabilityTest;

internal partial class Program
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
        double a;
        double b;
        double c;
        double cdf;
        int i;
        double pdf;
        double x;
        double x2;

        Console.WriteLine("");
        Console.WriteLine("NAKAGAMI_CDF_TEST");
        Console.WriteLine("  NAKAGAMI_CDF evaluates the Nakagami CDF;");
        Console.WriteLine("  NAKAGAMI_CDF_INV inverts the Nakagami CDF;");
        Console.WriteLine("  NAKAGAMI_PDF evaluates the Nakagami PDF;");

        a = 1.0;
        b = 2.0;
        c = 3.0;

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
            x = a + b * Math.Sqrt(i / c / 10.0);
            pdf = Nakagami.nakagami_pdf(x, a, b, c);
            cdf = Nakagami.nakagami_cdf(x, a, b, c);
            x2 = Nakagami.nakagami_cdf_inv(cdf, a, b, c);

            Console.WriteLine("  " + x.ToString().PadLeft(12)
                                   + "  " + pdf.ToString().PadLeft(12)
                                   + "  " + cdf.ToString().PadLeft(12)
                                   + "  " + x2.ToString().PadLeft(12) + "");
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
        double a;
        double b;
        double c;
        double mean;
        double variance;

        Console.WriteLine("");
        Console.WriteLine("NAKAGAMI_SAMPLE_TEST");
        Console.WriteLine("  NAKAGAMI_MEAN evaluates the Nakagami mean;");
        Console.WriteLine("  NAKAGAMI_VARIANCE evaluates the Nakagami variance;");

        a = 1.0;
        b = 2.0;
        c = 3.0;

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

        mean = Nakagami.nakagami_mean(a, b, c);
        variance = Nakagami.nakagami_variance(a, b, c);

        Console.WriteLine("");
        Console.WriteLine("  PDF mean =      " + mean + "");
        Console.WriteLine("  PDF variance =  " + variance + "");
    }

}