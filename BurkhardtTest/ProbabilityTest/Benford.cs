using System;
using Burkardt.Probability;

namespace ProbabilityTest;

internal partial class Program
{
    private static void benford_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    BENFORD_CDF_TEST tests BENFORD_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2016
//
//  Author:
//
//    John Burkardt
//
    {
        double cdf;
        double cdf2;
        int n;
        double pdf;

        Console.WriteLine("");
        Console.WriteLine("BENFORD_CDF_TEST");
        Console.WriteLine("  BENFORD_CDF evaluates the Benford CDF.");

        Console.WriteLine("");
        Console.WriteLine("       N          CDF(N)          CDF(N) by summing");
        Console.WriteLine("");

        cdf2 = 0.0;
        for (n = 1; n <= 9; n++)
        {
            cdf = Benford.benford_cdf(n);
            pdf = Benford.benford_pdf(n);
            cdf2 += pdf;
            Console.WriteLine("  " + n.ToString().PadLeft(6)
                                   + "  " + pdf.ToString().PadLeft(14)
                                   + "  " + cdf.ToString().PadLeft(14)
                                   + "  " + cdf2.ToString().PadLeft(14) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("       N          CDF(N)          CDF(N) by summing");
        Console.WriteLine("");

        cdf2 = 0.0;
        for (n = 10; n <= 99; n++)
        {
            cdf = Benford.benford_cdf(n);
            pdf = Benford.benford_pdf(n);
            cdf2 += pdf;
            Console.WriteLine("  " + n.ToString().PadLeft(6)
                                   + "  " + pdf.ToString().PadLeft(14)
                                   + "  " + cdf.ToString().PadLeft(14)
                                   + "  " + cdf2.ToString().PadLeft(14) + "");
        }
    }

    private static void benford_pdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    BENFORD_PDF_TEST tests BENFORD_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2016
//
//  Author:
//
//    John Burkardt
//
    {
        int n;
        double pdf;

        Console.WriteLine("");
        Console.WriteLine("BENFORD_PDF_TEST");
        Console.WriteLine("  BENFORD_PDF evaluates the Benford PDF.");

        Console.WriteLine("");
        Console.WriteLine("       N          PDF(N)");
        Console.WriteLine("");

        for (n = 1; n <= 9; n++)
        {
            pdf = Benford.benford_pdf(n);
            Console.WriteLine("  " + n.ToString().PadLeft(6)
                                   + "  " + pdf.ToString().PadLeft(14) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("       N          PDF(N)");
        Console.WriteLine("");

        for (n = 10; n <= 99; n++)
        {
            pdf = Benford.benford_pdf(n);
            Console.WriteLine("  " + n.ToString().PadLeft(6)
                                   + "  " + pdf.ToString().PadLeft(14) + "");
        }
    }

}