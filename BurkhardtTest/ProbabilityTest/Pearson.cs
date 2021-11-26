using System;
using Burkardt.Probability;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void pearson_05_pdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    PEARSON_05_PDF_TEST tests PEARSON_05_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 April 2016
//
//  Author:
//
//    John Burkardt
//
    {
        Console.WriteLine("");
        Console.WriteLine("PEARSON_05_PDF");
        Console.WriteLine("  PEARSON_05_PDF evaluates the Pearson 05 PDF.");

        const double x = 5.0;

        const double a = 1.0;
        const double b = 2.0;
        const double c = 3.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A = " + a + "");
        Console.WriteLine("  PDF parameter B = " + b + "");
        Console.WriteLine("  PDF parameter C = " + c + "");

        if (!Pearson.pearson_05_check(a, b, c))
        {
            Console.WriteLine("");
            Console.WriteLine("PEARSON_05_PDF - Fatal error!");
            Console.WriteLine("  The parameters are not legal.");
            return;
        }

        double pdf = Pearson.pearson_05_pdf(x, a, b, c);

        Console.WriteLine("");
        Console.WriteLine("  PDF argument X =  " + x + "");
        Console.WriteLine("  PDF value =       " + pdf + "");

    }
        
}