using System;
using Burkardt.Probability;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void angle_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    ANGLE_CDF_TEST tests ANGLE_CDF.
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
        Console.WriteLine("");
        Console.WriteLine("ANGLE_CDF_TEST");
        Console.WriteLine("  ANGLE_CDF evaluates the Angle CDF;");

        const int n = 5;
        const double x = 0.50E+00;

        double cdf = Angle.angle_cdf(x, n);

        Console.WriteLine("");
        Console.WriteLine("  Parameter N =     " + n + "");
        Console.WriteLine("  PDF argument X =   " + x + "");
        Console.WriteLine("  CDF value =       " + cdf + "");
    }

    private static void angle_pdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    ANGLE_PDF_TEST tests ANGLE_PDF.
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
        Console.WriteLine("");
        Console.WriteLine("ANGLE_PDF_TEST");
        Console.WriteLine("  ANGLE_PDF evaluates the Angle PDF;");

        const int n = 5;
        const double x = 0.50E+00;

        double pdf = Angle.angle_pdf(x, n);

        Console.WriteLine("");
        Console.WriteLine("  Parameter N =    " + n + "");
        Console.WriteLine("  PDF argument X =  " + x + "");
        Console.WriteLine("  PDF value =      " + pdf + "");
    }

    private static void angle_mean_test()

//****************************************************************************80
//
//  Purpose:
//
//    ANGLE_MEAN_TEST tests ANGLE_MEAN;
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
        Console.WriteLine("");
        Console.WriteLine("ANGLE_MEAN_TEST");
        Console.WriteLine("  ANGLIT_MEAN computes the Angle mean;");

        const int n = 5;
        double mean = Angle.angle_mean(n);

        Console.WriteLine("");
        Console.WriteLine("  Parameter N = " + n + "");
        Console.WriteLine("  PDF mean =    " + mean + "");
    }


}