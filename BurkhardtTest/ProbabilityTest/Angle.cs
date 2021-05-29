using System;
using Burkardt.Probability;

namespace Burkardt.ProbabilityTest
{
    partial class Program
    {
    static void angle_cdf_test()

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
        double cdf;
        int n;
        double x;

        Console.WriteLine("");
        Console.WriteLine("ANGLE_CDF_TEST");
        Console.WriteLine("  ANGLE_CDF evaluates the Angle CDF;");

        n = 5;
        x = 0.50E+00;

        cdf = Angle.angle_cdf(x, n);

        Console.WriteLine("");
        Console.WriteLine("  Parameter N =     " + n + "");
        Console.WriteLine("  PDF argument X =   " + x + "");
        Console.WriteLine("  CDF value =       " + cdf + "");

        return;
    }

    static void angle_pdf_test()

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
        int n;
        double pdf;
        double x;

        Console.WriteLine("");
        Console.WriteLine("ANGLE_PDF_TEST");
        Console.WriteLine("  ANGLE_PDF evaluates the Angle PDF;");

        n = 5;
        x = 0.50E+00;

        pdf = Angle.angle_pdf(x, n);

        Console.WriteLine("");
        Console.WriteLine("  Parameter N =    " + n + "");
        Console.WriteLine("  PDF argument X =  " + x + "");
        Console.WriteLine("  PDF value =      " + pdf + "");

        return;
    }

    static void angle_mean_test()

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
        double mean;
        int n;

        Console.WriteLine("");
        Console.WriteLine("ANGLE_MEAN_TEST");
        Console.WriteLine("  ANGLIT_MEAN computes the Angle mean;");

        n = 5;
        mean = Angle.angle_mean(n);

        Console.WriteLine("");
        Console.WriteLine("  Parameter N = " + n + "");
        Console.WriteLine("  PDF mean =    " + mean + "");

        return;
    }

        
    }
}