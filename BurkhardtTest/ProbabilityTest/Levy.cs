using System;
using Burkardt.Probability;

namespace ProbabilityTest;

internal partial class Program
{
    private static void levy_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    LEVY_CDF_TEST tests LEVY_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 April 2016
//
//  Author:
//
//    John Burkardt
//
    {
        double a;
        double b;
        double cdf;
        int i;
        double pdf;
        int seed = 123456789;
        double x;
        double x2;

        Console.WriteLine("");
        Console.WriteLine("LEVY_CDF_TEST");
        Console.WriteLine("  LEVY_CDF evaluates the Levy CDF;");
        Console.WriteLine("  LEVY_CDF_INV inverts the Levy CDF.");
        Console.WriteLine("  LEVY_PDF evaluates the Levy PDF;");

        a = 1.0;
        b = 2.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            x = Levy.levy_sample(a, b, ref seed);
            pdf = Levy.levy_pdf(x, a, b);
            cdf = Levy.levy_cdf(x, a, b);
            x2 = Levy.levy_cdf_inv(cdf, a, b);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }
        
}