using System;
using System.Globalization;
using Burkardt.Probability;

namespace ProbabilityTest;

internal static partial class Program
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
        int i;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("LEVY_CDF_TEST");
        Console.WriteLine("  LEVY_CDF evaluates the Levy CDF;");
        Console.WriteLine("  LEVY_CDF_INV inverts the Levy CDF.");
        Console.WriteLine("  LEVY_PDF evaluates the Levy PDF;");

        const double a = 1.0;
        const double b = 2.0;

        Console.WriteLine("");
        Console.WriteLine("  PDF parameter A =      " + a + "");
        Console.WriteLine("  PDF parameter B =      " + b + "");

        Console.WriteLine("");
        Console.WriteLine("       X            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (i = 1; i <= 10; i++)
        {
            double x = Levy.levy_sample(a, b, ref seed);
            double pdf = Levy.levy_pdf(x, a, b);
            double cdf = Levy.levy_cdf(x, a, b);
            double x2 = Levy.levy_cdf_inv(cdf, a, b);

            Console.WriteLine("  "
                              + x.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  "
                              + x2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }
        
}