using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal partial class Program
{
    private static void birthday_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    BIRTHDAY_CDF_TEST tests BIRTHDAY_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 March 2016
//
//  Author:
//
//    John Burkardt
//
    {
        double cdf;
        int n;
        int n2;
        double pdf;

        Console.WriteLine("");
        Console.WriteLine("BIRTHDAY_CDF_TEST");
        Console.WriteLine("  BIRTHDAY_CDF evaluates the CDF;");
        Console.WriteLine("  BIRTHDAY_CDF_INV inverts the CDF.");
        Console.WriteLine("  BIRTHDAY_PDF evaluates the PDF;");

        Console.WriteLine("");
        Console.WriteLine("       N            PDF           CDF            CDF_INV");
        Console.WriteLine("");

        for (n = 1; n <= 30; n++)
        {
            pdf = Birthday.birthday_pdf(n);

            cdf = Birthday.birthday_cdf(n);

            n2 = Birthday.birthday_cdf_inv(cdf);

            Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + cdf.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + n2.ToString(CultureInfo.InvariantCulture).PadLeft(8)+ "");
        }
    }

    private static void birthday_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    BIRTHDAY_SAMPLE_TEST tests BIRTHDAY_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 March 2016
//
//  Author:
//
//    John Burkardt
//
    {
        int i;
        double mean;
        int n;
        int nsample = 10000;
        double pdf;
        int seed;
        int[] x = new int [nsample];

        seed = 12345678;

        Console.WriteLine("");
        Console.WriteLine("BIRTHDAY_SAMPLE_TEST");
        Console.WriteLine("  BIRTHDAY_SAMPLE samples the Birthday distribution.");
        Console.WriteLine("");
        Console.WriteLine("   N            Mean           PDF");
        Console.WriteLine("");

        for (n = 10; n <= 40; n++)
        {
            for (i = 0; i < nsample; i++)
            {
                x[i] = Birthday.birthday_sample(n, ref seed);
            }

            mean = typeMethods.i4vec_mean(nsample, x);
            pdf = Birthday.birthday_pdf(n);
            Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + mean.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

}