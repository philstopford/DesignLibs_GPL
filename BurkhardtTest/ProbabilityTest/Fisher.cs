using System;
using System.Globalization;
using Burkardt.Probability;

namespace ProbabilityTest;

internal static partial class Program
{
    private static void fisher_pdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    FISHER_PDF_TEST tests FISHER_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 April 2016
//
//  Author:
//
//    John Burkardt
//
    {
        double kappa = 0;
        double[] mu = new double[3];
        const int n = 10;
        int test;
        const int test_num = 3;

        Console.WriteLine("");
        Console.WriteLine("FISHER_SAMPLE_TEST");
        Console.WriteLine("  FISHER_PDF evaluates the Fisher PDF.");

        for (test = 1; test <= test_num; test++)
        {
            switch (test)
            {
                case 1:
                    kappa = 0.0;
                    mu[0] = 1.0;
                    mu[1] = 0.0;
                    mu[2] = 0.0;
                    break;
                case 2:
                    kappa = 0.5;
                    mu[0] = 1.0;
                    mu[1] = 0.0;
                    mu[2] = 0.0;
                    break;
                case 3:
                    kappa = 10.0;
                    mu[0] = 1.0;
                    mu[1] = 0.0;
                    mu[2] = 0.0;
                    break;
            }

            Console.WriteLine("");
            Console.WriteLine("  PDF parameters:");
            Console.WriteLine("    Concentration parameter KAPPA = " + kappa + "");
            Console.WriteLine("    Direction MU(1:3) = "
                              + "  " + mu[0]
                              + "  " + mu[1]
                              + "  " + mu[2] + "");

            Console.WriteLine("");
            Console.WriteLine("      X                         PDF");
            Console.WriteLine("");

            int seed = 123456789;

            int j;
            for (j = 0; j < n; j++)
            {
                double[] x = Fisher.fisher_sample(kappa, mu, 1, ref seed);

                double pdf = Fisher.fisher_pdf(x, kappa, mu);

                Console.WriteLine("  " + x[0].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + x[1].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + x[2].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");


            }

        }
    }

}