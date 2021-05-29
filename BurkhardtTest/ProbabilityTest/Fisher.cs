using System;
using Burkardt.Probability;

namespace Burkardt.ProbabilityTest
{
    partial class Program
    {
        static void fisher_pdf_test()

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
            int j;
            double kappa = 0;
            double[] mu = new double[3];
            int n = 10;
            double pdf;
            int seed;
            int test;
            int test_num = 3;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("FISHER_SAMPLE_TEST");
            Console.WriteLine("  FISHER_PDF evaluates the Fisher PDF.");

            for (test = 1; test <= test_num; test++)
            {
                if (test == 1)
                {
                    kappa = 0.0;
                    mu[0] = 1.0;
                    mu[1] = 0.0;
                    mu[2] = 0.0;
                }
                else if (test == 2)
                {
                    kappa = 0.5;
                    mu[0] = 1.0;
                    mu[1] = 0.0;
                    mu[2] = 0.0;
                }
                else if (test == 3)
                {
                    kappa = 10.0;
                    mu[0] = 1.0;
                    mu[1] = 0.0;
                    mu[2] = 0.0;
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

                seed = 123456789;

                for (j = 0; j < n; j++)
                {
                    x = Fisher.fisher_sample(kappa, mu, 1, ref seed);

                    pdf = Fisher.fisher_pdf(x, kappa, mu);

                    Console.WriteLine("  " + x[0].ToString().PadLeft(10)
                                      + "  " + x[1].ToString().PadLeft(10)
                                      + "  " + x[2].ToString().PadLeft(10)
                                      + "  " + pdf.ToString().PadLeft(14) + "");


                }

            }

            return;
        }

    }
}