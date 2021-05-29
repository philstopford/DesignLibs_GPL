using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace Burkardt.ProbabilityTest
{
    partial class Program
    {
        static void poisson_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    POISSON_CDF_TEST tests POISSON_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 March 2016
//
//  Author:
//
//    John Burkardt
//
        {
            double a;
            double cdf;
            int i;
            double pdf;
            int seed = 123456789;
            int x;
            int x2;

            Console.WriteLine("");
            Console.WriteLine("POISSON_CDF_TEST");
            Console.WriteLine("  POISSON_CDF evaluates the Poisson CDF;");
            Console.WriteLine("  POISSON_CDF_INV inverts the Poisson CDF.");
            Console.WriteLine("  POISSON_PDF evaluates the Poisson PDF;");

            a = 10.0E+00;

            Console.WriteLine("");
            Console.WriteLine("  PDF parameter A =             " + a + "");

            if (!Poisson.poisson_check(a))
            {
                Console.WriteLine("");
                Console.WriteLine("POISSON_CDF_TEST - Fatal error!");
                Console.WriteLine("  The parameters are not legal.");
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("       X            PDF           CDF            CDF_INV");
            Console.WriteLine("");

            for (i = 1; i <= 10; i++)
            {
                x = Poisson.poisson_sample(a, ref seed);
                pdf = Poisson.poisson_pdf(x, a);
                cdf = Poisson.poisson_cdf(x, a);
                x2 = Poisson.poisson_cdf_inv(cdf, a);

                Console.WriteLine("  "
                                  + x.ToString().PadLeft(12) + "  "
                                  + pdf.ToString().PadLeft(12) + "  "
                                  + cdf.ToString().PadLeft(12) + "  "
                                  + x2.ToString().PadLeft(12) + "");
            }

            return;
        }

        static void poisson_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    POISSON_SAMPLE_TEST tests POISSON_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
        {
            int SAMPLE_NUM = 1000;

            double a;
            int i;
            double mean;
            int seed = 123456789;
            double variance;
            int[] x = new int[SAMPLE_NUM];
            int xmax;
            int xmin;

            Console.WriteLine("");
            Console.WriteLine("POISSON_SAMPLE_TEST");
            Console.WriteLine("  POISSON_SAMPLE samples the Poisson PDF.");
            Console.WriteLine("  POISSON_SAMPLE samples the Poisson PDF.");
            Console.WriteLine("  POISSON_SAMPLE samples the Poisson PDF.");

            a = 10.0;

            Console.WriteLine("");
            Console.WriteLine("  PDF parameter A =      " + a + "");

            if (!Poisson.poisson_check(a))
            {
                Console.WriteLine("");
                Console.WriteLine("POISSON_SAMPLE_TEST - Fatal error!");
                Console.WriteLine("  The parameters are not legal.");
                return;
            }

            mean = Poisson.poisson_mean(a);
            variance = Poisson.poisson_variance(a);

            Console.WriteLine("");
            Console.WriteLine("  PDF mean =     " + mean + "");
            Console.WriteLine("  PDF variance = " + variance + "");

            for (i = 0; i < SAMPLE_NUM; i++)
            {
                x[i] = Poisson.poisson_sample(a, ref seed);
            }

            mean = typeMethods.i4vec_mean(SAMPLE_NUM, x);
            variance = typeMethods.i4vec_variance(SAMPLE_NUM, x);
            xmax = typeMethods.i4vec_max(SAMPLE_NUM, x);
            xmin = typeMethods.i4vec_min(SAMPLE_NUM, x);

            Console.WriteLine("");
            Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
            Console.WriteLine("  Sample mean =     " + mean + "");
            Console.WriteLine("  Sample variance = " + variance + "");
            Console.WriteLine("  Sample maximum =  " + xmax + "");
            Console.WriteLine("  Sample minimum =  " + xmin + "");

        }

    }
}