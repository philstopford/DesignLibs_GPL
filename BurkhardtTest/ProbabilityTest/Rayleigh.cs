using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest
{
    partial class Program
    {
        static void rayleigh_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    RAYLEIGH_CDF_TEST tests RAYLEIGH_CDF.
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
            double a;
            double cdf;
            int i;
            double pdf;
            int seed = 123456789;
            double x;
            double x2;

            Console.WriteLine("");
            Console.WriteLine("RAYLEIGH_CDF_TEST");
            Console.WriteLine("  RAYLEIGH_CDF evaluates the Rayleigh CDF;");
            Console.WriteLine("  RAYLEIGH_CDF_INV inverts the Rayleigh CDF.");
            Console.WriteLine("  RAYLEIGH_PDF evaluates the Rayleigh PDF;");

            a = 2.0;

            Console.WriteLine("");
            Console.WriteLine("  PDF parameter A =             " + a + "");

            if (!Rayleigh.rayleigh_check(a))
            {
                Console.WriteLine("");
                Console.WriteLine("RAYLEIGH_CDF_TEST - Fatal error!");
                Console.WriteLine("  The parameters are not legal.");
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("       X            PDF           CDF            CDF_INV");
            Console.WriteLine("");

            for (i = 1; i <= 10; i++)
            {
                x = Rayleigh.rayleigh_sample(a, ref seed);
                pdf = Rayleigh.rayleigh_pdf(x, a);
                cdf = Rayleigh.rayleigh_cdf(x, a);
                x2 = Rayleigh.rayleigh_cdf_inv(cdf, a);

                Console.WriteLine("  "
                                  + x.ToString().PadLeft(12) + "  "
                                  + pdf.ToString().PadLeft(12) + "  "
                                  + cdf.ToString().PadLeft(12) + "  "
                                  + x2.ToString().PadLeft(12) + "");
            }

            return;
        }

        static void rayleigh_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    RAYLEIGH_SAMPLE_TEST tests RAYLEIGH_SAMPLE.
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
            int j;
            double mean;
            int seed = 123456789;
            double variance;
            double[] x = new double [SAMPLE_NUM];
            double xmax;
            double xmin;

            Console.WriteLine("");
            Console.WriteLine("RAYLEIGH_SAMPLE_TEST");
            Console.WriteLine("  RAYLEIGH_MEAN computes the Rayleigh mean;");
            Console.WriteLine("  RAYLEIGH_SAMPLE samples the Rayleigh distribution;");
            Console.WriteLine("  RAYLEIGH_VARIANCE computes the Rayleigh variance.");

            a = 2.0E+00;

            Console.WriteLine("");
            Console.WriteLine("  PDF parameter A =             " + a + "");

            if (!Rayleigh.rayleigh_check(a))
            {
                Console.WriteLine("");
                Console.WriteLine("RAYLEIGH_SAMPLE_TEST - Fatal error!");
                Console.WriteLine("  The parameters are not legal.");
                return;
            }

            mean = Rayleigh.rayleigh_mean(a);
            variance = Rayleigh.rayleigh_variance(a);

            Console.WriteLine("  PDF mean =                    " + mean + "");
            Console.WriteLine("  PDF variance =                " + variance + "");

            for (j = 0; j < SAMPLE_NUM; j++)
            {
                x[j] = Rayleigh.rayleigh_sample(a, ref seed);
            }

            mean = typeMethods.r8vec_mean(SAMPLE_NUM, x);
            variance = typeMethods.r8vec_variance(SAMPLE_NUM, x);
            xmax = typeMethods.r8vec_max(SAMPLE_NUM, x);
            xmin = typeMethods.r8vec_min(SAMPLE_NUM, x);

            Console.WriteLine("");
            Console.WriteLine("  Sample size =     " + SAMPLE_NUM + "");
            Console.WriteLine("  Sample mean =     " + mean + "");
            Console.WriteLine("  Sample variance = " + variance + "");
            Console.WriteLine("  Sample maximum =  " + xmax + "");
            Console.WriteLine("  Sample minimum =  " + xmin + "");

        }

    }
}