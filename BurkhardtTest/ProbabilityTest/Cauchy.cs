using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest
{
    partial class Program
    {
        static void cauchy_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    CAUCHY_CDF_TEST tests CAUCHY_CDF, CAUCHY_CDF_INV, CAUCHY_PDF;
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
            double a;
            double b;
            double cdf;
            int i;
            double pdf;
            int seed = 123456789;
            double x;
            double x2;

            Console.WriteLine("");
            Console.WriteLine("CAUCHY_CDF_TEST");
            Console.WriteLine("  CAUCHY_CDF evaluates the Cauchy CDF;");
            Console.WriteLine("  CAUCHY_CDF_INV inverts the Cauchy CDF.");
            Console.WriteLine("  CAUCHY_PDF evaluates the Cauchy PDF;");

            a = 2.0;
            b = 3.0;

            Console.WriteLine("");
            Console.WriteLine("  PDF parameter A =      " + a + "");
            Console.WriteLine("  PDF parameter B =      " + b + "");

            if (!Cauchy.cauchy_check(a, b))
            {
                Console.WriteLine("");
                Console.WriteLine("CAUCHY_CDF_TEST - Fatal error!");
                Console.WriteLine("  The parameters are not legal.");
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("       X            PDF           CDF            CDF_INV");
            Console.WriteLine("");

            for (i = 1; i <= 10; i++)
            {
                x = Cauchy.cauchy_sample(a, b, ref seed);
                pdf = Cauchy.cauchy_pdf(x, a, b);
                cdf = Cauchy.cauchy_cdf(x, a, b);
                x2 = Cauchy.cauchy_cdf_inv(cdf, a, b);

                Console.WriteLine("  "
                                  + x.ToString().PadLeft(12) + "  "
                                  + pdf.ToString().PadLeft(12) + "  "
                                  + cdf.ToString().PadLeft(12) + "  "
                                  + x2.ToString().PadLeft(12) + "");
            }

            return;
        }

        static void cauchy_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    CAUCHY_SAMPLE_TEST tests CAUCHY_MEAN, CAUCHY_SAMPLE, CAUCHY_VARIANCE;
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
            int SAMPLE_NUM = 1000;

            double a;
            double b;
            int i;
            double mean;
            int seed = 123456789;
            double variance;
            double[] x = new double [SAMPLE_NUM];
            double xmax;
            double xmin;

            Console.WriteLine("");
            Console.WriteLine("CAUCHY_SAMPLE_TEST");
            Console.WriteLine("  CAUCHY_MEAN computes the Cauchy mean;");
            Console.WriteLine("  CAUCHY_SAMPLE samples the Cauchy distribution;");
            Console.WriteLine("  CAUCHY_VARIANCE computes the Cauchy variance;");

            a = 2.0;
            b = 3.0;

            Console.WriteLine("");
            Console.WriteLine("  PDF parameter A =      " + a + "");
            Console.WriteLine("  PDF parameter B =      " + b + "");

            if (!Cauchy.cauchy_check(a, b))
            {
                Console.WriteLine("");
                Console.WriteLine("CAUCHY_SAMPLE_TEST - Fatal error!");
                Console.WriteLine("  The parameters are not legal.");
                return;
            }

            mean = Cauchy.cauchy_mean(a, b);
            variance = Cauchy.cauchy_variance(a, b);

            Console.WriteLine("");
            Console.WriteLine("  PDF mean =     " + mean + "");
            Console.WriteLine("  PDF variance = " + variance + "");

            for (i = 0; i < SAMPLE_NUM; i++)
            {
                x[i] = Cauchy.cauchy_sample(a, b, ref seed);
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