using System;
using Burkardt.Types;

namespace ProbabilityTest
{
    partial class Program
    {
        static void half_normal_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    HALF_NORMAL_CDF_TEST tests HALF_NORMAL_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 April 2016
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
            Console.WriteLine("HALF_NORMAL_CDF_TEST");
            Console.WriteLine("  HALF_NORMAL_CDF evaluates the Half Normal CDF;");
            Console.WriteLine("  HALF_NORMAL_CDF_INV inverts the Half Normal CDF.");
            Console.WriteLine("  HALF_NORMAL_PDF evaluates the Half Normal PDF;");

            a = 0.0;
            b = 2.0;

            Console.WriteLine("");
            Console.WriteLine("  PDF parameter A =      " + a + "");
            Console.WriteLine("  PDF parameter B =      " + b + "");

            if (!Burkardt.Probability.Half.half_normal_check(a, b))
            {
                Console.WriteLine("");
                Console.WriteLine("HALF_NORMAL_CDF_TEST - Fatal error!");
                Console.WriteLine("  The parameters are not legal.");
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("       X            PDF           CDF            CDF_INV");
            Console.WriteLine("");

            for (i = 1; i <= 10; i++)
            {
                x = Burkardt.Probability.Half.half_normal_sample(a, b, ref seed);
                pdf = Burkardt.Probability.Half.half_normal_pdf(x, a, b);
                cdf = Burkardt.Probability.Half.half_normal_cdf(x, a, b);
                x2 = Burkardt.Probability.Half.half_normal_cdf_inv(cdf, a, b);

                Console.WriteLine("  "
                                  + x.ToString().PadLeft(12) + "  "
                                  + pdf.ToString().PadLeft(12) + "  "
                                  + cdf.ToString().PadLeft(12) + "  "
                                  + x2.ToString().PadLeft(12) + "");
            }

            return;
        }

        static void half_normal_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    HALF_NORMAL_SAMPLE_TEST tests HALF_NORMAL_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 April 2016
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
            Console.WriteLine("HALF_NORMAL_SAMPLE_TEST");
            Console.WriteLine("  HALF_NORMAL_MEAN computes the Half Normal mean;");
            Console.WriteLine("  HALF_NORMAL_SAMPLE samples the Half Normal distribution;");
            Console.WriteLine("  HALF_NORMAL_VARIANCE computes the Half Normal variance;");

            a = 0.0;
            b = 10.0;

            Console.WriteLine("");
            Console.WriteLine("  PDF parameter A =      " + a + "");
            Console.WriteLine("  PDF parameter B =      " + b + "");

            if (!Burkardt.Probability.Half.half_normal_check(a, b))
            {
                Console.WriteLine("");
                Console.WriteLine("HALF_NORMAL_SAMPLE_TEST - Fatal error!");
                Console.WriteLine("  The parameters are not legal.");
                return;
            }

            mean = Burkardt.Probability.Half.half_normal_mean(a, b);
            variance = Burkardt.Probability.Half.half_normal_variance(a, b);

            Console.WriteLine("");
            Console.WriteLine("  PDF mean =     " + mean + "");
            Console.WriteLine("  PDF variance = " + variance + "");

            for (i = 0; i < SAMPLE_NUM; i++)
            {
                x[i] = Burkardt.Probability.Half.half_normal_sample(a, b, ref seed);
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