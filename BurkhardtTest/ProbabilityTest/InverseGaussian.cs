using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace Burkardt.ProbabilityTest
{
    partial class Program
    {
        static void inverse_gaussian_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    INVERSE_GAUSSIAN_CDF_TEST tests INVERSE_GAUSSIAN_CDF.
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

            Console.WriteLine("");
            Console.WriteLine("INVERSE_GAUSSIAN_CDF_TEST");
            Console.WriteLine("  INVERSE_GAUSSIAN_CDF evaluates the Inverse Gaussian CDF;");
            Console.WriteLine("  INVERSE_GAUSSIAN_PDF evaluates the Inverse Gaussian PDF;");

            a = 5.0;
            b = 2.0;

            Console.WriteLine("");
            Console.WriteLine("  PDF parameter A =      " + a + "");
            Console.WriteLine("  PDF parameter B =      " + b + "");

            if (!InverseGaussian.inverse_gaussian_check(a, b))
            {
                Console.WriteLine("");
                Console.WriteLine("INVERSE_GAUSSIAN_CDF_TEST - Fatal error!");
                Console.WriteLine("  The parameters are not legal.");
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("       X            PDF           CDF");
            Console.WriteLine("");

            for (i = 1; i <= 10; i++)
            {
                x = InverseGaussian.inverse_gaussian_sample(a, b, ref seed);
                pdf = InverseGaussian.inverse_gaussian_pdf(x, a, b);
                cdf = InverseGaussian.inverse_gaussian_cdf(x, a, b);

                Console.WriteLine("  "
                                  + x.ToString().PadLeft(12) + "  "
                                  + pdf.ToString().PadLeft(12) + "  "
                                  + cdf.ToString().PadLeft(12) + "");
            }

            return;
        }

        static void inverse_gaussian_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    INVERSE_GAUSSIAN_SAMPLE_TEST tests INVERSE_GAUSSIAN_SAMPLE.
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
            Console.WriteLine("INVERSE_GAUSSIAN_SAMPLE_TEST");
            Console.WriteLine("  INVERSE_GAUSSIAN_MEAN computes the Inverse Gaussian mean;");
            Console.WriteLine("  INVERSE_GAUSSIAN_SAMPLE samples the Inverse Gaussian distribution;");
            Console.WriteLine("  INVERSE_GAUSSIAN_VARIANCE computes the Inverse Gaussian variance;");

            a = 2.0;
            b = 3.0;

            Console.WriteLine("");
            Console.WriteLine("  PDF parameter A =      " + a + "");
            Console.WriteLine("  PDF parameter B =      " + b + "");

            if (!InverseGaussian.inverse_gaussian_check(a, b))
            {
                Console.WriteLine("");
                Console.WriteLine("INVERSE_GAUSSIAN_SAMPLE_TEST - Fatal error!");
                Console.WriteLine("  The parameters are not legal.");
                return;
            }

            mean = InverseGaussian.inverse_gaussian_mean(a, b);
            variance = InverseGaussian.inverse_gaussian_variance(a, b);

            Console.WriteLine("");
            Console.WriteLine("  PDF mean =     " + mean + "");
            Console.WriteLine("  PDF variance = " + variance + "");

            for (i = 0; i < SAMPLE_NUM; i++)
            {
                x[i] = InverseGaussian.inverse_gaussian_sample(a, b, ref seed);
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