using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace Burkardt.ProbabilityTest
{
    partial class Program
    {
        static void cardioid_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    CARDIOID_CDF_TEST tests CARDIOID_CDF, CARDIOID_CDF_INV, CARDIOID_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 July 2005
//
//  Author:
//
//    John Burkardt
//
        {
            double a = 0.0;
            double b = 0.25;
            double cdf;
            int i;
            double pdf;
            int seed = 123456789;
            double x;
            double x2;

            Console.WriteLine("");
            Console.WriteLine("CARDIOID_CDF_TEST");
            Console.WriteLine("  CARDIOID_CDF evaluates the Cardioid CDF;");
            Console.WriteLine("  CARDIOID_CDF_INV inverts the Cardioid CDF.");
            Console.WriteLine("  CARDIOID_PDF evaluates the Cardioid PDF;");

            Console.WriteLine("");
            Console.WriteLine("  PDF parameter A = " + a + "");
            Console.WriteLine("  PDF parameter B = " + b + "");

            if (!Cardioid.cardioid_check(a, b))
            {
                Console.WriteLine("");
                Console.WriteLine("CARDIOID_CDF_TEST - Fatal error!");
                Console.WriteLine("  The parameters are not legal.");
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("       X            PDF           CDF            CDF_INV");
            Console.WriteLine("");

            for (i = 0; i < 10; i++)
            {
                x = Cardioid.cardioid_sample(a, b, ref seed);
                pdf = Cardioid.cardioid_pdf(x, a, b);
                cdf = Cardioid.cardioid_cdf(x, a, b);
                x2 = Cardioid.cardioid_cdf_inv(cdf, a, b);
                Console.WriteLine("  " + x.ToString().PadLeft(12)
                                  + "  " + pdf.ToString().PadLeft(12)
                                  + "  " + cdf.ToString().PadLeft(12)
                                  + "  " + x2.ToString().PadLeft(12) + "");
            }
        }

        static void cardioid_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    CARDIOID_SAMPLE_TEST tests CARDIOID_MEAN, CARDIOID_SAMPLE, CARDIOID_VARIANCE.
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

            double a = 0.0;
            double b = 0.25;
            int i;
            double mean;
            int seed = 123456789;
            double variance;
            double[] x = new double [SAMPLE_NUM];
            double xmax;
            double xmin;

            Console.WriteLine("");
            Console.WriteLine("CARDIOID_SAMPLE_TEST");
            Console.WriteLine("  CARDIOID_MEAN computes the Cardioid mean;");
            Console.WriteLine("  CARDIOID_SAMPLE samples the Cardioid distribution;");
            Console.WriteLine("  CARDIOID_VARIANCE computes the Cardioid variance.");

            Console.WriteLine("");
            Console.WriteLine("  PDF parameter A = " + a + "");
            Console.WriteLine("  PDF parameter B = " + b + "");

            if (!Cardioid.cardioid_check(a, b))
            {
                Console.WriteLine("");
                Console.WriteLine("CARDIOID_SAMPLE_TEST - Fatal error!");
                Console.WriteLine("  The parameters are not legal.");
                return;
            }

            mean = Cardioid.cardioid_mean(a, b);
            variance = Cardioid.cardioid_variance(a, b);

            Console.WriteLine("");
            Console.WriteLine("  PDF mean =                    " + mean + "");
            Console.WriteLine("  PDF variance =                " + variance + "");

            for (i = 0; i < SAMPLE_NUM; i++)
            {
                x[i] = Cardioid.cardioid_sample(a, b, ref seed);
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