using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest
{
    partial class Program
    {
        static void triangle_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_CDF_TEST tests TRIANGLE_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
        {
            double a;
            double b;
            double c;
            double cdf;
            int i;
            double pdf;
            int seed = 123456789;
            double x;
            double x2;

            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_CDF_TEST");
            Console.WriteLine("  TRIANGLE_CDF evaluates the Triangle CDF;");
            Console.WriteLine("  TRIANGLE_CDF_INV inverts the Triangle CDF.");
            Console.WriteLine("  TRIANGLE_PDF evaluates the Triangle PDF;");

            a = 1.0;
            b = 3.0;
            c = 10.0;

            Console.WriteLine("");
            Console.WriteLine("  PDF parameter A =      " + a + "");
            Console.WriteLine("  PDF parameter B =      " + b + "");
            Console.WriteLine("  PDF parameter C =      " + c + "");

            if (!Triangle.triangle_check(a, b, c))
            {
                Console.WriteLine("");
                Console.WriteLine("TRIANGLE_CDF_TEST - Fatal error!");
                Console.WriteLine("  The parameters are not legal.");
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("       X            PDF           CDF            CDF_INV");
            Console.WriteLine("");

            for (i = 1; i <= 10; i++)
            {
                x = Triangle.triangle_sample(a, b, c, ref seed);
                pdf = Triangle.triangle_pdf(x, a, b, c);
                cdf = Triangle.triangle_cdf(x, a, b, c);
                x2 = Triangle.triangle_cdf_inv(cdf, a, b, c);

                Console.WriteLine("  "
                                  + x.ToString().PadLeft(12) + "  "
                                  + pdf.ToString().PadLeft(12) + "  "
                                  + cdf.ToString().PadLeft(12) + "  "
                                  + x2.ToString().PadLeft(12) + "");
            }

        }

        static void triangle_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_SAMPLE_TEST tests TRIANGLE_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 March 2016
//
//  Author:
//
//    John Burkardt
//
        {
            int SAMPLE_NUM = 1000;

            double a;
            double b;
            double c;
            int i;
            double mean;
            int seed = 123456789;
            double variance;
            double[] x = new double [SAMPLE_NUM];
            double xmax;
            double xmin;

            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_SAMPLE_TEST");
            Console.WriteLine("  TRIANGLE_MEAN computes the Triangle mean;");
            Console.WriteLine("  TRIANGLE_SAMPLE samples the Triangle distribution;");
            Console.WriteLine("  TRIANGLE_VARIANCE computes the Triangle variance;");

            a = 1.0;
            b = 3.0;
            c = 10.0;

            Console.WriteLine("");
            Console.WriteLine("  PDF parameter A =        " + a + "");
            Console.WriteLine("  PDF parameter B =        " + b + "");
            Console.WriteLine("  PDF parameter C =        " + c + "");

            if (!Triangle.triangle_check(a, b, c))
            {
                Console.WriteLine("");
                Console.WriteLine("TRIANGLE_SAMPLE_TEST - Fatal error!");
                Console.WriteLine("  The parameters are not legal.");
                return;
            }

            mean = Triangle.triangle_mean(a, b, c);
            variance = Triangle.triangle_variance(a, b, c);

            Console.WriteLine("");
            Console.WriteLine("  PDF mean =     " + mean + "");
            Console.WriteLine("  PDF variance = " + variance + "");

            for (i = 0; i < SAMPLE_NUM; i++)
            {
                x[i] = Triangle.triangle_sample(a, b, c, ref seed);
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