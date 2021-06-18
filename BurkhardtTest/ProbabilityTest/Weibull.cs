using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest
{
    partial class Program
    {
        static void weibull_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    WEIBULL_CDF_TEST tests WEIBULL_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 April 2016
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
            Console.WriteLine("WEIBULL_CDF_TEST");
            Console.WriteLine("  WEIBULL_CDF evaluates the Weibull CDF;");
            Console.WriteLine("  WEIBULL_CDF_INV inverts the Weibull CDF.");
            Console.WriteLine("  WEIBULL_PDF evaluates the Weibull PDF;");

            a = 2.0;
            b = 3.0;
            c = 4.0;

            Console.WriteLine("");
            Console.WriteLine("  PDF parameter A =      " + a + "");
            Console.WriteLine("  PDF parameter B =      " + b + "");
            Console.WriteLine("  PDF parameter C =      " + c + "");

            if (!Weibull.weibull_check(a, b, c))
            {
                Console.WriteLine("");
                Console.WriteLine("WEIBULL_CDF_TEST - Fatal error!");
                Console.WriteLine("  The parameters are not legal.");
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("       X            PDF           CDF            CDF_INV");
            Console.WriteLine("");

            for (i = 1; i <= 10; i++)
            {
                x = Weibull.weibull_sample(a, b, c, ref seed);
                pdf = Weibull.weibull_pdf(x, a, b, c);
                cdf = Weibull.weibull_cdf(x, a, b, c);
                x2 = Weibull.weibull_cdf_inv(cdf, a, b, c);

                Console.WriteLine("  "
                                  + x.ToString().PadLeft(12) + "  "
                                  + pdf.ToString().PadLeft(12) + "  "
                                  + cdf.ToString().PadLeft(12) + "  "
                                  + x2.ToString().PadLeft(12) + "");
            }

            return;
        }

        static void weibull_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    WEIBULL_SAMPLE_TEST tests WEIBULL_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 April 2016
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
            Console.WriteLine("WEIBULL_SAMPLE_TEST");
            Console.WriteLine("  WEIBULL_MEAN computes the Weibull mean;");
            Console.WriteLine("  WEIBULL_SAMPLE samples the Weibull distribution;");
            Console.WriteLine("  WEIBULL_VARIANCE computes the Weibull variance.");

            a = 2.0;
            b = 3.0;
            c = 4.0;

            Console.WriteLine("");
            Console.WriteLine("  PDF parameter A =      " + a + "");
            Console.WriteLine("  PDF parameter B =      " + b + "");
            Console.WriteLine("  PDF parameter C =      " + c + "");

            if (!Weibull.weibull_check(a, b, c))
            {
                Console.WriteLine("");
                Console.WriteLine("WEIBULL_SAMPLE_TEST - Fatal error!");
                Console.WriteLine("  The parameters are not legal.");
                return;
            }

            mean = Weibull.weibull_mean(a, b, c);
            variance = Weibull.weibull_variance(a, b, c);

            Console.WriteLine("");
            Console.WriteLine("  PDF mean =     " + mean + "");
            Console.WriteLine("  PDF variance = " + variance + "");

            for (i = 0; i < SAMPLE_NUM; i++)
            {
                x[i] = Weibull.weibull_sample(a, b, c, ref seed);
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

        static void weibull_discrete_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    WEIBULL_DISCRETE_CDF_TEST tests WEIBULL_DISCRETE_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 April 2016
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
            int x;
            int x2;

            Console.WriteLine("");
            Console.WriteLine("WEIBULL_DISCRETE_CDF_TEST");
            Console.WriteLine("  WEIBULL_DISCRETE_CDF evaluates the Weibull Discrete CDF;");
            Console.WriteLine("  WEIBULL_DISCRETE_CDF_INV inverts the Weibull Discrete CDF.");
            Console.WriteLine("  WEIBULL_DISCRETE_PDF evaluates the Weibull Discrete PDF;");

            a = 0.5;
            b = 1.5;

            Console.WriteLine("");
            Console.WriteLine("  PDF parameter A =      " + a + "");
            Console.WriteLine("  PDF parameter B =      " + b + "");

            if (!Weibull.weibull_discrete_check(a, b))
            {
                Console.WriteLine("");
                Console.WriteLine("WEIBULL_DISCRETE_CDF_TEST - Fatal error!");
                Console.WriteLine("  The parameters are not legal.");
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("       X            PDF           CDF            CDF_INV");
            Console.WriteLine("");

            for (i = 1; i <= 10; i++)
            {
                x = Weibull.weibull_discrete_sample(a, b, ref seed);
                pdf = Weibull.weibull_discrete_pdf(x, a, b);
                cdf = Weibull.weibull_discrete_cdf(x, a, b);
                x2 = Weibull.weibull_discrete_cdf_inv(cdf, a, b);

                Console.WriteLine("  "
                                  + x.ToString().PadLeft(12) + "  "
                                  + pdf.ToString().PadLeft(12) + "  "
                                  + cdf.ToString().PadLeft(12) + "  "
                                  + x2.ToString().PadLeft(12) + "");
            }

            return;
        }

        static void weibull_discrete_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    WEIBULL_DISCRETE_SAMPLE_TEST tests WEIBULL_DISCRETE_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 April 2016
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
            Console.WriteLine("WEIBULL_DISCRETE_SAMPLE_TEST");
            Console.WriteLine("  WEIBULL_DISCRETE_SAMPLE samples the Weibull Discrete distribution;");

            a = 0.5;
            b = 1.5;

            Console.WriteLine("");
            Console.WriteLine("  PDF parameter A =      " + a + "");
            Console.WriteLine("  PDF parameter B =      " + b + "");

            if (!Weibull.weibull_discrete_check(a, b))
            {
                Console.WriteLine("");
                Console.WriteLine("WEIBULL_DISCRETE_SAMPLE_TEST - Fatal error!");
                Console.WriteLine("  The parameters are not legal.");
                return;
            }

            for (i = 0; i < SAMPLE_NUM; i++)
            {
                x[i] = Weibull.weibull_discrete_sample(a, b, ref seed);
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