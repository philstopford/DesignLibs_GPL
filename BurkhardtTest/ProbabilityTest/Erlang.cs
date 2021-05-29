using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace Burkardt.ProbabilityTest
{
    partial class Program
    {
        static void erlang_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    ERLANG_CDF_TEST tests ERLANG_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 March 2016
//
//  Author:
//
//    John Burkardt
//
        {
            double a;
            double b;
            int c;
            double cdf;
            int i;
            double pdf;
            int seed = 123456789;
            double x;
            double x2;

            Console.WriteLine("");
            Console.WriteLine("ERLANG_CDF_TEST");
            Console.WriteLine("  ERLANG_CDF evaluates the Erlang CDF;");
            Console.WriteLine("  ERLANG_CDF_INV inverts the Erlang CDF.");
            Console.WriteLine("  ERLANG_PDF evaluates the Erlang PDF;");

            a = 1.0;
            b = 2.0;
            c = 3;

            Console.WriteLine("");
            Console.WriteLine("  PDF parameter A =      " + a + "");
            Console.WriteLine("  PDF parameter B =      " + b + "");
            Console.WriteLine("  PDF parameter C =      " + c + "");

            if (!Erlang.erlang_check(a, b, c))
            {
                Console.WriteLine("");
                Console.WriteLine("ERLANG_CDF_TEST - Fatal error!");
                Console.WriteLine("  The parameters are not legal.");
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("       X            PDF           CDF            CDF_INV");
            Console.WriteLine("");

            for (i = 1; i <= 10; i++)
            {
                x = Erlang.erlang_sample(a, b, c, ref seed);
                pdf = Erlang.erlang_pdf(x, a, b, c);
                cdf = Erlang.erlang_cdf(x, a, b, c);
                x2 = Erlang.erlang_cdf_inv(cdf, a, b, c);

                Console.WriteLine("  "
                                  + x.ToString().PadLeft(12) + "  "
                                  + pdf.ToString().PadLeft(12) + "  "
                                  + cdf.ToString().PadLeft(12) + "  "
                                  + x2.ToString().PadLeft(12) + "");
            }

        }

        static void erlang_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    ERLANG_SAMPLE_TEST tests ERLANG_SAMPLE.
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
            double b;
            int c;
            int i;
            double mean;
            int seed = 123456789;
            double variance;
            double[] x = new double [SAMPLE_NUM];
            double xmax;
            double xmin;

            Console.WriteLine("");
            Console.WriteLine("ERLANG_SAMPLE_TEST");
            Console.WriteLine("  ERLANG_MEAN computes the Erlang mean;");
            Console.WriteLine("  ERLANG_SAMPLE samples the Erlang distribution;");
            Console.WriteLine("  ERLANG_VARIANCE computes the Erlang variance;");

            a = 1.0;
            b = 2.0;
            c = 3;

            Console.WriteLine("");
            Console.WriteLine("  PDF parameter A =      " + a + "");
            Console.WriteLine("  PDF parameter B =      " + b + "");
            Console.WriteLine("  PDF parameter C =      " + c + "");

            if (!Erlang.erlang_check(a, b, c))
            {
                Console.WriteLine("");
                Console.WriteLine("ERLANG_SAMPLE_TEST - Fatal error!");
                Console.WriteLine("  The parameters are not legal.");
                return;
            }

            mean = Erlang.erlang_mean(a, b, c);
            variance = Erlang.erlang_variance(a, b, c);

            Console.WriteLine("");
            Console.WriteLine("  PDF mean =     " + mean + "");
            Console.WriteLine("  PDF variance = " + variance + "");

            for (i = 0; i < SAMPLE_NUM; i++)
            {
                x[i] = Erlang.erlang_sample(a, b, c, ref seed);
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