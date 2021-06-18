using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest
{
    partial class Program
    {
            static void zipf_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    ZIPF_CDF_TEST tests ZIPF_CDF, ZIPF_CDF_INV, ZIPF_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 March 2016
//
//  Author:
//
//    John Burkardt
//
            {
                double a;
                double cdf;
                double pdf;
                int x;
                int x2;

                Console.WriteLine("");
                Console.WriteLine("ZIPF_CDF_TEST");
                Console.WriteLine("  ZIPF_CDF evaluates the Zipf CDF;");
                Console.WriteLine("  ZIPF_CDF_INV inverts the Zipf CDF;");
                Console.WriteLine("  ZIPF_PDF evaluates the Zipf PDF;");

                a = 2.0E+00;

                Console.WriteLine("");
                Console.WriteLine("  PDF parameter A =             " + a + "");

                if (!Zipf.zipf_check(a))
                {
                    Console.WriteLine("");
                    Console.WriteLine("ZIPF_CDF_TEST - Fatal error!");
                    Console.WriteLine("  The parameters are not legal.");
                    return;
                }

                Console.WriteLine("");
                Console.WriteLine("       X            PDF           CDF    CDF_INV()");
                Console.WriteLine("");

                for (x = 1; x <= 20; x++)
                {
                    pdf = Zipf.zipf_pdf(x, a);
                    cdf = Zipf.zipf_cdf(x, a);
                    x2 = Zipf.zipf_cdf_inv(a, cdf);

                    Console.WriteLine("  "
                                      + x.ToString().PadLeft(12) + "  "
                                      + pdf.ToString().PadLeft(12) + "  "
                                      + cdf.ToString().PadLeft(12) + "  "
                                      + x2.ToString().PadLeft(12) + "");
                }
            }

            static void zipf_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    ZIPF_SAMPLE_TEST tests ZIPF_MEAN, ZIPF_SAMPLE, ZIPF_VARIANCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 March 2016
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
                int[] x = new int[SAMPLE_NUM];
                int xmax;
                int xmin;

                Console.WriteLine("");
                Console.WriteLine("ZIPF_SAMPLE_TEST");
                Console.WriteLine("  ZIPF_MEAN computes the Zipf mean;");
                Console.WriteLine("  ZIPF_SAMPLE samples the Zipf distribution;");
                Console.WriteLine("  ZIPF_VARIANCE computes the Zipf variance.");

                a = 4.0E+00;

                Console.WriteLine("");
                Console.WriteLine("  PDF parameter A =             " + a + "");

                if (!Zipf.zipf_check(a))
                {
                    Console.WriteLine("");
                    Console.WriteLine("ZIPF_SAMPLE_TEST - Fatal error!");
                    Console.WriteLine("  The parameters are not legal.");
                    return;
                }

                mean = Zipf.zipf_mean(a);
                variance = Zipf.zipf_variance(a);

                Console.WriteLine("  PDF mean =                    " + mean + "");
                Console.WriteLine("  PDF variance =                " + variance + "");

                for (j = 0; j < SAMPLE_NUM; j++)
                {
                    x[j] = Zipf.zipf_sample(a, ref seed);
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