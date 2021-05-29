using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace Burkardt.ProbabilityTest
{
    partial class Programs
    {
        static void runs_pdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    RUNS_PDF_TEST tests RUNS_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 April 2016
//
//  Author:
//
//    John Burkardt
//
        {
            int m;
            int n;
            double pdf;
            double pdf_total;
            int r;

            Console.WriteLine("");
            Console.WriteLine("RUNS_PDF_TEST");
            Console.WriteLine("  RUNS_PDF evaluates the Runs PDF;");
            Console.WriteLine("");
            Console.WriteLine("  M is the number of symbols of one kind,");
            Console.WriteLine("  N is the number of symbols of the other kind,");
            Console.WriteLine("  R is the number of runs (sequences of one symbol)");
            Console.WriteLine("");
            Console.WriteLine("         M         N         R      PDF");
            Console.WriteLine("");

            m = 6;

            for (n = 0; n <= 9; n++)
            {
                Console.WriteLine("");
                pdf_total = 0.0;

                for (r = 1; r <= 2 * Math.Min(m, n) + 2; r++)
                {
                    pdf = Runs.runs_pdf(m, n, r);

                    Console.WriteLine("  " + m.ToString().PadLeft(8)
                                      + "  " + n.ToString().PadLeft(8)
                                      + "  " + r.ToString().PadLeft(8)
                                      + "  " + pdf.ToString().PadLeft(14) + "");

                    pdf_total = pdf_total + pdf;
                }

                Console.WriteLine("  " + m.ToString().PadLeft(8)
                                  + "  " + "        "
                                  + "  " + "        "
                                  + "  " + pdf_total.ToString().PadLeft(14) + "");

            }

            return;
        }

        static void runs_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    RUNS_SAMPLE_TEST tests RUNS_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 April 2016
//
//  Author:
//
//    John Burkardt
//
        {
            int SAMPLE_NUM = 1000;

            int i;
            int m;
            double mean;
            int n;
            int seed = 123456789;
            double variance;
            int[] x = new int[SAMPLE_NUM];
            int xmax;
            int xmin;

            Console.WriteLine("");
            Console.WriteLine("RUNS_SAMPLE_TEST");
            Console.WriteLine("  RUNS_MEAN computes the Runs mean;");
            Console.WriteLine("  RUNS_SAMPLE samples the Runs distribution;");
            Console.WriteLine("  RUNS_VARIANCE computes the Runs variance");

            m = 10;
            n = 5;

            Console.WriteLine("");
            Console.WriteLine("  PDF parameter M = " + m + "");
            Console.WriteLine("  PDF parameter N = " + n + "");

            mean = Runs.runs_mean(m, n);
            variance = Runs.runs_variance(m, n);

            Console.WriteLine("  PDF mean =        " + mean + "");
            Console.WriteLine("  PDF variance =    " + variance + "");

            for (i = 0; i < SAMPLE_NUM; i++)
            {
                x[i] = Runs.runs_sample(m, n, ref seed);
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