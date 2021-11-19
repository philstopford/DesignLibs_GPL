using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal partial class Program
{
    private static void buffon_box_pdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    BUFFON_BOX_PDF_TEST tests BUFFON_BOX_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 April 2016
//
//  Author:
//
//    John Burkardt
//
    {
        double a;
        double b;
        int i;
        int j;
        int k;
        double l;
        double pdf;

        Console.WriteLine("");
        Console.WriteLine("BUFFON_BOX_PDF_TEST tests BUFFON_BOX_PDF.");
        Console.WriteLine("  BUFFON_BOX_PDF evaluates the Buffon-Laplace PDF, the probability");
        Console.WriteLine("  that, on a grid of cells of width A and height B,");
        Console.WriteLine("  a needle of length L, dropped at random, will cross");
        Console.WriteLine("  at least one grid line.");
        Console.WriteLine("");
        Console.WriteLine("      A         B         L        PDF");
        Console.WriteLine("");

        for (i = 1; i <= 5; i++)
        {
            a = i;
            for (j = 1; j <= 5; j++)
            {
                b = j;
                for (k = 0; k <= 5; k++)
                {
                    l = k * Math.Min(a, b) / 5.0;
                    pdf = Buffon.buffon_box_pdf(a, b, l);
                    Console.WriteLine("  " + a.ToString().PadLeft(8)
                                           + "  " + b.ToString().PadLeft(8)
                                           + "  " + l.ToString().PadLeft(8)
                                           + "  " + pdf.ToString().PadLeft(14) + "");
                }

                Console.WriteLine("");
            }
        }
    }

    private static void buffon_box_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    BUFFON_BOX_SAMPLE_TEST tests BUFFON_BOX_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 April 2016
//
//  Author:
//
//    John Burkardt
//
    {
        int TEST_NUM = 4;

        double a;
        double b;
        double err;
        int hits;
        double l;
        double pi = 3.141592653589793238462643;
        double pi_est;
        int test;
        int trial_num;
        int[] trial_num_test =  {
                10, 100, 10000, 1000000
            }
            ;

        a = 1.0;
        b = 1.0;
        l = 1.0;

        Console.WriteLine("");
        Console.WriteLine("BUFFON_BOX_SAMPLE_TEST");
        Console.WriteLine("  BUFFON_BOX_SAMPLE simulates a Buffon-Laplace needle dropping");
        Console.WriteLine("  experiment.  On a grid of cells of width A and height B,");
        Console.WriteLine("  a needle of length L is dropped at random.  We count");
        Console.WriteLine("  the number of times it crosses at least one grid line,");
        Console.WriteLine("  and use this to estimate the value of PI.");

        Console.WriteLine("");
        Console.WriteLine("  Cell width A =    " + a + "");
        Console.WriteLine("  Cell height B =   " + b + "");
        Console.WriteLine("  Needle length L = " + l + "");
        Console.WriteLine("");
        Console.WriteLine("    Trials      Hits          Est(Pi)     Err");
        Console.WriteLine("");

        for (test = 0; test < TEST_NUM; test++)
        {
            trial_num = trial_num_test[test];

            hits = Buffon.buffon_box_sample(a, b, l, trial_num);

            pi_est = hits switch
            {
                > 0 => (2.0 * l * (a + b) - l * l) * trial_num / (a * b * hits),
                _ => typeMethods.r8_huge()
            };

            err = Math.Abs(pi_est - pi);

            Console.WriteLine("  " + trial_num.ToString().PadLeft(8)
                                   + "  " + hits.ToString().PadLeft(8)
                                   + "  " + pi_est.ToString().PadLeft(14)
                                   + "  " + err.ToString().PadLeft(14) + "");
        }

    }

    private static void buffon_pdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    BUFFON_PDF_TEST tests BUFFON_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 April 2016
//
//  Author:
//
//    John Burkardt
//
    {
        double a;
        int i;
        int k;
        double l;
        double pdf;

        Console.WriteLine("");
        Console.WriteLine("BUFFON_PDF_TEST");
        Console.WriteLine("  BUFFON_PDF evaluates the Buffon PDF, the probability");
        Console.WriteLine("  that, on a grid of cells of width A,");
        Console.WriteLine("  a needle of length L, dropped at random, will cross");
        Console.WriteLine("  at least one grid line.");
        Console.WriteLine("");
        Console.WriteLine("      A         L        PDF");
        Console.WriteLine("");

        for (i = 1; i <= 5; i++)
        {
            a = i;
            for (k = 0; k <= 5; k++)
            {
                l = k * a / 5.0;
                pdf = Buffon.buffon_pdf(a, l);
                Console.WriteLine("  " + a.ToString().PadLeft(8)
                                       + "  " + l.ToString().PadLeft(8)
                                       + "  " + pdf.ToString().PadLeft(14)+ "");
            }

            Console.WriteLine("");
        }
    }

    private static void buffon_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    BUFFON_SAMPLE_TEST tests BUFFON_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 April 2016
//
//  Author:
//
//    John Burkardt
//
    {
        int TEST_NUM = 4;

        double a;
        double err;
        int hits;
        double l;
        double pi = 3.141592653589793238462643;
        double pi_est;
        int test;
        int trial_num;
        int[] trial_num_test =  {
            10, 100, 10000, 1000000
        };

        a = 1.0;
        l = 1.0;

        Console.WriteLine("");
        Console.WriteLine("BUFFON_SAMPLE_TEST");
        Console.WriteLine("  BUFFON_SAMPLE simulates a Buffon needle dropping");
        Console.WriteLine("  experiment.  On a grid of cells of width A,");
        Console.WriteLine("  a needle of length L is dropped at random.  We count");
        Console.WriteLine("  the number of times it crosses at least one grid line,");
        Console.WriteLine("  and use this to estimate the value of PI.");

        Console.WriteLine("");
        Console.WriteLine("  Cell width A =    " + a + "");
        Console.WriteLine("  Needle length L = " + l + "");
        Console.WriteLine("");
        Console.WriteLine("    Trials      Hits          Est(Pi)     Err");
        Console.WriteLine("");

        for (test = 0; test < TEST_NUM; test++)
        {
            trial_num = trial_num_test[test];

            hits = Buffon.buffon_sample(a, l, trial_num);

            pi_est = hits switch
            {
                > 0 => 2.0 * l * trial_num / (a * hits),
                _ => typeMethods.r8_huge()
            };

            err = Math.Abs(pi_est - pi);

            Console.WriteLine("  " + trial_num.ToString().PadLeft(8)
                                   + "  " + hits.ToString().PadLeft(8)
                                   + "  " + pi_est.ToString().PadLeft(14)
                                   + "  " + err.ToString().PadLeft(14) + "");
        }
    }

}