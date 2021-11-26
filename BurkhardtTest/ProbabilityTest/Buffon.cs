using System;
using System.Globalization;
using Burkardt.Probability;
using Burkardt.Types;

namespace ProbabilityTest;

internal static partial class Program
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
        int i;

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
            double a = i;
            int j;
            for (j = 1; j <= 5; j++)
            {
                double b = j;
                int k;
                for (k = 0; k <= 5; k++)
                {
                    double l = k * Math.Min(a, b) / 5.0;
                    double pdf = Buffon.buffon_box_pdf(a, b, l);
                    Console.WriteLine("  " + a.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                           + "  " + b.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                           + "  " + l.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                           + "  " + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
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
        const int TEST_NUM = 4;

        int test;
        int[] trial_num_test =  {
                10, 100, 10000, 1000000
            }
            ;

        const double a = 1.0;
        const double b = 1.0;
        const double l = 1.0;

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
            int trial_num = trial_num_test[test];

            int hits = Buffon.buffon_box_sample(a, b, l, trial_num);

            double pi_est = hits switch
            {
                > 0 => (2.0 * l * (a + b) - l * l) * trial_num / (a * b * hits),
                _ => typeMethods.r8_huge()
            };

            double err = Math.Abs(pi_est - Math.PI);

            Console.WriteLine("  " + trial_num.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + hits.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + pi_est.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
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
        int i;

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
            double a = i;
            int k;
            for (k = 0; k <= 5; k++)
            {
                double l = k * a / 5.0;
                double pdf = Buffon.buffon_pdf(a, l);
                Console.WriteLine("  " + a.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + l.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + pdf.ToString(CultureInfo.InvariantCulture).PadLeft(14)+ "");
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
        const int TEST_NUM = 4;

        int test;
        int[] trial_num_test =  {
            10, 100, 10000, 1000000
        };

        const double a = 1.0;
        const double l = 1.0;

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
            int trial_num = trial_num_test[test];

            int hits = Buffon.buffon_sample(a, l, trial_num);

            double pi_est = hits switch
            {
                > 0 => 2.0 * l * trial_num / (a * hits),
                _ => typeMethods.r8_huge()
            };

            double err = Math.Abs(pi_est - Math.PI);

            Console.WriteLine("  " + trial_num.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + hits.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + pi_est.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

}