using System;
using Burkardt.Probability;
using Burkardt.Types;

namespace Burkardt.ProbabilityTest
{
    partial class Program
    {
        static void english_letter_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    ENGLISH_LETTER_CDF_TEST tests ENGLISH_LETTER_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 March 2016
//
//  Author:
//
//    John Burkardt
//
        {
            char c;
            char c2;
            double cdf;
            int i;
            double pdf;
            int seed;

            Console.WriteLine("");
            Console.WriteLine("ENGLISH_LETTER_CDF_TEST");
            Console.WriteLine("  ENGLISH_LETTER_CDF evaluates the English Letter CDF;");
            Console.WriteLine("  ENGLISH_LETTER_CDF_INV inverts the English Letter CDF.");
            Console.WriteLine("  ENGLISH_LETTER_PDF evaluates the English Letter PDF;");

            seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("   C              PDF             CDF    CDF_INV");
            Console.WriteLine("");

            for (i = 1; i <= 10; i++)
            {
                c = English.english_letter_sample(ref seed);
                pdf = English.english_letter_pdf(c);
                cdf = English.english_letter_cdf(c);
                c2 = English.english_letter_cdf_inv(cdf);

                Console.WriteLine("  '" + c + "'"
                                  + "  " + pdf.ToString().PadLeft(14)
                                  + "  " + cdf.ToString().PadLeft(14)
                                  + "  '" + c2 + "'");
            }
        }

        static void english_sentence_length_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    ENGLISH_SENTENCE_LENGTH_CDF_TEST tests ENGLISH_SENTENCE_LENGTH_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 April 2016
//
//  Author:
//
//    John Burkardt
//
        {
            double cdf;
            int i;
            double pdf;
            int seed = 123456789;
            int x;
            int x2;

            Console.WriteLine("");
            Console.WriteLine("ENGLISH_SENTENCE_LENGTH_CDF_TEST");
            Console.WriteLine("  ENGLISH_SENTENCE_LENGTH_CDF evaluates the English Sentence Length CDF;");
            Console.WriteLine("  ENGLISH_SENTENCE_LENGTH_CDF_INV inverts the English Sentence Length CDF.");
            Console.WriteLine("  ENGLISH_SENTENCE_LENGTH_PDF evaluates the English Sentence Length PDF;");

            Console.WriteLine("");
            Console.WriteLine("       X            PDF           CDF            CDF_INV");
            Console.WriteLine("");

            for (i = 1; i <= 10; i++)
            {
                x = English.english_sentence_length_sample(ref seed);

                pdf = English.english_sentence_length_pdf(x);

                cdf = English.english_sentence_length_cdf(x);

                x2 = English.english_sentence_length_cdf_inv(cdf);

                Console.WriteLine("  " + x.ToString().PadLeft(12)
                                       + "  " + pdf.ToString().PadLeft(12)
                                       + "  " + cdf.ToString().PadLeft(12)
                                       + "  " + x2.ToString().PadLeft(12) + "");
            }

            return;
        }

        static void english_sentence_length_sample_test()

//****************************************************************************80
//
//  Purpose:
//
//    ENGLISH_SENTENCE_LENGTH_SAMPLE_TEST tests ENGLISH_SENTENCE_LENGTH_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 April 2016
//
//  Author:
//
//    John Burkardt
//
        {
            int SAMPLE_NUM = 1000;

            int i;
            double mean;
            int seed = 123456789;
            double variance;
            int[] x = new int[SAMPLE_NUM];
            int xmax;
            int xmin;

            Console.WriteLine("");
            Console.WriteLine("ENGLISH_SENTENCE_LENGTH_SAMPLE_TEST");
            Console.WriteLine("  ENGLISH_SENTENCE_LENGTH_MEAN computes the English Sentence Length mean;");
            Console.WriteLine("  ENGLISH_SENTENCE_LENGTH_SAMPLE samples the English Sentence Length distribution;");
            Console.WriteLine("  ENGLISH_SENTENCE_LENGTH_VARIANCE computes the English Sentence Length variance.");

            mean = English.english_sentence_length_mean();
            variance = English.english_sentence_length_variance();

            Console.WriteLine("");
            Console.WriteLine("  PDF mean =                    " + mean + "");
            Console.WriteLine("  PDF variance =                " + variance + "");

            for (i = 0; i < SAMPLE_NUM; i++)
            {
                x[i] = English.english_sentence_length_sample(ref seed);
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

        static void english_word_length_cdf_test()

//****************************************************************************80
//
//  Purpose:
//
//    ENGLISH_WORD_LENGTH_CDF_TEST tests ENGLISH_WORD_LENGTH_CDF.
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
            double cdf;
            int i;
            double pdf;
            int seed = 123456789;
            int x;
            int x2;

            Console.WriteLine("");
            Console.WriteLine("ENGLISH_WORD_LENGTH_CDF_TEST");
            Console.WriteLine("  ENGLISH_WORD_LENGTH_CDF evaluates the English Word LengthCDF;");
            Console.WriteLine("  ENGLISH_WORD_LENGTH_CDF_INV inverts the English Word LengthCDF.");
            Console.WriteLine("  ENGLISH_WORD_LENGTH_PDF evaluates the English Word LengthPDF;");

            Console.WriteLine("");
            Console.WriteLine("       X            PDF           CDF            CDF_INV");
            Console.WriteLine("");

            for (i = 1; i <= 10; i++)
            {
                x = English.english_word_length_sample(ref seed);

                pdf = English.english_word_length_pdf(x);

                cdf = English.english_word_length_cdf(x);

                x2 = English.english_word_length_cdf_inv(cdf);

                Console.WriteLine("  " + x.ToString().PadLeft(12)
                                       + "  " + pdf.ToString().PadLeft(12)
                                       + "  " + cdf.ToString().PadLeft(12)
                                       + "  " + x2.ToString().PadLeft(12) + "");
            }

            return;
        }

        static void english_word_length_sample_test()
//****************************************************************************80
//
//  Purpose:
//
//    ENGLISH_WORD_LENGTH_SAMPLE_TEST tests ENGLISH_WORD_LENGTH_SAMPLE.
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

            int i;
            double mean;
            int seed = 123456789;
            double variance;
            int[] x = new int[SAMPLE_NUM];
            int xmax;
            int xmin;

            Console.WriteLine("");
            Console.WriteLine("ENGLISH_WORD_LENGTH_SAMPLE_TEST");
            Console.WriteLine("  ENGLISH_WORD_LENGTH_MEAN computes the English Word Lengthmean;");
            Console.WriteLine("  ENGLISH_WORD_LENGTH_SAMPLE samples the English Word Lengthdistribution;");
            Console.WriteLine("  ENGLISH_WORD_LENGTH_VARIANCE computes the English Word Lengthvariance.");

            mean = English.english_word_length_mean();
            variance = English.english_word_length_variance();

            Console.WriteLine("");
            Console.WriteLine("  PDF mean =                    " + mean + "");
            Console.WriteLine("  PDF variance =                " + variance + "");

            for (i = 0; i < SAMPLE_NUM; i++)
            {
                x[i] = English.english_word_length_sample(ref seed);
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