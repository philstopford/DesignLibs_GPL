using System;
using Burkardt.Uniform;

namespace Burkardt.Probability;

public static class English
{
    public static double english_letter_cdf(char c)
        //*****************************************************************************/
        //
        //  Purpose:
        //
        //    ENGLISH_LETTER_CDF evaluates the English Letter CDF.
        //
        //  Discussion:
        //
        //    CDF('c') = 0.12441
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
        //  Reference:
        //
        //    Robert Lewand,
        //    Cryptological Mathematics,
        //    Mathematics Association of America, 2000,
        //    ISBN13: 978-0883857199
        //
        //  Parameters:
        //
        //    Input, char C, the letter whose probability is desired.
        //    'a' <= c <= 'z', but case is ignored.
        //
        //    Output, double ENGLISH_LETTER_CDF, the probability that a random letter is less
        //    than or equal to C.
        //
    {
        double cdf;
        double[] cdf_vec = {
                0.00000,
                0.08167, 0.09659, 0.12441, 0.16694, 0.29396,
                0.31624, 0.33639, 0.39733, 0.46699, 0.46852,
                0.47624, 0.51649, 0.54055, 0.60804, 0.68311,
                0.70240, 0.70335, 0.76322, 0.82649, 0.91705,
                0.94463, 0.95441, 0.97802, 0.97952, 0.99926,
                1.00000
            }
            ;
        int i;

        switch (c)
        {
            case >= 'a' and <= 'z':
                i = c - 'a' + 1;
                cdf = cdf_vec[i];
                break;
            case >= 'A' and <= 'Z':
                i = c - 'A' + 1;
                cdf = cdf_vec[i];
                break;
            default:
                cdf = 0.0;
                break;
        }

        return cdf;
    }

    public static char english_letter_cdf_inv(double cdf)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ENGLISH_LETTER_CDF_INV inverts the English Letter CDF.
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
        //  Reference:
        //
        //    Robert Lewand,
        //    Cryptological Mathematics,
        //    Mathematics Association of America, 2000,
        //    ISBN13: 978-0883857199
        //
        //  Parameters:
        //
        //    Input, double CDF, a cumulative probability between 0 and 1.
        //
        //    Input, char ENGLISH_LETTER_CDF_INV, the corresponding letter.
        //
    {
        double[] cdf_vec = {
            0.00000,
            0.08167, 0.09659, 0.12441, 0.16694, 0.29396,
            0.31624, 0.33639, 0.39733, 0.46699, 0.46852,
            0.47624, 0.51649, 0.54055, 0.60804, 0.68311,
            0.70240, 0.70335, 0.76322, 0.82649, 0.91705,
            0.94463, 0.95441, 0.97802, 0.97952, 0.99926,
            1.00000
        };

        char c = ' ';

        for (int i = 1; i < 27; i++)
        {
            if (cdf <= cdf_vec[i])
            {
                c = (char)('a' + (i - 1));
                break;
            }
        }

        return c;
    }

    public static double english_letter_pdf(char c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ENGLISH_LETTER_PDF evaluates the English Letter PDF.
        //
        //  Discussion:
        //
        //    PDF('c') = 0.02782
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
        //  Reference:
        //
        //    Robert Lewand,
        //    Cryptological Mathematics,
        //    Mathematics Association of America, 2000,
        //    ISBN13: 978-0883857199
        //
        //  Parameters:
        //
        //    Input, char C, the letter whose probability is desired.
        //    'a' <= c <= 'z', but case is ignored.
        //
        //    Output, double ENGLISH_LETTER_PDF, the value of the PDF.
        //
    {
        int i;
        double pdf;
        double[] pdf_vec =  {
                0.08167, 0.01492, 0.02782, 0.04253, 0.12702,
                0.02228, 0.02015, 0.06094, 0.06966, 0.00153,
                0.00772, 0.04025, 0.02406, 0.06749, 0.07507,
                0.01929, 0.00095, 0.05987, 0.06327, 0.09056,
                0.02758, 0.00978, 0.02361, 0.00150, 0.01974,
                0.00074
            }
            ;

        switch (c)
        {
            case >= 'a' and <= 'z':
                i = c - 'a';
                pdf = pdf_vec[i];
                break;
            case >= 'A' and <= 'Z':
                i = c - 'A';
                pdf = pdf_vec[i];
                break;
            default:
                pdf = 0.0;
                break;
        }

        return pdf;
    }

    public static char english_letter_sample(ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ENGLISH_LETTER_SAMPLE samples the English Letter PDF.
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
        //  Parameters:
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, character C, a sample of the PDF.
        //
    {
        double cdf = UniformRNG.r8_uniform_01(ref seed);

        char c = english_letter_cdf_inv(cdf);

        return c;
    }

    public static double english_sentence_length_cdf(int x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ENGLISH_SENTENCE_LENGTH_CDF evaluates the English Sentence Length CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 August 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Henry Kucera, Winthrop Francis,
        //    Computational Analysis of Present-Day American English,
        //    Brown University Press, 1967.
        //
        //  Parameters:
        //
        //    Input, int X, the sentence length whose CDF is desired.
        //
        //    Output, double ENGLISH_SENTENCE_LENGTH_CDF, the value of the CDF.
        //
    {
        int SENTENCE_LENGTH_MAX = 79;

        double cdf = 0;
        double[] pdf_vec =  {
                0.00806,
                0.01370,
                0.01862,
                0.02547,
                0.03043,
                0.03189,
                0.03516,
                0.03545,
                0.03286,
                0.03533,
                0.03562,
                0.03788,
                0.03669,
                0.03751,
                0.03518,
                0.03541,
                0.03434,
                0.03305,
                0.03329,
                0.03103,
                0.02867,
                0.02724,
                0.02647,
                0.02526,
                0.02086,
                0.02178,
                0.02128,
                0.01801,
                0.01690,
                0.01556,
                0.01512,
                0.01326,
                0.01277,
                0.01062,
                0.01051,
                0.00901,
                0.00838,
                0.00764,
                0.00683,
                0.00589,
                0.00624,
                0.00488,
                0.00477,
                0.00406,
                0.00390,
                0.00350,
                0.00318,
                0.00241,
                0.00224,
                0.00220,
                0.00262,
                0.00207,
                0.00174,
                0.00174,
                0.00128,
                0.00121,
                0.00103,
                0.00117,
                0.00124,
                0.00082,
                0.00088,
                0.00061,
                0.00061,
                0.00075,
                0.00063,
                0.00056,
                0.00052,
                0.00057,
                0.00031,
                0.00029,
                0.00021,
                0.00017,
                0.00021,
                0.00034,
                0.00031,
                0.00011,
                0.00011,
                0.00008,
                0.00006
            }
            ;
        double pdf_sum = 0.99768;

        switch (x)
        {
            case < 1:
                cdf = 0.0;
                break;
            default:
            {
                if (x < SENTENCE_LENGTH_MAX)
                {
                    cdf = 0.0;
                    for (int i = 0; i < x; i++)
                    {
                        cdf += pdf_vec[i];
                    }

                    cdf /= pdf_sum;
                }
                else if (SENTENCE_LENGTH_MAX <= x)
                {
                    cdf = 1.0;
                }

                break;
            }
        }

        return cdf;
    }

    public static int english_sentence_length_cdf_inv(double cdf)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ENGLISH_SENTENCE_LENGTH_CDF_INV inverts the English Sentence Length CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 August 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Henry Kucera, Winthrop Francis,
        //    Computational Analysis of Present-Day American English,
        //    Brown University Press, 1967.
        //
        //  Parameters:
        //
        //    Input, double CDF, the value of the CDF.
        //    0.0 <= CDF <= 1.0.
        //
        //    Output, int ENGLISH_SENTENCE_LENGTH_CDF_INV, the corresponding sentence
        //    length for which CDF(X-1) < CDF <= CDF(X)
        //
    {
        int SENTENCE_LENGTH_MAX = 79;

        double[] pdf_vec =  {
            0.00806,
            0.01370,
            0.01862,
            0.02547,
            0.03043,
            0.03189,
            0.03516,
            0.03545,
            0.03286,
            0.03533,
            0.03562,
            0.03788,
            0.03669,
            0.03751,
            0.03518,
            0.03541,
            0.03434,
            0.03305,
            0.03329,
            0.03103,
            0.02867,
            0.02724,
            0.02647,
            0.02526,
            0.02086,
            0.02178,
            0.02128,
            0.01801,
            0.01690,
            0.01556,
            0.01512,
            0.01326,
            0.01277,
            0.01062,
            0.01051,
            0.00901,
            0.00838,
            0.00764,
            0.00683,
            0.00589,
            0.00624,
            0.00488,
            0.00477,
            0.00406,
            0.00390,
            0.00350,
            0.00318,
            0.00241,
            0.00224,
            0.00220,
            0.00262,
            0.00207,
            0.00174,
            0.00174,
            0.00128,
            0.00121,
            0.00103,
            0.00117,
            0.00124,
            0.00082,
            0.00088,
            0.00061,
            0.00061,
            0.00075,
            0.00063,
            0.00056,
            0.00052,
            0.00057,
            0.00031,
            0.00029,
            0.00021,
            0.00017,
            0.00021,
            0.00034,
            0.00031,
            0.00011,
            0.00011,
            0.00008,
            0.00006
        };
        double pdf_sum = 0.99768;
        int x;

        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine("");
                Console.WriteLine("ENGLISH_SENTENCE_LENGTH_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
        }

        double cum = 0.0;

        for (int j = 0; j < SENTENCE_LENGTH_MAX; j++)
        {
            cum += pdf_vec[j];

            if (cdf <= cum / pdf_sum)
            {
                x = j + 1;
                return x;
            }
        }

        x = SENTENCE_LENGTH_MAX;

        return x;
    }

    public static double english_sentence_length_mean()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ENGLISH_SENTENCE_LENGTH_MEAN evaluates the mean of the English Sentence Length PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 August 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Henry Kucera, Winthrop Francis,
        //    Computational Analysis of Present-Day American English,
        //    Brown University Press, 1967.
        //
        //  Parameters:
        //
        //    Output, double ENGLISH_SENTENCE_LENGTH_MEAN, the mean of the PDF.
        //
    {
        int SENTENCE_LENGTH_MAX = 79;

        double[] pdf_vec =  {
            0.00806,
            0.01370,
            0.01862,
            0.02547,
            0.03043,
            0.03189,
            0.03516,
            0.03545,
            0.03286,
            0.03533,
            0.03562,
            0.03788,
            0.03669,
            0.03751,
            0.03518,
            0.03541,
            0.03434,
            0.03305,
            0.03329,
            0.03103,
            0.02867,
            0.02724,
            0.02647,
            0.02526,
            0.02086,
            0.02178,
            0.02128,
            0.01801,
            0.01690,
            0.01556,
            0.01512,
            0.01326,
            0.01277,
            0.01062,
            0.01051,
            0.00901,
            0.00838,
            0.00764,
            0.00683,
            0.00589,
            0.00624,
            0.00488,
            0.00477,
            0.00406,
            0.00390,
            0.00350,
            0.00318,
            0.00241,
            0.00224,
            0.00220,
            0.00262,
            0.00207,
            0.00174,
            0.00174,
            0.00128,
            0.00121,
            0.00103,
            0.00117,
            0.00124,
            0.00082,
            0.00088,
            0.00061,
            0.00061,
            0.00075,
            0.00063,
            0.00056,
            0.00052,
            0.00057,
            0.00031,
            0.00029,
            0.00021,
            0.00017,
            0.00021,
            0.00034,
            0.00031,
            0.00011,
            0.00011,
            0.00008,
            0.00006
        };
        double pdf_sum = 0.99768;

        double mean = 0.0;
        for (int j = 1; j <= SENTENCE_LENGTH_MAX; j++)
        {
            mean += j * pdf_vec[j - 1];
        }

        mean /= pdf_sum;

        return mean;
    }

    public static double english_sentence_length_pdf(int x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ENGLISH_SENTENCE_LENGTH_PDF evaluates the English Sentence Length PDF.
        //
        //  Discussion:
        //
        //    PDF(A,B;X) = B(X) if 1 <= X <= A
        //                = 0    otherwise
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 August 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Henry Kucera, Winthrop Francis,
        //    Computational Analysis of Present-Day American English,
        //    Brown University Press, 1967.
        //
        //  Parameters:
        //
        //    Input, int X, the sentence length whose probability is desired.
        //
        //    Output, double ENGLISH_SENTENCE_LENGTH_PDF, the value of the PDF.
        //
    {
        int SENTENCE_LENGTH_MAX = 79;

        double pdf;
        double[] pdf_vec =  {
            0.00806,
            0.01370,
            0.01862,
            0.02547,
            0.03043,
            0.03189,
            0.03516,
            0.03545,
            0.03286,
            0.03533,
            0.03562,
            0.03788,
            0.03669,
            0.03751,
            0.03518,
            0.03541,
            0.03434,
            0.03305,
            0.03329,
            0.03103,
            0.02867,
            0.02724,
            0.02647,
            0.02526,
            0.02086,
            0.02178,
            0.02128,
            0.01801,
            0.01690,
            0.01556,
            0.01512,
            0.01326,
            0.01277,
            0.01062,
            0.01051,
            0.00901,
            0.00838,
            0.00764,
            0.00683,
            0.00589,
            0.00624,
            0.00488,
            0.00477,
            0.00406,
            0.00390,
            0.00350,
            0.00318,
            0.00241,
            0.00224,
            0.00220,
            0.00262,
            0.00207,
            0.00174,
            0.00174,
            0.00128,
            0.00121,
            0.00103,
            0.00117,
            0.00124,
            0.00082,
            0.00088,
            0.00061,
            0.00061,
            0.00075,
            0.00063,
            0.00056,
            0.00052,
            0.00057,
            0.00031,
            0.00029,
            0.00021,
            0.00017,
            0.00021,
            0.00034,
            0.00031,
            0.00011,
            0.00011,
            0.00008,
            0.00006
        };
        double pdf_sum = 0.99768;

        pdf = x switch
        {
            >= 1 when x <= SENTENCE_LENGTH_MAX => pdf_vec[x - 1] / pdf_sum,
            _ => 0.0
        };

        return pdf;
    }

    public static int english_sentence_length_sample(ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ENGLISH_SENTENCE_LENGTH_SAMPLE samples the English Sentence Length PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 August 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Henry Kucera, Winthrop Francis,
        //    Computational Analysis of Present-Day American English,
        //    Brown University Press, 1967.
        //
        //  Parameters:
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, int ENGLISH_SENTENCE_LENGTH_SAMPLE, a sample of the PDF.
        //
    {
        double cdf = UniformRNG.r8_uniform_01(ref seed);

        int x = english_sentence_length_cdf_inv(cdf);

        return x;
    }

    public static double english_sentence_length_variance()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ENGLISH_SENTENCE_LENGTH_VARIANCE: variance of the English Sentence Length PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 August 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Henry Kucera, Winthrop Francis,
        //    Computational Analysis of Present-Day American English,
        //    Brown University Press, 1967.
        //
        //  Parameters:
        //
        //    Output, double ENGLISH_SENTENCE_LENGTH_VARIANCE, the variance of the PDF.
        //
    {
        int SENTENCE_LENGTH_MAX = 79;

        double[] pdf_vec =  {
            0.00806,
            0.01370,
            0.01862,
            0.02547,
            0.03043,
            0.03189,
            0.03516,
            0.03545,
            0.03286,
            0.03533,
            0.03562,
            0.03788,
            0.03669,
            0.03751,
            0.03518,
            0.03541,
            0.03434,
            0.03305,
            0.03329,
            0.03103,
            0.02867,
            0.02724,
            0.02647,
            0.02526,
            0.02086,
            0.02178,
            0.02128,
            0.01801,
            0.01690,
            0.01556,
            0.01512,
            0.01326,
            0.01277,
            0.01062,
            0.01051,
            0.00901,
            0.00838,
            0.00764,
            0.00683,
            0.00589,
            0.00624,
            0.00488,
            0.00477,
            0.00406,
            0.00390,
            0.00350,
            0.00318,
            0.00241,
            0.00224,
            0.00220,
            0.00262,
            0.00207,
            0.00174,
            0.00174,
            0.00128,
            0.00121,
            0.00103,
            0.00117,
            0.00124,
            0.00082,
            0.00088,
            0.00061,
            0.00061,
            0.00075,
            0.00063,
            0.00056,
            0.00052,
            0.00057,
            0.00031,
            0.00029,
            0.00021,
            0.00017,
            0.00021,
            0.00034,
            0.00031,
            0.00011,
            0.00011,
            0.00008,
            0.00006
        };
        double pdf_sum = 0.99768;

        double mean = 0.0;
        for (int j = 1; j <= SENTENCE_LENGTH_MAX; j++)
        {
            mean += j * pdf_vec[j - 1];
        }

        mean /= pdf_sum;

        double variance = 0.0;
        for (int j = 1; j <= SENTENCE_LENGTH_MAX; j++)
        {
            variance += pdf_vec[j - 1] * Math.Pow(j - mean, 2);
        }

        variance /= pdf_sum;

        return variance;
    }

    public static double english_word_length_cdf(int x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ENGLISH_WORD_LENGTH_CDF evaluates the English Word Length CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 July 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Henry Kucera, Winthrop Francis,
        //    Computational Analysis of Present-Day American English,
        //    Brown University Press, 1967.
        //
        //  Parameters:
        //
        //    Input, int X, the word length whose CDF is desired.
        //
        //    Output, double ENGLISH_WORD_LENGTH_CDF, the value of the CDF.
        //
    {
        int WORD_LENGTH_MAX = 27;

        double cdf = 0;
        double[] pdf_vec =  {
                0.03160,
                0.16975,
                0.21192,
                0.15678,
                0.10852,
                0.08524,
                0.07724,
                0.05623,
                0.04032,
                0.02766,
                0.01582,
                0.00917,
                0.00483,
                0.00262,
                0.00099,
                0.00050,
                0.00027,
                0.00022,
                0.00011,
                0.00006,
                0.00005,
                0.00002,
                0.00001,
                0.00001,
                0.00001,
                0.00001,
                0.00001
            }
            ;
        double pdf_sum = 0.99997;

        switch (x)
        {
            case < 1:
                cdf = 0.0;
                break;
            default:
            {
                if (x < WORD_LENGTH_MAX)
                {
                    cdf = 0.0;
                    for (int i = 0; i < x; i++)
                    {
                        cdf += pdf_vec[i];
                    }

                    cdf /= pdf_sum;
                }
                else if (WORD_LENGTH_MAX <= x)
                {
                    cdf = 1.0;
                }

                break;
            }
        }

        return cdf;
    }

    public static int english_word_length_cdf_inv(double cdf)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ENGLISH_WORD_LENGTH_CDF_INV inverts the English Word Length CDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 July 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Henry Kucera, Winthrop Francis,
        //    Computational Analysis of Present-Day American English,
        //    Brown University Press, 1967.
        //
        //  Parameters:
        //
        //    Input, double CDF, the value of the CDF.
        //    0.0 <= CDF <= 1.0.
        //
        //    Output, int ENGLISH_WORD_LENGTH_CDF_INV, the corresponding word
        //    length for which CDF(X-1) < CDF <= CDF(X)
        //
    {
        int WORD_LENGTH_MAX = 27;

        double[] pdf_vec =  {
                0.03160,
                0.16975,
                0.21192,
                0.15678,
                0.10852,
                0.08524,
                0.07724,
                0.05623,
                0.04032,
                0.02766,
                0.01582,
                0.00917,
                0.00483,
                0.00262,
                0.00099,
                0.00050,
                0.00027,
                0.00022,
                0.00011,
                0.00006,
                0.00005,
                0.00002,
                0.00001,
                0.00001,
                0.00001,
                0.00001,
                0.00001
            }
            ;
        double pdf_sum = 0.99997;
        int x;

        switch (cdf)
        {
            case < 0.0:
            case > 1.0:
                Console.WriteLine("");
                Console.WriteLine("ENGLISH_WORD_LENGTH_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
        }

        double cum = 0.0;

        for (int j = 0; j < WORD_LENGTH_MAX; j++)
        {
            cum += pdf_vec[j];

            if (cdf <= cum / pdf_sum)
            {
                x = j + 1;
                return x;
            }
        }

        x = WORD_LENGTH_MAX;

        return x;
    }

    public static double english_word_length_mean()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ENGLISH_WORD_LENGTH_MEAN evaluates the mean of the English Word Length PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 July 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Henry Kucera, Winthrop Francis,
        //    Computational Analysis of Present-Day American English,
        //    Brown University Press, 1967.
        //
        //  Parameters:
        //
        //    Output, double ENGLISH_WORD_LENGTH_MEAN, the mean of the PDF.
        //
    {
        int WORD_LENGTH_MAX = 27;

        double[] pdf_vec =  {
                0.03160,
                0.16975,
                0.21192,
                0.15678,
                0.10852,
                0.08524,
                0.07724,
                0.05623,
                0.04032,
                0.02766,
                0.01582,
                0.00917,
                0.00483,
                0.00262,
                0.00099,
                0.00050,
                0.00027,
                0.00022,
                0.00011,
                0.00006,
                0.00005,
                0.00002,
                0.00001,
                0.00001,
                0.00001,
                0.00001,
                0.00001
            }
            ;
        double pdf_sum = 0.99997;

        double mean = 0.0;
        for (int j = 1; j <= WORD_LENGTH_MAX; j++)
        {
            mean += j * pdf_vec[j - 1];
        }

        mean /= pdf_sum;

        return mean;
    }

    public static double english_word_length_pdf(int x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ENGLISH_WORD_LENGTH_PDF evaluates the English Word Length PDF.
        //
        //  Discussion:
        //
        //    PDF(A,B;X) = B(X) if 1 <= X <= A
        //                = 0    otherwise
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 July 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Henry Kucera, Winthrop Francis,
        //    Computational Analysis of Present-Day American English,
        //    Brown University Press, 1967.
        //
        //  Parameters:
        //
        //    Input, int X, the word length whose probability is desired.
        //
        //    Output, double ENGLISH_WORD_LENGTH_PDF, the value of the PDF.
        //
    {
        int WORD_LENGTH_MAX = 27;

        double pdf;
        double[] pdf_vec =  {
                0.03160,
                0.16975,
                0.21192,
                0.15678,
                0.10852,
                0.08524,
                0.07724,
                0.05623,
                0.04032,
                0.02766,
                0.01582,
                0.00917,
                0.00483,
                0.00262,
                0.00099,
                0.00050,
                0.00027,
                0.00022,
                0.00011,
                0.00006,
                0.00005,
                0.00002,
                0.00001,
                0.00001,
                0.00001,
                0.00001,
                0.00001
            }
            ;
        double pdf_sum = 0.99997;

        pdf = x switch
        {
            >= 1 when x <= WORD_LENGTH_MAX => pdf_vec[x - 1] / pdf_sum,
            _ => 0.0
        };

        return pdf;
    }

    public static int english_word_length_sample(ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ENGLISH_WORD_LENGTH_SAMPLE samples the English Word Length PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 July 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Henry Kucera, Winthrop Francis,
        //    Computational Analysis of Present-Day American English,
        //    Brown University Press, 1967.
        //
        //  Parameters:
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, int ENGLISH_WORD_LENGTH_SAMPLE, a sample of the PDF.
        //
    {
        double cdf = UniformRNG.r8_uniform_01(ref seed);

        int x = english_word_length_cdf_inv(cdf);

        return x;
    }

    public static double english_word_length_variance()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ENGLISH_WORD_LENGTH_VARIANCE: variance of the English Word Length PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 July 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Henry Kucera, Winthrop Francis,
        //    Computational Analysis of Present-Day American English,
        //    Brown University Press, 1967.
        //
        //  Parameters:
        //
        //    Output, double ENGLISH_WORD_LENGTH_VARIANCE, the variance of the PDF.
        //
    {
        int WORD_LENGTH_MAX = 27;

        double[] pdf_vec =  {
                0.03160,
                0.16975,
                0.21192,
                0.15678,
                0.10852,
                0.08524,
                0.07724,
                0.05623,
                0.04032,
                0.02766,
                0.01582,
                0.00917,
                0.00483,
                0.00262,
                0.00099,
                0.00050,
                0.00027,
                0.00022,
                0.00011,
                0.00006,
                0.00005,
                0.00002,
                0.00001,
                0.00001,
                0.00001,
                0.00001,
                0.00001
            }
            ;
        double pdf_sum = 0.99997;
        double variance;

        double mean = 0.0;
        for (int j = 1; j <= WORD_LENGTH_MAX; j++)
        {
            mean += j * pdf_vec[j - 1];
        }

        mean /= pdf_sum;

        variance = 0.0;
        for (int j = 1; j <= WORD_LENGTH_MAX; j++)
        {
            variance += pdf_vec[j - 1] * Math.Pow(j - mean, 2);
        }

        variance /= pdf_sum;

        return variance;
    }
}