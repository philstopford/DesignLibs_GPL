using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.Probability;

public static class NegativeBinomial
{
    public static double negative_binomial_cdf(int x, int a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NEGATIVE_BINOMIAL_CDF evaluates the Negative Binomial CDF.
        //
        //  Discussion:
        //
        //    A simple summing approach is used.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int X, the argument of the CDF.
        //
        //    Input, int A, double B, parameters of the PDF.
        //    0 <= A,
        //    0 < B <= 1.
        //
        //    Output, double NEGATIVE_BINOMIAL_CDF, the value of the CDF.
        //
    {
        double cdf = 0.0;

        for (int y = a; y <= x; y++)
        {
            int cnk = typeMethods.i4_choose(y - 1, a - 1);

            double pdf = cnk * Math.Pow(b, a) * Math.Pow(1.0 - b, y - a);

            cdf += pdf;
        }

        return cdf;
    }

    public static int negative_binomial_cdf_inv(double cdf, int a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NEGATIVE_BINOMIAL_CDF_INV inverts the Negative Binomial CDF.
        //
        //  Discussion:
        //
        //    A simple discrete approach is used.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double CDF, the value of the CDF.
        //
        //    Input, int A, double B, parameters of the PDF.
        //    0 <= A,
        //    0 < B <= 1.
        //
        //    Output, int NEGATIVE_BINOMIAL_CDF_INV, the smallest X whose cumulative
        //    density function is greater than or equal to CDF.
        //
    {
        double cum;
        double pdf;
        int x;
        int x_max = 1000;
        switch (cdf)
        {
    
//
            case < 0.0:
            case > 1.0:
                Console.WriteLine(" ");
                Console.WriteLine("NEGATIVE_BINOMIAL_CDF_INV - Fatal error!");
                Console.WriteLine("  CDF < 0 or 1 < CDF.");
                return 1;
        }

        cum = 0.0;

        x = a;

        for (;;)
        {
            pdf = negative_binomial_pdf(x, a, b);

            cum += pdf;

            if (cdf <= cum || x_max <= x)
            {
                break;
            }

            x += 1;
        }

        return x;
    }

    public static void negative_binomial_cdf_values(ref int n_data, ref int f, ref int s, ref double p,
            ref double cdf )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NEGATIVE_BINOMIAL_CDF_VALUES returns values of the negative binomial CDF.
        //
        //  Discussion:
        //
        //    Assume that a coin has a probability P of coming up heads on
        //    any one trial.  Suppose that we plan to flip the coin until we
        //    achieve a total of S heads.  If we let F represent the number of
        //    tails that occur in this process, then the value of F satisfies
        //    a negative binomial PDF:
        //
        //      PDF(F,S,P) = Choose ( F from F+S-1 ) * P^S * (1-P)^F
        //
        //    The negative binomial CDF is the probability that there are F or
        //    fewer failures upon the attainment of the S-th success.  Thus,
        //
        //      CDF(F,S,P) = sum ( 0 <= G <= F ) PDF(G,S,P)
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Needs["Statistics`DiscreteDistributions`]
        //      dist = NegativeBinomialDistribution [ s, p ]
        //      CDF [ dist, f ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 August 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    F C Powell,
        //    Statistical Tables for Sociology, Biology and Physical Sciences,
        //    Cambridge University Press, 1982.
        //
        //    Stephen Wolfram,
        //    The Mathematica Book,
        //    Fourth Edition,
        //    Cambridge University Press, 1999,
        //    ISBN: 0-521-64314-7,
        //    LC: QA76.95.W65.
        //
        //  Parameters:
        //
        //    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, int &F, the maximum number of failures.
        //
        //    Output, int &S, the number of successes.
        //
        //    Output, double &P, the probability of a success on one trial.
        //
        //    Output, double &CDF, the probability of at most F failures
        //    before the S-th success.
        //
    {
        const int N_MAX = 27;

        double[] cdf_vec =
        {
            0.6367187500000000E+00,
            0.3632812500000000E+00,
            0.1445312500000000E+00,
            0.5000000000000000E+00,
            0.2265625000000000E+00,
            0.6250000000000000E-01,
            0.3437500000000000E+00,
            0.1093750000000000E+00,
            0.1562500000000000E-01,
            0.1792000000000000E+00,
            0.4096000000000000E-01,
            0.4096000000000000E-02,
            0.7047000000000000E-01,
            0.1093500000000000E-01,
            0.7290000000000000E-03,
            0.9861587127990000E+00,
            0.9149749500510000E+00,
            0.7471846521450000E+00,
            0.8499053647030009E+00,
            0.5497160941090026E+00,
            0.2662040052146710E+00,
            0.6513215599000000E+00,
            0.2639010709000000E+00,
            0.7019082640000000E-01,
            0.1000000000000000E+01,
            0.1990000000000000E-01,
            0.1000000000000000E-03
        };

        int[] f_vec =
        {
            4, 3, 2,
            3, 2, 1,
            2, 1, 0,
            2, 1, 0,
            2, 1, 0,
            11, 10, 9,
            17, 16, 15,
            9, 8, 7,
            2, 1, 0
        };

        double[] p_vec =
        {
            0.50E+00,
            0.50E+00,
            0.50E+00,
            0.50E+00,
            0.50E+00,
            0.50E+00,
            0.50E+00,
            0.50E+00,
            0.50E+00,
            0.40E+00,
            0.40E+00,
            0.40E+00,
            0.30E+00,
            0.30E+00,
            0.30E+00,
            0.30E+00,
            0.30E+00,
            0.30E+00,
            0.10E+00,
            0.10E+00,
            0.10E+00,
            0.10E+00,
            0.10E+00,
            0.10E+00,
            0.10E-01,
            0.10E-01,
            0.10E-01
        };

        int[] s_vec =
        {
            4, 5, 6,
            4, 5, 6,
            4, 5, 6,
            4, 5, 6,
            4, 5, 6,
            1, 2, 3,
            1, 2, 3,
            1, 2, 3,
            0, 1, 2
        };

        n_data = n_data switch
        {
            < 0 => 0,
            _ => n_data
        };

        n_data += 1;

        if (N_MAX < n_data)
        {
            n_data = 0;
            f = 0;
            s = 0;
            p = 0.0;
            cdf = 0.0;
        }
        else
        {
            f = f_vec[n_data - 1];
            s = s_vec[n_data - 1];
            p = p_vec[n_data - 1];
            cdf = cdf_vec[n_data - 1];
        }
    }

    public static bool negative_binomial_check(int a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NEGATIVE_BINOMIAL_CHECK checks parameters of the Negative Binomial PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int A, double B, parameters of the PDF.
        //    0 <= A,
        //    0 < B <= 1.
        //
        //    Output, bool NEGATIVE_BINOMIAL_CHECK, is true if the
        //    parameters are legal.
        //
    {
        switch (a)
        {
            case < 0:
                Console.WriteLine(" ");
                Console.WriteLine("NEGATIVE_BINOMIAL_CHECK - Warning!");
                Console.WriteLine("  A < 0.");
                return false;
        }

        switch (b)
        {
            case <= 0.0:
            case > 1.0:
                Console.WriteLine(" ");
                Console.WriteLine("NEGATIVE_BINOMIAL_CHECK - Warning!");
                Console.WriteLine("  B <= 0 or 1 < B.");
                return false;
            default:
                return true;
        }
    }

    public static double negative_binomial_mean(int a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NEGATIVE_BINOMIAL_MEAN returns the mean of the Negative Binomial PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int A, double B, parameters of the PDF.
        //    0 <= A,
        //    0 < B <= 1.
        //
        //    Output, double NEGATIVE_BINOMIAL_MEAN, the mean of the PDF.
        //
    {
        double mean = a / b;

        return mean;
    }

    public static double negative_binomial_pdf(int x, int a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NEGATIVE_BINOMIAL_PDF evaluates the Negative Binomial PDF.
        //
        //  Discussion:
        //
        //    PDF(A,B;X) = C(X-1,A-1) * B^A * ( 1 - B )^(X-A)
        //
        //    PDF(A,B;X) is the probability that the A-th success will
        //    occur on the X-th trial, given that the probability
        //    of a success on a single trial is B.
        //
        //    The Negative Binomial PDF is also known as the Pascal PDF or
        //    the "Polya" PDF.
        //
        //    NEGATIVE_BINOMIAL_PDF(1,B;X) = GEOMETRIC_PDF(B;X)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int X, the number of trials.
        //    A <= X.
        //
        //    Input, int A, the number of successes required.
        //    0 <= A <= X, normally.
        //
        //    Input, double B, the probability of a success on a single trial.
        //    0.0 < B <= 1.0.
        //
        //    Output, double NEGATIVE_BINOMIAL_PDF, the value of the PDF.
        //
    {
        double pdf;

        if (x < a)
        {
            pdf = 0.0;
        }
        else
        {
            int cnk = typeMethods.i4_choose(x - 1, a - 1);

            pdf = cnk * Math.Pow(b, a) * Math.Pow(1.0 - b, x - a);
        }

        return pdf;
    }

    public static int negative_binomial_sample(int a, double b, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NEGATIVE_BINOMIAL_SAMPLE samples the Negative Binomial PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int A, double B, parameters of the PDF.
        //    0 <= A,
        //    0 < B <= 1.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, int NEGATIVE_BINOMIAL_SAMPLE, a sample of the PDF.
        //
    {
        double r;
        int x;

        switch (b)
        {
            case 1.0:
                x = a;
                return x;
            case 0.0:
                x = typeMethods.i4_huge();
                return x;
        }

        x = 0;
        int num_success = 0;

        while (num_success < a)
        {
            x += 1;
            r = UniformRNG.r8_uniform_01(ref seed);

            if (r <= b)
            {
                num_success += 1;
            }

        }

        return x;
    }

    public static double negative_binomial_variance(int a, double b)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NEGATIVE_BINOMIAL_VARIANCE returns the variance of the Negative Binomial PDF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int A, double B, parameters of the PDF.
        //    0 <= A,
        //    0 < B <= 1.
        //
        //    Output, double NEGATIVE_BINOMIAL_VARIANCE, the variance of the PDF.
        //
    {
        double variance = a * (1.0 - b) / (b * b);

        return variance;
    }
}