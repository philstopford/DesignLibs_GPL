namespace Burkardt.CDFLib;

public static partial class CDF
{
    public static void negative_binomial_cdf_values(ref int n_data, ref int f, ref int s, ref double p,
            ref double cdf)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    negative_binomial_cdf_values() returns values of the negative binomial CDF.
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
        //    Thanks to Sebastiono Vigna for reporting an error in which the
        //    expression "if ( n_data < 0 )" was used, 05 May 2021.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 May 2021
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
        //  Input:
        //
        //    int *N_DATA.  The user sets N_DATA to 0 before the first call.
        //
        //  Output:
        //
        //    int *N_DATA.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    int *F, the maximum number of failures.
        //
        //    int *S, the number of successes.
        //
        //    double *P, the probability of a success on one trial.
        //
        //    double *CDF, the probability of at most F failures before the
        //    S-th success.
        //
    {
        const int N_MAX = 27;

        double[] cdf_vec =
            {
                0.6367, 0.3633, 0.1445,
                0.5000, 0.2266, 0.0625,
                0.3438, 0.1094, 0.0156,
                0.1792, 0.0410, 0.0041,
                0.0705, 0.0109, 0.0007,
                0.9862, 0.9150, 0.7472,
                0.8499, 0.5497, 0.2662,
                0.6513, 0.2639, 0.0702,
                1.0000, 0.0199, 0.0001
            }
            ;
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
            }
            ;
        double[] p_vec =
            {
                0.50, 0.50, 0.50,
                0.50, 0.50, 0.50,
                0.50, 0.50, 0.50,
                0.40, 0.40, 0.40,
                0.30, 0.30, 0.30,
                0.30, 0.30, 0.30,
                0.10, 0.10, 0.10,
                0.10, 0.10, 0.10,
                0.01, 0.01, 0.01
            }
            ;
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
            }
            ;

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
            p = 0.0E+00;
            cdf = 0.0E+00;
        }
        else
        {
            f = f_vec[n_data - 1];
            s = s_vec[n_data - 1];
            p = p_vec[n_data - 1];
            cdf = cdf_vec[n_data - 1];
        }
    }
}