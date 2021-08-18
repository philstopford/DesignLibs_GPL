using System;

namespace Burkardt.Values
{
    public static class Binomial
    {
        public static int choose ( int n, int k )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CHOOSE computes the binomial coefficient C(N,K).
            //
            //  Discussion:
            //
            //    The value is calculated in such a way as to avoid overflow and
            //    roundoff.  The calculation is done in integer arithmetic.
            //
            //    The formula used is:
            //
            //      C(N,K) = N! / ( K! * (N-K)! )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    22 May 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    ML Wolfson, HV Wright,
            //    Algorithm 160:
            //    Combinatorial of M Things Taken N at a Time,
            //    Communications of the ACM,
            //    Volume 6, Number 4, April 1963, page 161.
            //
            //  Parameters:
            //
            //    Input, int N, K, are the values of N and K.
            //
            //    Output, int CHOOSE, the number of combinations of N
            //    things taken K at a time.
            //
        {
            int i;
            int mn;
            int mx;
            int value;

            mn = Math.Min ( k, n - k );

            if ( mn < 0 )
            {
                value = 0;
            }
            else if ( mn == 0 )
            {
                value = 1;
            }
            else
            {
                mx = Math.Max ( k, n - k );
                value = mx + 1;

                for ( i = 2; i <= mn; i++ )
                {
                    value = ( value * ( mx + i ) ) / i;
                }
            }

            return value;
        }
        public static void binomial_values(ref int n_data, ref int a, ref int b, ref int fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BINOMIAL_VALUES returns some values of the binomial coefficients.
            //
            //  Discussion:
            //
            //    The formula for the binomial coefficient is
            //
            //      C(N,K) = N! / ( K! * (N-K)! )
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Binomial[n,k]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 August 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Milton Abramowitz, Irene Stegun,
            //    Handbook of Mathematical Functions,
            //    National Bureau of Standards, 1964,
            //    ISBN: 0-486-61272-4,
            //    LC: QA47.A34.
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
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref int A, &B, the arguments of the function.
            //
            //    Output, ref int FX, the value of the function.
            //
        {
            int N_MAX = 20;

            int[] a_vec =
            {
                1, 6, 6, 6, 15,
                15, 15, 15, 15, 15,
                15, 25, 25, 25, 25,
                25, 25, 25, 25, 25
            };

            int[] b_vec =
            {
                0, 1, 3, 5, 1,
                3, 5, 7, 9, 11,
                13, 1, 3, 5, 7,
                9, 11, 13, 15, 17
            };

            int[] fx_vec =
            {
                1,
                6,
                20,
                6,
                15,
                455,
                3003,
                6435,
                5005,
                1365,
                105,
                25,
                2300,
                53130,
                480700,
                2042975,
                4457400,
                5200300,
                3268760,
                1081575
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                a = 0;
                b = 0;
                fx = 0;
            }
            else
            {
                a = a_vec[n_data - 1];
                b = b_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void binomial_cdf_values(ref int n_data, ref int a, ref double b, ref int x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BINOMIAL_CDF_VALUES returns some values of the binomial CDF.
            //
            //  Discussion:
            //
            //    CDF(X)(A,B) is the probability of at most X successes in A trials,
            //    given that the probability of success on a single trial is B.
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Needs["Statistics`DiscreteDistributions`]
            //      dist = BinomialDistribution [ n, p ]
            //      CDF [ dist, x ]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    15 August 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Milton Abramowitz, Irene Stegun,
            //    Handbook of Mathematical Functions,
            //    National Bureau of Standards, 1964,
            //    ISBN: 0-486-61272-4,
            //    LC: QA47.A34.
            //
            //    Stephen Wolfram,
            //    The Mathematica Book,
            //    Fourth Edition,
            //    Cambridge University Press, 1999,
            //    ISBN: 0-521-64314-7,
            //    LC: QA76.95.W65.
            //
            //    Daniel Zwillinger,
            //    CRC Standard Mathematical Tables and Formulae,
            //    30th Edition, CRC Press, 1996, pages 651-652.
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref int A, a parameter of the function.
            //
            //    Output, ref double B, a parameter of the function.
            //
            //    Output, ref int X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 17;

            int[] a_vec =
            {
                2, 2, 2, 2,
                2, 4, 4, 4,
                4, 10, 10, 10,
                10, 10, 10, 10,
                10
            };

            double[] b_vec =
            {
                0.05E+00,
                0.05E+00,
                0.05E+00,
                0.50E+00,
                0.50E+00,
                0.25E+00,
                0.25E+00,
                0.25E+00,
                0.25E+00,
                0.05E+00,
                0.10E+00,
                0.15E+00,
                0.20E+00,
                0.25E+00,
                0.30E+00,
                0.40E+00,
                0.50E+00
            };

            double[] fx_vec =
            {
                0.9025000000000000E+00,
                0.9975000000000000E+00,
                0.1000000000000000E+01,
                0.2500000000000000E+00,
                0.7500000000000000E+00,
                0.3164062500000000E+00,
                0.7382812500000000E+00,
                0.9492187500000000E+00,
                0.9960937500000000E+00,
                0.9999363101685547E+00,
                0.9983650626000000E+00,
                0.9901259090013672E+00,
                0.9672065024000000E+00,
                0.9218730926513672E+00,
                0.8497316674000000E+00,
                0.6331032576000000E+00,
                0.3769531250000000E+00
            };

            int[] x_vec =
            {
                0, 1, 2, 0,
                1, 0, 1, 2,
                3, 4, 4, 4,
                4, 4, 4, 4,
                4
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                a = 0;
                b = 0.0;
                x = 0;
                fx = 0.0;
            }
            else
            {
                a = a_vec[n_data - 1];
                b = b_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void binomial_pdf_values(ref int n_data, ref int a, ref double b, ref int x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BINOMIAL_PDF_VALUES returns some values of the binomial PDF.
            //
            //  Discussion:
            //
            //    PDF(X)(A,B) is the probability of X successes in A trials,
            //    given that the probability of success on a single trial is B.
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Needs["Statistics`DiscreteDistributions`]
            //      dist = BinomialDistribution [ n, p ]
            //      PDF [ dist, x ]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    27 July 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Milton Abramowitz, Irene Stegun,
            //    Handbook of Mathematical Functions,
            //    National Bureau of Standards, 1964,
            //    ISBN: 0-486-61272-4,
            //    LC: QA47.A34.
            //
            //    Stephen Wolfram,
            //    The Mathematica Book,
            //    Fourth Edition,
            //    Cambridge University Press, 1999,
            //    ISBN: 0-521-64314-7,
            //    LC: QA76.95.W65.
            //
            //    Daniel Zwillinger,
            //    CRC Standard Mathematical Tables and Formulae,
            //    30th Edition, CRC Press, 1996, pages 651-652.
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref int A, a parameter of the function.
            //
            //    Output, ref double B, a parameter of the function.
            //
            //    Output, ref int X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 10;

            int[] a_vec =
            {
                5, 12, 6, 13, 9,
                1, 2, 17, 6, 8
            };

            double[] b_vec =
            {
                0.8295092339327006,
                0.06611873491603133,
                0.0438289977791071,
                0.4495389603071763,
                0.7972869541062562,
                0.3507523379805466,
                0.8590968552798568,
                0.007512364073964213,
                0.1136640464424993,
                0.2671322702601793
            };

            double[] fx_vec =
            {
                0.3927408939646697,
                0.0006199968732461383,
                0.764211224733124,
                0.0004260353334364943,
                0.302948289145794,
                0.3507523379805466,
                0.01985369619202562,
                0.006854388879646552,
                0.000002156446446382985,
                0.0005691150511772053
            };

            int[] x_vec =
            {
                5, 5, 0, 0, 7,
                1, 0, 2, 6, 7
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                a = 0;
                b = 0.0;
                x = 0;
                fx = 0.0;
            }
            else
            {
                a = a_vec[n_data - 1];
                b = b_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

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
            //    Stephen Wolfram,
            //    The Mathematica Book,
            //    Fourth Edition,
            //    Cambridge University Press, 1999,
            //    ISBN: 0-521-64314-7,
            //    LC: QA76.95.W65.
            //
            //  Input:
            //
            //    ref int N_DATA.  The user sets N_DATA to 0 before the first call.
            //
            //  Output:
            //
            //    ref int N_DATA.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    ref int F, the maximum number of failures.
            //
            //    ref int S, the number of successes.
            //
            //    ref double P, the probability of a success on one trial.
            //
            //    ref double CDF, the probability of at most F failures
            //    before the S-th success.
            //
        {
            int N_MAX = 27;

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

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

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

    }
}