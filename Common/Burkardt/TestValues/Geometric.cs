namespace Burkardt.TestValues
{
    public static class Geometric
    {
        public static void geometric_cdf_values(ref int n_data, ref int x, ref double p, ref double cdf)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GEOMETRIC_CDF_VALUES returns values of the geometric CDF.
            //
            //  Discussion:
            //
            //    The geometric or Pascal probability density function gives the
            //    probability that the first success will happen on the X-th Bernoulli
            //    trial, given that the probability of a success on a single trial is P.
            //
            //    The value of CDF ( X, P ) is the probability that the first success
            //    will happen on or before the X-th trial.
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Needs["Statistics`DiscreteDistributions`]
            //      dist = GeometricDistribution [ p ]
            //      CDF [ dist, x ]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    21 August 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Stephen Wolfram,
            //    The Mathematica Book,
            //    Fourth Edition,
            //    Cambridge University Press, 1999,
            //    ISBN: 0-521-64314-7,
            //    LC: QA76.95.W65.
            //
            //    Daniel Zwillinger and Stephen Kokoska,
            //    CRC Standard Probability and Statistics Tables and Formulae,
            //    Chapman and Hall / CRC Press, 2000.
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref int X, the number of trials.
            //
            //    Output, ref double P, the probability of success
            //    on one trial.
            //
            //    Output, ref double CDF, the cumulative density function.
            //
        {
            int N_MAX = 14;

            double[] cdf_vec =
            {
                0.1900000000000000E+00,
                0.2710000000000000E+00,
                0.3439000000000000E+00,
                0.6861894039100000E+00,
                0.3600000000000000E+00,
                0.4880000000000000E+00,
                0.5904000000000000E+00,
                0.9141006540800000E+00,
                0.7599000000000000E+00,
                0.8704000000000000E+00,
                0.9375000000000000E+00,
                0.9843750000000000E+00,
                0.9995117187500000E+00,
                0.9999000000000000E+00
            };

            double[] p_vec =
            {
                0.1E+00,
                0.1E+00,
                0.1E+00,
                0.1E+00,
                0.2E+00,
                0.2E+00,
                0.2E+00,
                0.2E+00,
                0.3E+00,
                0.4E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.9E+00
            };

            int[] x_vec =
            {
                1, 2, 3, 10, 1,
                2, 3, 10, 3, 3,
                3, 5, 10, 3
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                x = 0;
                p = 0.0;
                cdf = 0.0;
            }
            else
            {
                x = x_vec[n_data - 1];
                p = p_vec[n_data - 1];
                cdf = cdf_vec[n_data - 1];
            }
        }
    }
}