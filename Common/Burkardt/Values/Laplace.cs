﻿namespace Burkardt.Values;

public static class Laplace
{

    public static void laplace_cdf_values(ref int n_data, ref double mu, ref double beta, ref double x,
            ref double fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAPLACE_CDF_VALUES returns some values of the Laplace CDF.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Needs["Statistics`ContinuousDistributions`"]
        //      dist = LaplaceDistribution [ mu, beta ]
        //      CDF [ dist, x ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 August 2004
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
        //    Output, ref double MU, the mean of the distribution.
        //
        //    Output, ref double BETA, the shape parameter.
        //
        //    Output, ref double X, the argument of the function.
        //
        //    Output, ref double FX, the value of the function.
        //
    {
        const int N_MAX = 12;

        double[] beta_vec =
        {
            0.1000000000000000E+01,
            0.1000000000000000E+01,
            0.1000000000000000E+01,
            0.1000000000000000E+01,
            0.2000000000000000E+01,
            0.3000000000000000E+01,
            0.4000000000000000E+01,
            0.5000000000000000E+01,
            0.2000000000000000E+01,
            0.2000000000000000E+01,
            0.2000000000000000E+01,
            0.2000000000000000E+01
        };

        double[] fx_vec =
        {
            0.5000000000000000E+00,
            0.8160602794142788E+00,
            0.9323323583816937E+00,
            0.9751064658160680E+00,
            0.6967346701436833E+00,
            0.6417343447131054E+00,
            0.6105996084642976E+00,
            0.5906346234610091E+00,
            0.5000000000000000E+00,
            0.3032653298563167E+00,
            0.1839397205857212E+00,
            0.1115650800742149E+00
        };

        double[] mu_vec =
        {
            0.0000000000000000E+01,
            0.0000000000000000E+01,
            0.0000000000000000E+01,
            0.0000000000000000E+01,
            0.0000000000000000E+01,
            0.0000000000000000E+01,
            0.0000000000000000E+01,
            0.0000000000000000E+01,
            0.1000000000000000E+01,
            0.2000000000000000E+01,
            0.3000000000000000E+01,
            0.4000000000000000E+01
        };

        double[] x_vec =
        {
            0.0000000000000000E+01,
            0.1000000000000000E+01,
            0.2000000000000000E+01,
            0.3000000000000000E+01,
            0.1000000000000000E+01,
            0.1000000000000000E+01,
            0.1000000000000000E+01,
            0.1000000000000000E+01,
            0.1000000000000000E+01,
            0.1000000000000000E+01,
            0.1000000000000000E+01,
            0.1000000000000000E+01
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
            mu = 0.0;
            beta = 0.0;
            x = 0.0;
            fx = 0.0;
        }
        else
        {
            mu = mu_vec[n_data - 1];
            beta = beta_vec[n_data - 1];
            x = x_vec[n_data - 1];
            fx = fx_vec[n_data - 1];
        }
    }

}