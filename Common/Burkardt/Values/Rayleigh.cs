namespace Burkardt.Values;

public static class Rayleigh
{

    public static void rayleigh_cdf_values(ref int n_data, ref double sigma, ref double x, ref double fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RAYLEIGH_CDF_VALUES returns some values of the Rayleigh CDF.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Needs["Statistics`ContinuousDistributions`"]
        //      dist = RayleighDistribution [ sigma ]
        //      CDF [ dist, x ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 September 2004
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
        //    Output, ref double SIGMA, the shape parameter of the distribution.
        //
        //    Output, ref double X, the argument of the function.
        //
        //    Output, ref double FX, the value of the function.
        //
    {
        const int N_MAX = 9;

        double[] fx_vec =
        {
            0.8646647167633873E+00,
            0.9996645373720975E+00,
            0.9999999847700203E+00,
            0.999999999999987E+00,
            0.8646647167633873E+00,
            0.3934693402873666E+00,
            0.1992625970831920E+00,
            0.1175030974154046E+00,
            0.7688365361336422E-01
        };

        double[] sigma_vec =
        {
            0.5000000000000000E+00,
            0.5000000000000000E+00,
            0.5000000000000000E+00,
            0.5000000000000000E+00,
            0.1000000000000000E+01,
            0.2000000000000000E+01,
            0.3000000000000000E+01,
            0.4000000000000000E+01,
            0.5000000000000000E+01
        };

        double[] x_vec =
        {
            0.1000000000000000E+01,
            0.2000000000000000E+01,
            0.3000000000000000E+01,
            0.4000000000000000E+01,
            0.2000000000000000E+01,
            0.2000000000000000E+01,
            0.2000000000000000E+01,
            0.2000000000000000E+01,
            0.2000000000000000E+01
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
            sigma = 0.0;
            x = 0.0;
            fx = 0.0;
        }
        else
        {
            sigma = sigma_vec[n_data - 1];
            x = x_vec[n_data - 1];
            fx = fx_vec[n_data - 1];
        }
    }

}