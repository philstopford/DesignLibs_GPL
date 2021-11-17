﻿namespace Burkardt.Values;

public static class Gegenbauer
{

    public static void gegenbauer_poly_values(ref int n_data, ref int n, ref double a, ref double x,
            ref double fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEGENBAUER_POLY_VALUES returns some values of the Gegenbauer polynomials.
        //
        //  Discussion:
        //
        //    The Gegenbauer polynomials are also known as the "spherical
        //    polynomials" or "ultraspherical polynomials".
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      GegenbauerC[n,m,x]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 August 2004
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
        //    Output, ref int N, the order parameter of the function.
        //
        //    Output, ref double A, the real parameter of the function.
        //
        //    Output, ref double X, the argument of the function.
        //
        //    Output, ref double FX, the value of the function.
        //
    {
        const int N_MAX = 38;

        double[] a_vec =
        {
            0.5E+00,
            0.5E+00,
            0.5E+00,
            0.5E+00,
            0.5E+00,
            0.5E+00,
            0.5E+00,
            0.5E+00,
            0.5E+00,
            0.5E+00,
            0.5E+00,
            0.0E+00,
            1.0E+00,
            2.0E+00,
            3.0E+00,
            4.0E+00,
            5.0E+00,
            6.0E+00,
            7.0E+00,
            8.0E+00,
            9.0E+00,
            10.0E+00,
            3.0E+00,
            3.0E+00,
            3.0E+00,
            3.0E+00,
            3.0E+00,
            3.0E+00,
            3.0E+00,
            3.0E+00,
            3.0E+00,
            3.0E+00,
            3.0E+00,
            3.0E+00,
            3.0E+00,
            3.0E+00,
            3.0E+00,
            3.0E+00
        };

        double[] fx_vec =
        {
            1.0000000000E+00,
            0.2000000000E+00,
            -0.4400000000E+00,
            -0.2800000000E+00,
            0.2320000000E+00,
            0.3075200000E+00,
            -0.0805760000E+00,
            -0.2935168000E+00,
            -0.0395648000E+00,
            0.2459712000E+00,
            0.1290720256E+00,
            0.0000000000E+00,
            -0.3600000000E+00,
            -0.0800000000E+00,
            0.8400000000E+00,
            2.4000000000E+00,
            4.6000000000E+00,
            7.4400000000E+00,
            10.9200000000E+00,
            15.0400000000E+00,
            19.8000000000E+00,
            25.2000000000E+00,
            -9.0000000000E+00,
            -0.1612800000E+00,
            -6.6729600000E+00,
            -8.3750400000E+00,
            -5.5267200000E+00,
            0.0000000000E+00,
            5.5267200000E+00,
            8.3750400000E+00,
            6.6729600000E+00,
            0.1612800000E+00,
            -9.0000000000E+00,
            -15.4252800000E+00,
            -9.6969600000E+00,
            22.4409600000E+00,
            100.8892800000E+00,
            252.0000000000E+00
        };

        int[] n_vec =
        {
            0, 1, 2,
            3, 4, 5,
            6, 7, 8,
            9, 10, 2,
            2, 2, 2,
            2, 2, 2,
            2, 2, 2,
            2, 5, 5,
            5, 5, 5,
            5, 5, 5,
            5, 5, 5,
            5, 5, 5,
            5, 5
        };

        double[] x_vec =
        {
            0.20E+00,
            0.20E+00,
            0.20E+00,
            0.20E+00,
            0.20E+00,
            0.20E+00,
            0.20E+00,
            0.20E+00,
            0.20E+00,
            0.20E+00,
            0.20E+00,
            0.40E+00,
            0.40E+00,
            0.40E+00,
            0.40E+00,
            0.40E+00,
            0.40E+00,
            0.40E+00,
            0.40E+00,
            0.40E+00,
            0.40E+00,
            0.40E+00,
            -0.50E+00,
            -0.40E+00,
            -0.30E+00,
            -0.20E+00,
            -0.10E+00,
            0.00E+00,
            0.10E+00,
            0.20E+00,
            0.30E+00,
            0.40E+00,
            0.50E+00,
            0.60E+00,
            0.70E+00,
            0.80E+00,
            0.90E+00,
            1.00E+00
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
            n = 0;
            a = 0.0;
            x = 0.0;
            fx = 0.0;
        }
        else
        {
            n = n_vec[n_data - 1];
            a = a_vec[n_data - 1];
            x = x_vec[n_data - 1];
            fx = fx_vec[n_data - 1];
        }
    }

}