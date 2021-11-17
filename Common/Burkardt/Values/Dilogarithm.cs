﻿namespace Burkardt.Values;

public static class Dilogarithm
{
    public static void dilogarithm_values(ref int n_data, ref double x, ref double fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DILOGARITHM_VALUES returns some values of the dilogarithm function.
        //
        //  Discussion:
        //
        //    The dilogarithm is defined as
        //
        //      Li_2(X) = - integral ( 0 <= T <= X ) ln ( 1 - T ) / T dT
        //
        //    The dilogarithm is also known as Spence's integral.
        //
        //    In Abramowitz and Stegun form of the function is different,
        //    and is equivalent to evaluated Li_2(1-X).
        //
        //    The dilogarithm is the special case, with N = 2, of the
        //    polylogarithm Li_N(X).
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      PolyLog[2,X]
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
        //    Wolfram Media / Cambridge University Press, 1999.
        //
        //  Parameters:
        //
        //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, ref double X, the argument of the function.
        //
        //    Output, ref double FX, the value of the function.
        //
    {
        const int N_MAX = 21;

        double[] fx_vec =
        {
            0.0000000000000000E+00,
            0.5063929246449603E-01,
            0.1026177910993911E+00,
            0.1560350339454831E+00,
            0.2110037754397048E+00,
            0.2676526390827326E+00,
            0.3261295100754761E+00,
            0.3866059411605865E+00,
            0.4492829744712817E+00,
            0.5143989891542119E+00,
            0.5822405264650125E+00,
            0.6531576315069018E+00,
            0.7275863077163334E+00,
            0.8060826895177240E+00,
            0.8893776242860387E+00,
            0.9784693929303061E+00,
            0.1074794600008248E+01,
            0.1180581123830255E+01,
            0.1299714723004959E+01,
            0.1440633796970039E+01,
            0.1644934066848226E+01
        };

        double[] x_vec =
        {
            0.00E+00,
            0.05E+00,
            0.10E+00,
            0.15E+00,
            0.20E+00,
            0.25E+00,
            0.30E+00,
            0.35E+00,
            0.40E+00,
            0.45E+00,
            0.50E+00,
            0.55E+00,
            0.60E+00,
            0.65E+00,
            0.70E+00,
            0.75E+00,
            0.80E+00,
            0.85E+00,
            0.90E+00,
            0.95E+00,
            0.10E+01
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
            x = 0.0;
            fx = 0.0;
        }
        else
        {
            x = x_vec[n_data - 1];
            fx = fx_vec[n_data - 1];
        }
    }
}