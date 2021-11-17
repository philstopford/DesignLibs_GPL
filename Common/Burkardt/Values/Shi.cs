﻿namespace Burkardt.Values;

public static class Shi
{
    public static void shi_values(ref int n_data, ref double x, ref double fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SHI_VALUES returns some values of the hyperbolic sine integral function.
        //
        //  Discussion:
        //
        //    SHI(X) = integral ( 0 <= T <= X ) sinh ( T ) / T dt
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      SinhIntegral[x]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 June 2007
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
        //    Output, ref double X, the argument of the function.
        //
        //    Output, ref double FX, the value of the function.
        //
    {
        const int N_MAX = 16;

        double[] fx_vec =
        {
            0.5069967498196672,
            0.6121303965633808,
            0.7193380189288998,
            0.8289965633789345,
            0.9414978265114335,
            1.057250875375729,
            1.300250361022057,
            1.561713388361002,
            1.845814141358504,
            2.157290343425901,
            2.501567433354976,
            3.549340406224435,
            4.973440475859807,
            6.966162067504942,
            9.817326911233034,
            13.96788504934715
        };

        double[] x_vec =
        {
            0.5E+00,
            0.6E+00,
            0.7E+00,
            0.8E+00,
            0.9E+00,
            1.0E+00,
            1.2E+00,
            1.4E+00,
            1.6E+00,
            1.8E+00,
            2.0E+00,
            2.5E+00,
            3.0E+00,
            3.5E+00,
            4.0E+00,
            4.5E+00
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