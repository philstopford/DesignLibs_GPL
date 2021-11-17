﻿namespace Burkardt.Values;

public static class Clausen
{
    public static void clausen_values(ref int n_data, ref double x, ref double fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CLAUSEN_VALUES returns some values of the Clausen's integral.
        //
        //  Discussion:
        //
        //    The function is defined by:
        //
        //      CLAUSEN(x) = integral ( 0 <= t <= x ) -ln ( 2 * sin ( t / 2 ) ) dt
        //
        //    The data was reported by McLeod.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 August 2004
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
        //    Allan McLeod,
        //    Algorithm 757:
        //    MISCFUN: A software package to compute uncommon special functions,
        //    ACM Transactions on Mathematical Software,
        //    Volume 22, Number 3, September 1996, pages 288-301.
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
        const int N_MAX = 20;

        double[] fx_vec =
        {
            0.14137352886760576684E-01,
            0.13955467081981281934E+00,
            -0.38495732156574238507E+00,
            0.84831187770367927099E+00,
            0.10139591323607685043E+01,
            -0.93921859275409211003E+00,
            0.72714605086327924743E+00,
            0.43359820323553277936E+00,
            -0.98026209391301421161E-01,
            -0.56814394442986978080E+00,
            -0.70969701784448921625E+00,
            0.99282013254695671871E+00,
            -0.98127747477447367875E+00,
            -0.64078266570172320959E+00,
            0.86027963733231192456E+00,
            0.39071647608680211043E+00,
            0.47574793926539191502E+00,
            0.10105014481412878253E+01,
            0.96332089044363075154E+00,
            -0.61782699481929311757E+00
        };

        double[] x_vec =
        {
            0.0019531250E+00,
            0.0312500000E+00,
            -0.1250000000E+00,
            0.5000000000E+00,
            1.0000000000E+00,
            -1.5000000000E+00,
            2.0000000000E+00,
            2.5000000000E+00,
            -3.0000000000E+00,
            4.0000000000E+00,
            4.2500000000E+00,
            -5.0000000000E+00,
            5.5000000000E+00,
            6.0000000000E+00,
            8.0000000000E+00,
            -10.0000000000E+00,
            15.0000000000E+00,
            20.0000000000E+00,
            -30.0000000000E+00,
            50.0000000000E+00
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