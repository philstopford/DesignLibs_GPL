namespace Burkardt.Values;

public static class ErrorFunc
{
    public static void erf_values(ref int n_data, ref double x, ref double fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ERF_VALUES returns some values of the ERF or "error" function.
        //
        //  Discussion:
        //
        //    The error function is defined by:
        //
        //      ERF(X) = ( 2 / sqrt ( PI ) * integral ( 0 <= T <= X ) exp ( - T^2 ) dT
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Erf[x]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 August 2004
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
        const int N_MAX = 21;

        double[] fx_vec =
        {
            0.0000000000000000E+00,
            0.1124629160182849E+00,
            0.2227025892104785E+00,
            0.3286267594591274E+00,
            0.4283923550466685E+00,
            0.5204998778130465E+00,
            0.6038560908479259E+00,
            0.6778011938374185E+00,
            0.7421009647076605E+00,
            0.7969082124228321E+00,
            0.8427007929497149E+00,
            0.8802050695740817E+00,
            0.9103139782296354E+00,
            0.9340079449406524E+00,
            0.9522851197626488E+00,
            0.9661051464753107E+00,
            0.9763483833446440E+00,
            0.9837904585907746E+00,
            0.9890905016357307E+00,
            0.9927904292352575E+00,
            0.9953222650189527E+00
        };

        double[] x_vec =
        {
            0.0E+00,
            0.1E+00,
            0.2E+00,
            0.3E+00,
            0.4E+00,
            0.5E+00,
            0.6E+00,
            0.7E+00,
            0.8E+00,
            0.9E+00,
            1.0E+00,
            1.1E+00,
            1.2E+00,
            1.3E+00,
            1.4E+00,
            1.5E+00,
            1.6E+00,
            1.7E+00,
            1.8E+00,
            1.9E+00,
            2.0E+00
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

    public static void erfc_values(ref int n_data, ref double x, ref double fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ERFC_VALUES returns some values of the ERFC or "complementary error" function.
        //
        //  Discussion:
        //
        //    The complementary error function is defined by:
        //
        //      ERFC(X) = 1 - ( 2 / sqrt ( PI ) * integral ( 0 <= T <= X ) exp ( - T^2 ) dT
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Erfc[x]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 May 2007
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
        const int N_MAX = 21;

        double[] fx_vec =
        {
            1.000000000000000E+00,
            0.7772974107895215E+00,
            0.5716076449533315E+00,
            0.3961439091520741E+00,
            0.2578990352923395E+00,
            0.1572992070502851E+00,
            0.08968602177036462E+00,
            0.04771488023735119E+00,
            0.02365161665535599E+00,
            0.01090949836426929E+00,
            0.004677734981047266E+00,
            0.001862846297981891E+00,
            0.0006885138966450786E+00,
            0.0002360344165293492E+00,
            0.00007501319466545902E+00,
            0.00002209049699858544E+00,
            6.025761151762095E-06,
            1.521993362862285E-06,
            3.558629930076853E-07,
            7.700392745696413E-08,
            1.541725790028002E-08
        };

        double[] x_vec =
        {
            0.0E+00,
            0.2E+00,
            0.4E+00,
            0.6E+00,
            0.8E+00,
            1.0E+00,
            1.2E+00,
            1.4E+00,
            1.6E+00,
            1.8E+00,
            2.0E+00,
            2.2E+00,
            2.4E+00,
            2.6E+00,
            2.8E+00,
            3.0E+00,
            3.2E+00,
            3.4E+00,
            3.6E+00,
            3.8E+00,
            4.0E+00
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