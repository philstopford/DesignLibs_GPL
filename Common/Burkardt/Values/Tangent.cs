namespace Burkardt.Values;

public static class Tangent
{
    public static void tan_values(ref int n_data, ref double x, ref double fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TAN_VALUES returns some values of the tangent function.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Tan[x]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 June 2007
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
        const int N_MAX = 15;

        double[] fx_vec =
        {
            0.00000000000000000000,
            0.26794919243112270647,
            0.54630248984379051326,
            0.57735026918962576451,
            1.0000000000000000000,
            1.5574077246549022305,
            1.7320508075688772935,
            3.7320508075688772935,
            7.5957541127251504405,
            15.257051688265539110,
            -2.1850398632615189916,
            -0.14254654307427780530,
            0.0000000000000000000,
            1.1578212823495775831,
            -3.3805150062465856370
        };

        double[] x_vec =
        {
            0.00000000000000000000,
            0.26179938779914943654,
            0.50000000000000000000,
            0.52359877559829887308,
            0.78539816339744830962,
            1.0000000000000000000,
            1.0471975511965977462,
            1.3089969389957471827,
            1.4398966328953219010,
            1.5053464798451092601,
            2.0000000000000000000,
            3.0000000000000000000,
            3.1415926535897932385,
            4.0000000000000000000,
            5.0000000000000000000
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

    public static void tanh_values(ref int n_data, ref double x, ref double fx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TANH_VALUES returns some values of the hyperbolic tangent function.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Tanh[x]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 June 2007
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
        const int N_MAX = 18;

        double[] fx_vec =
        {
            -0.99990920426259513121,
            -0.76159415595576488812,
            0.00000000000000000000,
            0.099667994624955817118,
            0.19737532022490400074,
            0.29131261245159090582,
            0.37994896225522488527,
            0.46211715726000975850,
            0.53704956699803528586,
            0.60436777711716349631,
            0.66403677026784896368,
            0.71629787019902442081,
            0.76159415595576488812,
            0.96402758007581688395,
            0.99505475368673045133,
            0.99932929973906704379,
            0.99990920426259513121,
            0.99999999587769276362
        };

        double[] x_vec =
        {
            -5.0,
            -1.0,
            0.0,
            0.1,
            0.2,
            0.3,
            0.4,
            0.5,
            0.6,
            0.7,
            0.8,
            0.9,
            1.0,
            2.0,
            3.0,
            4.0,
            5.0,
            10.0
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