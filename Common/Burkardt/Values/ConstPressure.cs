﻿namespace Burkardt.Values;

public static class ConstPressure
{
    public static void cp_values(ref int n_data, ref double tc, ref double p, ref double cp)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CP_VALUES returns some values of the specific heat at constant pressure.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 August 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Lester Haar, John Gallagher and George Kell,
        //    NBS/NRC Steam Tables:
        //    Thermodynamic and Transport Properties and Computer Programs
        //    for Vapor and Liquid States of Water in SI Units,
        //    Hemisphere Publishing Corporation, Washington, 1984,
        //    TJ270.H3, pages 229-237.
        //
        //  Parameters:
        //
        //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, ref double TC, the temperature, in degrees Celsius.
        //
        //    Output, ref double P, the pressure, in bar.
        //
        //    Output, ref double CP, the specific heat at constant pressure,
        //    in KJ/(kg K).
        //
    {
        const int N_MAX = 24;

        double[] cp_vec =
        {
            4.228E+00,
            2.042E+00,
            1.975E+00,
            2.013E+00,
            2.040E+00,
            2.070E+00,
            2.135E+00,
            2.203E+00,
            2.378E+00,
            2.541E+00,
            2.792E+00,
            2.931E+00,
            4.226E+00,
            4.223E+00,
            4.202E+00,
            4.177E+00,
            4.130E+00,
            4.089E+00,
            4.053E+00,
            4.021E+00,
            3.909E+00,
            3.844E+00,
            3.786E+00,
            2.890E+00
        };

        double[] p_vec =
        {
            1.0E+00,
            1.0E+00,
            1.0E+00,
            1.0E+00,
            1.0E+00,
            1.0E+00,
            1.0E+00,
            1.0E+00,
            1.0E+00,
            1.0E+00,
            1.0E+00,
            1.0E+00,
            5.0E+00,
            10.0E+00,
            50.0E+00,
            100.0E+00,
            200.0E+00,
            300.0E+00,
            400.0E+00,
            500.0E+00,
            1000.0E+00,
            1500.0E+00,
            2000.0E+00,
            5000.0E+00
        };

        double[] tc_vec =
        {
            0.0E+00,
            100.0E+00,
            200.0E+00,
            300.0E+00,
            350.0E+00,
            400.0E+00,
            500.0E+00,
            600.0E+00,
            850.0E+00,
            1100.0E+00,
            1600.0E+00,
            2000.0E+00,
            0.0E+00,
            0.0E+00,
            0.0E+00,
            0.0E+00,
            0.0E+00,
            0.0E+00,
            0.0E+00,
            0.0E+00,
            0.0E+00,
            0.0E+00,
            0.0E+00,
            0.0E+00
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
            tc = 0.0;
            p = 0.0;
            cp = 0.0;
        }
        else
        {
            tc = tc_vec[n_data - 1];
            p = p_vec[n_data - 1];
            cp = cp_vec[n_data - 1];
        }
    }

}