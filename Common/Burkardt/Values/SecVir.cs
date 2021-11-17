﻿namespace Burkardt.Values;

public static class SecVir
{
    public static void secvir_values(ref int n_data, ref double tc, ref double vir)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SECVIR_VALUES returns some values of the second virial coefficient.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 February 2002
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
        //    TJ270.H3, pages 24-25.
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
        //    Output, ref double VIR, the second virial coefficient, in
        //    m^3/kg.
        //
    {
        const int N_MAX = 19;

        double[] tc_vec =
        {
            0.0E+00,
            5.0E+00,
            10.0E+00,
            20.0E+00,
            30.0E+00,
            40.0E+00,
            60.0E+00,
            90.0E+00,
            120.0E+00,
            150.0E+00,
            180.0E+00,
            210.0E+00,
            240.0E+00,
            300.0E+00,
            400.0E+00,
            500.0E+00,
            700.0E+00,
            1000.0E+00,
            2000.0E+00
        };

        double[] vir_vec =
        {
            -98.96E+00,
            -90.08E+00,
            -82.29E+00,
            -69.36E+00,
            -59.19E+00,
            -51.07E+00,
            -39.13E+00,
            -27.81E+00,
            -20.83E+00,
            -16.21E+00,
            -12.98E+00,
            -10.63E+00,
            -8.85E+00,
            -6.39E+00,
            -4.03E+00,
            -2.71E+00,
            -1.32E+00,
            -0.39E+00,
            0.53E+00
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
            vir = 0.0;
        }
        else
        {
            tc = tc_vec[n_data - 1];
            vir = vir_vec[n_data - 1];
        }
    }

}