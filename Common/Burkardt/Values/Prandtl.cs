﻿namespace Burkardt.Values;

public static class Prandtl
{

    public static void prandtl_values(ref int n_data, ref double tc, ref double p, ref double pr)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PRANDTL_VALUES returns some values of the Prandtl number.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 February 2002
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
        //    TJ270.H3, page 265.
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
        //    Output, ref double PR, the Prandtl number, dimensionless.
        //
    {
        const int N_MAX = 35;

        double[] pr_vec =
        {
            13.50E+00,
            13.48E+00,
            13.46E+00,
            13.39E+00,
            13.27E+00,
            13.15E+00,
            13.04E+00,
            12.93E+00,
            12.83E+00,
            12.73E+00,
            12.63E+00,
            12.53E+00,
            12.43E+00,
            12.34E+00,
            12.25E+00,
            12.08E+00,
            11.92E+00,
            11.77E+00,
            11.62E+00,
            11.48E+00,
            11.36E+00,
            11.23E+00,
            11.12E+00,
            10.91E+00,
            10.72E+00,
            10.55E+00,
            6.137E+00,
            3.555E+00,
            2.378E+00,
            1.000E+00,
            0.974E+00,
            0.960E+00,
            0.924E+00,
            0.899E+00,
            0.882E+00
        };

        double[] p_vec =
        {
            1.0E+00,
            5.0E+00,
            10.0E+00,
            25.0E+00,
            50.0E+00,
            75.0E+00,
            100.0E+00,
            125.0E+00,
            150.0E+00,
            175.0E+00,
            200.0E+00,
            225.0E+00,
            250.0E+00,
            275.0E+00,
            300.0E+00,
            350.0E+00,
            400.0E+00,
            450.0E+00,
            500.0E+00,
            550.0E+00,
            600.0E+00,
            650.0E+00,
            700.0E+00,
            800.0E+00,
            900.0E+00,
            1000.0E+00,
            1.0E+00,
            1.0E+00,
            1.0E+00,
            1.0E+00,
            1.0E+00,
            1.0E+00,
            1.0E+00,
            1.0E+00,
            1.0E+00
        };

        double[] tc_vec =
        {
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
            0.0E+00,
            0.0E+00,
            0.0E+00,
            0.0E+00,
            25.0E+00,
            50.0E+00,
            75.0E+00,
            100.0E+00,
            150.0E+00,
            200.0E+00,
            400.0E+00,
            600.0E+00,
            800.0E+00
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
            pr = 0.0;
        }
        else
        {
            tc = tc_vec[n_data - 1];
            p = p_vec[n_data - 1];
            pr = pr_vec[n_data - 1];
        }
    }

}