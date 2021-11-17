﻿namespace Burkardt.Values;

public static class Viscosity
{
    public static void viscosity_values(ref int n_data, ref double tc, ref double p, ref double eta)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VISCOSITY_VALUES returns some values of the viscosity function.
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
        //    TJ270.H3, page 263.
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
        //    Output, ref double ETA, the viscosity, in MegaPascal seconds.
        //
    {
        const int N_MAX = 34;

        double[] eta_vec =
        {
            1792.0E+00,
            1791.0E+00,
            1790.0E+00,
            1786.0E+00,
            1780.0E+00,
            1775.0E+00,
            1769.0E+00,
            1764.0E+00,
            1759.0E+00,
            1754.0E+00,
            1749.0E+00,
            1744.0E+00,
            1739.0E+00,
            1735.0E+00,
            1731.0E+00,
            1722.0E+00,
            1714.0E+00,
            1707.0E+00,
            1700.0E+00,
            1694.0E+00,
            1687.0E+00,
            1682.0E+00,
            1676.0E+00,
            1667.0E+00,
            1659.0E+00,
            1653.0E+00,
            890.8E+00,
            547.1E+00,
            378.4E+00,
            12.28E+00,
            16.18E+00,
            24.45E+00,
            32.61E+00,
            40.38E+00
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
            eta = 0.0;
        }
        else
        {
            tc = tc_vec[n_data - 1];
            p = p_vec[n_data - 1];
            eta = eta_vec[n_data - 1];
        }
    }

}