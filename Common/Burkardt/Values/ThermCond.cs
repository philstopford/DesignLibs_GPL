namespace Burkardt.Values;

public static class ThermCond
{
    public static void thercon_values(ref int n_data, ref double tc, ref double p, ref double lambda)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    THERCON_VALUES returns some values of the thermal conductivity.
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
        //    TJ270.H3, page 264.
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
        //    Output, ref double LAMBDA, the thermal conductivity, in
        //    mW/(m degrees Kelvin).
        //
    {
        const int N_MAX = 35;

        double[] lambda_vec =
        {
            561.00E+00,
            561.30E+00,
            561.50E+00,
            562.40E+00,
            563.70E+00,
            565.10E+00,
            566.50E+00,
            567.90E+00,
            569.30E+00,
            570.60E+00,
            572.00E+00,
            573.40E+00,
            574.80E+00,
            576.10E+00,
            577.50E+00,
            580.20E+00,
            582.90E+00,
            585.50E+00,
            588.10E+00,
            590.70E+00,
            593.30E+00,
            595.80E+00,
            598.30E+00,
            603.10E+00,
            607.80E+00,
            612.20E+00,
            607.20E+00,
            643.60E+00,
            666.80E+00,
            25.08E+00,
            28.85E+00,
            33.28E+00,
            54.76E+00,
            79.89E+00,
            107.30E+00
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
            lambda = 0.0;
        }
        else
        {
            tc = tc_vec[n_data - 1];
            p = p_vec[n_data - 1];
            lambda = lambda_vec[n_data - 1];
        }
    }

}