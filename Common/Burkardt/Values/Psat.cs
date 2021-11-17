namespace Burkardt.Values;

public static class Psat
{
    public static void psat_values(ref int n_data, ref double tc, ref double p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PSAT_VALUES returns some values of the saturation pressure.
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
        //    TJ270.H3, pages 9-15.
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
        //    Output, ref double P, the saturation pressure, in bar.
        //
    {
        const int N_MAX = 12;

        double[] p_vec =
        {
            0.0061173E+00,
            0.0065716E+00,
            0.0087260E+00,
            0.12344E+00,
            1.0132E+00,
            2.3201E+00,
            4.7572E+00,
            15.537E+00,
            39.737E+00,
            85.838E+00,
            165.21E+00,
            220.55E+00
        };

        double[] tc_vec =
        {
            0.100000E-01,
            0.100000E+01,
            0.500000E+01,
            0.500000E+02,
            0.100000E+03,
            0.125000E+03,
            0.150000E+03,
            0.200000E+03,
            0.250000E+03,
            0.300000E+03,
            0.350000E+03,
            0.373976E+03
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
        }
        else
        {
            tc = tc_vec[n_data - 1];
            p = p_vec[n_data - 1];
        }
    }
}