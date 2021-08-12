namespace Burkardt.TestValues
{
    public static class Tsat
    {
        public static void tsat_values(ref int n_data, ref double p, ref double tc)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TSAT_VALUES returns some values of the saturation temperature.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 February 2002
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
            //    TJ270.H3, pages 16-22.
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref double P, the pressure, in bar.
            //
            //    Output, ref double TC, the saturation temperature, in
            //    degrees Celsius.
            //
        {
            int N_MAX = 20;

            double[] p_vec =
            {
                0.0061173E+00,
                0.012E+00,
                0.025E+00,
                0.055E+00,
                0.080E+00,
                0.110E+00,
                0.160E+00,
                0.250E+00,
                0.500E+00,
                0.750E+00,
                1.000E+00,
                1.500E+00,
                2.000E+00,
                5.000E+00,
                10.000E+00,
                20.000E+00,
                50.000E+00,
                100.000E+00,
                200.000E+00,
                220.550E+00
            };

            double[] tc_vec =
            {
                0.010E+00,
                9.655E+00,
                21.080E+00,
                34.589E+00,
                41.518E+00,
                47.695E+00,
                55.327E+00,
                64.980E+00,
                81.339E+00,
                91.783E+00,
                99.632E+00,
                111.378E+00,
                120.443E+00,
                151.866E+00,
                179.916E+00,
                212.417E+00,
                263.977E+00,
                311.031E+00,
                365.800E+00,
                373.976E+00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                p = 0.0;
                tc = 0.0;
            }
            else
            {
                p = p_vec[n_data - 1];
                tc = tc_vec[n_data - 1];
            }
        }

    }
}