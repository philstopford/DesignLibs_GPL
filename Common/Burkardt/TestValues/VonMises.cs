namespace Burkardt.TestValues
{
    public static class VonMises
    {
        public static void von_mises_cdf_values(ref int n_data, ref double a, ref double b, ref double x,
                ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    VON_MISES_CDF_VALUES returns some values of the von Mises CDF.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    08 December 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Kanti Mardia, Peter Jupp,
            //    Directional Statistics,
            //    Wiley, 2000, QA276.M335
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref double A, &B, the parameters of the function.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 23;

            double[] a_vec =
            {
                0.0E+00,
                0.0E+00,
                0.0E+00,
                0.0E+00,
                0.0E+00,
                0.1E+01,
                0.1E+01,
                0.1E+01,
                0.1E+01,
                0.1E+01,
                0.1E+01,
                -0.2E+01,
                -0.1E+01,
                0.0E+01,
                0.1E+01,
                0.2E+01,
                0.3E+01,
                0.0E+00,
                0.0E+00,
                0.0E+00,
                0.0E+00,
                0.0E+00,
                0.0E+00
            };

            double[] b_vec =
            {
                0.1E+01,
                0.1E+01,
                0.1E+01,
                0.1E+01,
                0.1E+01,
                0.2E+01,
                0.2E+01,
                0.2E+01,
                0.2E+01,
                0.2E+01,
                0.2E+01,
                0.3E+01,
                0.3E+01,
                0.3E+01,
                0.3E+01,
                0.3E+01,
                0.3E+01,
                0.0E+00,
                0.1E+01,
                0.2E+01,
                0.3E+01,
                0.4E+01,
                0.5E+01
            };

            double[] fx_vec =
            {
                0.2535089956281180E-01,
                0.1097539041177346E+00,
                0.5000000000000000E+00,
                0.8043381312498558E+00,
                0.9417460124555197E+00,
                0.5000000000000000E+00,
                0.6018204118446155E+00,
                0.6959356933122230E+00,
                0.7765935901304593E+00,
                0.8410725934916615E+00,
                0.8895777369550366E+00,
                0.9960322705517925E+00,
                0.9404336090170247E+00,
                0.5000000000000000E+00,
                0.5956639098297530E-01,
                0.3967729448207649E-02,
                0.2321953958111930E-03,
                0.6250000000000000E+00,
                0.7438406999109122E+00,
                0.8369224904294019E+00,
                0.8941711407897124E+00,
                0.9291058600568743E+00,
                0.9514289900655436E+00
            };

            double[] x_vec =
            {
                -0.2617993977991494E+01,
                -0.1570796326794897E+01,
                0.0000000000000000E+00,
                0.1047197551196598E+01,
                0.2094395102393195E+01,
                0.1000000000000000E+01,
                0.1200000000000000E+01,
                0.1400000000000000E+01,
                0.1600000000000000E+01,
                0.1800000000000000E+01,
                0.2000000000000000E+01,
                0.0000000000000000E+00,
                0.0000000000000000E+00,
                0.0000000000000000E+00,
                0.0000000000000000E+00,
                0.0000000000000000E+00,
                0.0000000000000000E+00,
                0.7853981633974483E+00,
                0.7853981633974483E+00,
                0.7853981633974483E+00,
                0.7853981633974483E+00,
                0.7853981633974483E+00,
                0.7853981633974483E+00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                a = 0.0;
                b = 0.0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                a = a_vec[n_data - 1];
                b = b_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

    }
}