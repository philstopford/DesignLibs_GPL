namespace Burkardt.TestValues
{
    public static class SphericalHarmonic
    {
        public static void spherical_harmonic_values(ref int n_data, ref int l, ref int m, ref double theta,
                ref double phi, ref double yr, ref double yi)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERICAL_HARMONIC_VALUES returns values of spherical harmonic functions.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by
            //
            //      SphericalHarmonicY [ l, m, theta, phi ]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 March 2005
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
            //    Eric Weisstein,
            //    CRC Concise Encyclopedia of Mathematics,
            //    CRC Press, 1998.
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
            //    Output, ref int L, ref int M, ref double THETA, &PHI, the arguments
            //    of the function.
            //
            //    Output, ref double YR, &YI, the real and imaginary parts of
            //    the function.
            //
        {
            int N_MAX = 20;

            int[] l_vec =
            {
                0, 1, 2,
                3, 4, 5,
                5, 5, 5,
                5, 4, 4,
                4, 4, 4,
                3, 3, 3,
                3, 3
            };
            int[] m_vec =
            {
                0, 0, 1,
                2, 3, 5,
                4, 3, 2,
                1, 2, 2,
                2, 2, 2,
                -1, -1, -1,
                -1, -1
            };
            double[] phi_vec =
            {
                0.1047197551196598E+01, 0.1047197551196598E+01, 0.1047197551196598E+01,
                0.1047197551196598E+01, 0.1047197551196598E+01, 0.6283185307179586E+00,
                0.6283185307179586E+00, 0.6283185307179586E+00, 0.6283185307179586E+00,
                0.6283185307179586E+00, 0.7853981633974483E+00, 0.7853981633974483E+00,
                0.7853981633974483E+00, 0.7853981633974483E+00, 0.7853981633974483E+00,
                0.4487989505128276E+00, 0.8975979010256552E+00, 0.1346396851538483E+01,
                0.1795195802051310E+01, 0.2243994752564138E+01
            };
            double[] theta_vec =
            {
                0.5235987755982989E+00, 0.5235987755982989E+00, 0.5235987755982989E+00,
                0.5235987755982989E+00, 0.5235987755982989E+00, 0.2617993877991494E+00,
                0.2617993877991494E+00, 0.2617993877991494E+00, 0.2617993877991494E+00,
                0.2617993877991494E+00, 0.6283185307179586E+00, 0.1884955592153876E+01,
                0.3141592653589793E+01, 0.4398229715025711E+01, 0.5654866776461628E+01,
                0.3926990816987242E+00, 0.3926990816987242E+00, 0.3926990816987242E+00,
                0.3926990816987242E+00, 0.3926990816987242E+00
            };
            double[] yi_vec =
            {
                0.0000000000000000E+00, 0.0000000000000000E+00, -0.2897056515173922E+00,
                0.1916222768312404E+00, 0.0000000000000000E+00, 0.0000000000000000E+00,
                0.3739289485283311E-02, -0.4219517552320796E-01, 0.1876264225575173E+00,
                -0.3029973424491321E+00, 0.4139385503112256E+00, -0.1003229830187463E+00,
                0.0000000000000000E+00, -0.1003229830187463E+00, 0.4139385503112256E+00,
                -0.1753512375142586E+00, -0.3159720118970196E+00, -0.3940106541811563E+00,
                -0.3940106541811563E+00, -0.3159720118970196E+00
            };
            double[] yr_vec =
            {
                0.2820947917738781E+00, 0.4231421876608172E+00, -0.1672616358893223E+00,
                -0.1106331731112457E+00, 0.1354974113737760E+00, 0.5390423109043568E-03,
                -0.5146690442951909E-02, 0.1371004361349490E-01, 0.6096352022265540E-01,
                -0.4170400640977983E+00, 0.0000000000000000E+00, 0.0000000000000000E+00,
                0.0000000000000000E+00, 0.0000000000000000E+00, 0.0000000000000000E+00,
                0.3641205966137958E+00, 0.2519792711195075E+00, 0.8993036065704300E-01,
                -0.8993036065704300E-01, -0.2519792711195075E+00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                l = 0;
                m = 0;
                theta = 0.0;
                phi = 0.0;
                yr = 0.0;
                yi = 0.0;
            }
            else
            {
                l = l_vec[n_data - 1];
                m = m_vec[n_data - 1];
                theta = theta_vec[n_data - 1];
                phi = phi_vec[n_data - 1];
                yr = yr_vec[n_data - 1];
                yi = yi_vec[n_data - 1];
            }
        }

    }
}