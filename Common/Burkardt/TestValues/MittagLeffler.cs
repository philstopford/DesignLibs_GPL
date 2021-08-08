namespace TestValues
{
    public static class MittagLeffler
    {

        public static void mittag_leffler_ea_values(ref int n_data, ref double a, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MITTAG_LEFFLER_EA_VALUES: values of one-parameter Mittag Leffler function.
            //
            //  Discussion:
            //
            //    E(alpha;z) = sum ( 0 <= k < oo ) z^k / Gamma ( alpha * k + 1 )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    02 February 2017
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref double A, the parameter of the function.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 25;

            double[] a_vec =
            {
                0.25,
                0.25,
                0.25,
                0.25,
                0.25,
                1.75,
                1.75,
                1.75,
                1.75,
                1.75,
                2.25,
                2.25,
                2.25,
                2.25,
                2.25,
                1.0,
                2.0,
                3.0,
                4.0,
                5.0,
                6.0,
                7.0,
                8.0,
                9.0,
                10.0
            };

            double[] fx_vec =
            {
                1.0,
                0.2525646348870172,
                0.1427989464258737,
                0.994231091392936E-01,
                0.762370352397216E-01,
                1.0,
                -0.2579027070618285,
                0.2716853059670252,
                0.1846579916604009E-01,
                -0.139707389642194,
                1.0,
                -0.1022121843823497E+01,
                0.1860844522589611E+01,
                0.2644615445996891E+01,
                0.7762512036620307,
                0.6737946999045729E-02,
                -0.6172728764571668,
                0.2010457248089053,
                0.7922864454196143,
                0.958340222567225,
                0.993055607747429,
                0.999007936794713,
                0.999875992064687,
                0.999986221340384,
                0.999998622134038
            };

            double[] x_vec =
            {
                1.0,
                0.2525646348870172,
                0.1427989464258737,
                0.994231091392936E-01,
                0.762370352397216E-01,
                1.0,
                -0.2579027070618285,
                0.2716853059670252,
                0.1846579916604009E-01,
                -0.139707389642194,
                1.0,
                -0.1022121843823497E+01,
                0.1860844522589611E+01,
                0.2644615445996891E+01,
                0.7762512036620307,
                0.6737946999045729E-02,
                -0.6172728764571668,
                0.2010457248089053,
                0.7922864454196143,
                0.958340222567225,
                0.993055607747429,
                0.999007936794713,
                0.999875992064687,
                0.999986221340384,
                0.999998622134038
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
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                a = a_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void mittag_leffler_eab_values(ref int n_data, ref double a, ref double b, ref double x,
                ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MITTAG_LEFFLER_EAB_VALUES: values of two-parameter Mittag Leffler function.
            //
            //  Discussion:
            //
            //    E(alpha,beta;z) = sum ( 0 <= k < oo ) z^k / Gamma ( alpha * k + beta )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 February 2017
            //
            //  Author:
            //
            //    John Burkardt
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
            int N_MAX = 18;

            double[] a_vec =
            {
                1.0E+00,
                1.0E+00,
                1.0E+00,
                1.5E+00,
                1.5E+00,
                1.5E+00,
                3.0E+00,
                3.0E+00,
                3.0E+00,
                1.5E+00,
                1.5E+00,
                1.5E+00,
                1.5E+00,
                2.0E+00,
                2.5E+00,
                2.0E+00,
                2.0E+00,
                2.0E+00
            };

            double[] b_vec =
            {
                0.0E+00,
                2.5E+00,
                5.0E+00,
                0.0E+00,
                5.0E+00,
                10.0E+00,
                1.0E+00,
                1.1E+00,
                1.2E+00,
                5.0E+00,
                10.0E+00,
                15.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                1.0E+00,
                2.0E+00,
                3.0E+00
            };

            double[] fx_vec =
            {
                -0.1947001957678512E+00,
                0.6821152851737027E+00,
                0.3966713294631294E-01,
                -0.7783327874361643E+00,
                0.376189576668972E-01,
                0.2653864659692168E-05,
                0.5521129475536116E+00,
                0.6561318862220054E+00,
                0.7417555514147703E+00,
                0.283627999653323E-01,
                0.2382766080122566E-05,
                0.1057337628882522E-10,
                0.6524069073077504E-01,
                0.4926693884523065E-01,
                0.4440848516337653E-01,
                0.217818355660857E+01,
                0.1368298872008591E+01,
                0.5890917783042855E+00
            };

            double[] x_vec =
            {
                -0.25E+00,
                -0.25E+00,
                -0.25E+00,
                -1.25E+00,
                -1.25E+00,
                -1.25E+00,
                -2.75E+00,
                -2.75E+00,
                -2.75E+00,
                -5.0E+00,
                -5.0E+00,
                -5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                2.0E+00,
                2.0E+00,
                2.0E+00
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