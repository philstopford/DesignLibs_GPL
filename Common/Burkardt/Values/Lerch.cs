namespace Burkardt.Values
{
    public static class Lerch
    {

        public static void lerch_values(ref int n_data, ref double z, ref int s, ref double a, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LERCH_VALUES returns some values of the Lerch transcendent function.
            //
            //  Discussion:
            //
            //    The Lerch function is defined as
            //
            //      Phi(z,s,a) = Sum ( 0 <= k < +oo ) z^k / ( a + k )^s
            //
            //    omitting any terms with ( a + k ) = 0.
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      LerchPhi[z,s,a]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 September 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
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
            //    Output, ref double Z, the parameters of the function.
            //
            //    Output, ref int S, the parameters of the function.
            //
            //    Output, ref double A, the parameters of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 12;

            double[] a_vec =
            {
                0.0E+00,
                0.0E+00,
                0.0E+00,
                1.0E+00,
                1.0E+00,
                1.0E+00,
                2.0E+00,
                2.0E+00,
                2.0E+00,
                3.0E+00,
                3.0E+00,
                3.0E+00
            };

            double[] fx_vec =
            {
                0.1644934066848226E+01,
                0.1202056903159594E+01,
                0.1000994575127818E+01,
                0.1164481052930025E+01,
                0.1074426387216080E+01,
                0.1000492641212014E+01,
                0.2959190697935714E+00,
                0.1394507503935608E+00,
                0.9823175058446061E-03,
                0.1177910993911311E+00,
                0.3868447922298962E-01,
                0.1703149614186634E-04
            };

            int[] s_vec =
            {
                2, 3, 10,
                2, 3, 10,
                2, 3, 10,
                2, 3, 10
            };

            double[] z_vec =
            {
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.5000000000000000E+00,
                0.5000000000000000E+00,
                0.5000000000000000E+00,
                0.3333333333333333E+00,
                0.3333333333333333E+00,
                0.3333333333333333E+00,
                0.1000000000000000E+00,
                0.1000000000000000E+00,
                0.1000000000000000E+00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                z = 0.0;
                s = 0;
                a = 0.0;
                fx = 0.0;
            }
            else
            {
                z = z_vec[n_data - 1];
                s = s_vec[n_data - 1];
                a = a_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

    }
}