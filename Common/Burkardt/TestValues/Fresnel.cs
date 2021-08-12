namespace Burkardt.TestValues
{
    public static class Fresnel
    {

        public static void fresnel_cos_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FRESNEL_COS_VALUES returns values of the Fresnel cosine integral function.
            //
            //  Discussion:
            //
            //    The Fresnel cosine integral is defined by:
            //
            //      C(X) = integral ( 0 <= T <= X ) cos ( PI * T^2 / 2 ) dT
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      FresnelC[x]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 November 2015
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
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 16;

            double[] fx_vec =
            {
                0.0000000000000000E+00,
                0.1999210575944531E+00,
                0.3974807591723594E+00,
                0.5810954469916523E+00,
                0.7228441718963561E+00,
                0.7798934003768228E+00,
                0.7154377229230734E+00,
                0.5430957835462564E+00,
                0.3654616834404877E+00,
                0.3336329272215571E+00,
                0.4882534060753408E+00,
                0.6362860449033195E+00,
                0.5549614058564281E+00,
                0.3889374961919690E+00,
                0.4674916516989059E+00,
                0.6057207892976856E+00
            };

            double[] x_vec =
            {
                0.0E+00,
                0.2E+00,
                0.4E+00,
                0.6E+00,
                0.8E+00,
                1.0E+00,
                1.2E+00,
                1.4E+00,
                1.6E+00,
                1.8E+00,
                2.0E+00,
                2.2E+00,
                2.4E+00,
                2.6E+00,
                2.8E+00,
                3.0E+00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void fresnel_sin_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FRESNEL_SIN_VALUES returns some values of the Fresnel sine integral function.
            //
            //  Discussion:
            //
            //    The Fresnel sine integral is defined by
            //
            //      S(X) = integral ( 0 <= T <= X ) sin ( pi * T^2 / 2 ) dT
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      FresnelS[x]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 November 2015
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
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 16;

            double[] fx_vec =
            {
                0.0000000000000000E+00,
                0.4187609161656762E-02,
                0.3335943266061318E-01,
                0.1105402073593870E+00,
                0.2493413930539178E+00,
                0.4382591473903548E+00,
                0.6234009185462497E+00,
                0.7135250773634121E+00,
                0.6388876835093809E+00,
                0.4509387692675831E+00,
                0.3434156783636982E+00,
                0.4557046121246569E+00,
                0.6196899649456836E+00,
                0.5499893231527195E+00,
                0.3915284435431718E+00,
                0.4963129989673750E+00
            };

            double[] x_vec =
            {
                0.0E+00,
                0.2E+00,
                0.4E+00,
                0.6E+00,
                0.8E+00,
                1.0E+00,
                1.2E+00,
                1.4E+00,
                1.6E+00,
                1.8E+00,
                2.0E+00,
                2.2E+00,
                2.4E+00,
                2.6E+00,
                2.8E+00,
                3.0E+00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

    }
}