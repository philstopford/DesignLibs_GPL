namespace Burkardt.Values
{
    public static class ExpIntegral
    {
        public static void e1_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    E1_VALUES returns some values of the exponential integral function E1(X).
            //
            //  Definition:
            //
            //    The exponential integral E1(X) is defined by the formula:
            //
            //      E1(X) = integral ( 1 <= T <= +oo ) exp ( -X*T ) / T dT
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      ExpIntegralE[1,x]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    20 August 2004
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
                0.5597735947761608E+00,
                0.4543795031894021E+00,
                0.3737688432335091E+00,
                0.3105965785455430E+00,
                0.2601839393259996E+00,
                0.2193839343955203E+00,
                0.1859909045360402E+00,
                0.1584084368514626E+00,
                0.1354509578491291E+00,
                0.1162193125713579E+00,
                0.1000195824066327E+00,
                0.8630833369753979E-01,
                0.7465464440125305E-01,
                0.6471312936386886E-01,
                0.5620437817453485E-01,
                0.4890051070806112E-01
            };

            double[] x_vec =
            {
                0.5E+00,
                0.6E+00,
                0.7E+00,
                0.8E+00,
                0.9E+00,
                1.0E+00,
                1.1E+00,
                1.2E+00,
                1.3E+00,
                1.4E+00,
                1.5E+00,
                1.6E+00,
                1.7E+00,
                1.8E+00,
                1.9E+00,
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
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void ei_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    EI_VALUES returns some values of the exponential integral function EI(X).
            //
            //  Definition:
            //
            //    The exponential integral EI(X) has the formula:
            //
            //      EI(X) = - integral ( -X <= T < +oo ) exp ( -T ) / T dT
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      ExpIntegralEi[x]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    19 August 2004
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
                0.4542199048631736E+00,
                0.7698812899373594E+00,
                0.1064907194624291E+01,
                0.1347396548212326E+01,
                0.1622811713696867E+01,
                0.1895117816355937E+01,
                0.2167378279563403E+01,
                0.2442092285192652E+01,
                0.2721398880232024E+01,
                0.3007207464150646E+01,
                0.3301285449129798E+01,
                0.3605319949019469E+01,
                0.3920963201354904E+01,
                0.4249867557487934E+01,
                0.4593713686953585E+01,
                0.4954234356001890E+01
            };

            double[] x_vec =
            {
                0.5E+00,
                0.6E+00,
                0.7E+00,
                0.8E+00,
                0.9E+00,
                1.0E+00,
                1.1E+00,
                1.2E+00,
                1.3E+00,
                1.4E+00,
                1.5E+00,
                1.6E+00,
                1.7E+00,
                1.8E+00,
                1.9E+00,
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