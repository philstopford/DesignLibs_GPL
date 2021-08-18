namespace Burkardt.Values
{
    public static class Kelvin
    {

        public static void kei0_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KEI0_VALUES returns some values of the Kelvin KEI function of order NU = 0.
            //
            //  Discussion:
            //
            //    The function is defined by:
            //
            //      KER(NU,X) + i * KEI(NU,X) = exp(-nu*Pi*I/2) * K(NU,X*exp(PI*I/4))
            //
            //    where K(NU,X) is the K Bessel function.
            //
            //    In Mathematica, KEI(NU,X) can be defined by:
            //
            //      Im [ Exp [ -NU * Pi * I / 2 ] * BesselK [ NU, X * Exp[ Pi * I / 4 ] ] ]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 June 2006
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
            int N_MAX = 11;

            double[] fx_vec =
            {
                0.0000000000000000,
                -0.6715816950943676,
                -0.4949946365187199,
                -0.3313955623385585,
                -0.2024000677647043,
                -0.1106960991556749,
                -0.05112188404598678,
                -0.01600256851827124,
                0.002198399294972520,
                0.009720918540151990,
                0.01118758650986964
            };
            double[] x_vec =
            {
                0.0,
                0.5,
                1.0,
                1.5,
                2.0,
                2.5,
                3.0,
                3.5,
                4.0,
                4.5,
                5.0
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

        public static void kei1_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KEI1_VALUES returns some values of the Kelvin KEI function of order NU = 1.
            //
            //  Discussion:
            //
            //    The function is defined by:
            //
            //      KER(NU,X) + i * KEI(NU,X) = exp(-nu*Pi*I/2) * K(NU,X*exp(PI*I/4))
            //
            //    where K(NU,X) is the K Bessel function.
            //
            //    In Mathematica, KEI(NU,X) can be defined by:
            //
            //      Im [ Exp [ -NU * Pi * I / 2 ] * BesselK [ NU, X * Exp[ Pi * I / 4 ] ] ]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 June 2006
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
            int N_MAX = 10;

            double[] fx_vec =
            {
                -1.051182085412523,
                -0.2419959664297382,
                0.001008680985009855,
                0.08004939780706674,
                0.09331378813535750,
                0.08027022252392219,
                0.05937625647622691,
                0.03916601076917133,
                0.02300216024690250,
                0.01157775439325247
            };
            double[] x_vec =
            {
                0.5,
                1.0,
                1.5,
                2.0,
                2.5,
                3.0,
                3.5,
                4.0,
                4.5,
                5.0
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

        public static void ker0_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KER0_VALUES returns some values of the Kelvin KER function of order NU = 0.
            //
            //  Discussion:
            //
            //    The function is defined by:
            //
            //      KER(NU,X) + i * KEI(NU,X) = exp(-nu*Pi*I/2) * K(NU,X*exp(PI*I/4))
            //
            //    where K(NU,X) is the K Bessel function.
            //
            //    In Mathematica, KER(NU,X) can be defined by:
            //
            //      Re [ Exp [ -NU * Pi * I / 2 ] * BesselK [ NU, X * Exp[ Pi * I / 4 ] ] ]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 June 2006
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
            int N_MAX = 10;

            double[] fx_vec =
            {
                0.8559058721186342,
                0.2867062087283160,
                0.05293491548771044,
                -0.04166451399150953,
                -0.06968797258904534,
                -0.06702923330379870,
                -0.05263927724224119,
                -0.03617884789954761,
                -0.02199987504667382,
                -0.01151172719949066
            };
            double[] x_vec =
            {
                0.5,
                1.0,
                1.5,
                2.0,
                2.5,
                3.0,
                3.5,
                4.0,
                4.5,
                5.0
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

        public static void ker1_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KER1_VALUES returns some values of the Kelvin KER function of order NU = 1.
            //
            //  Discussion:
            //
            //    The function is defined by:
            //
            //      KER(NU,X) + i * KEI(NU,X) = exp(-nu*Pi*I/2) * K(NU,X*exp(PI*I/4))
            //
            //    where K(NU,X) is the K Bessel function.
            //
            //    In Mathematica, KER(NU,X) can be defined by:
            //
            //      Re [ Exp [ -NU * Pi * I / 2 ] * BesselK [ NU, X * Exp[ Pi * I / 4 ] ] ]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 June 2006
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
            int N_MAX = 10;

            double[] fx_vec =
            {
                -1.522403406532090,
                -0.7403222768419827,
                -0.4170442851662574,
                -0.2308059295181230,
                -0.1172561358598705,
                -0.04989830778751491,
                -0.01272324936181659,
                0.005351296460277448,
                0.01209090413515866,
                0.01273739048421857
            };
            double[] x_vec =
            {
                0.5,
                1.0,
                1.5,
                2.0,
                2.5,
                3.0,
                3.5,
                4.0,
                4.5,
                5.0
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