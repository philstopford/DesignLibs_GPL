namespace Burkardt.Values
{
    public static class Chebyshev
    {
        public static void cheby_t_poly_values(ref int n_data, ref int n, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CHEBY_T_POLY_VALUES returns values of Chebyshev polynomials T(n,x).
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      ChebyshevT[n,x]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    10 July 2015
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
            //    Output, ref int N, the order of the function.
            //
            //    Output, ref double X, the point where the function is evaluated.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 14;

            double[] fx_vec =
            {
                0.0000000000000000E+00,
                0.1000000000000000E+01,
                0.8000000000000000E+00,
                0.2800000000000000E+00,
                -0.3520000000000000E+00,
                -0.8432000000000000E+00,
                -0.9971200000000000E+00,
                -0.7521920000000000E+00,
                -0.2063872000000000E+00,
                0.4219724800000000E+00,
                0.8815431680000000E+00,
                0.9884965888000000E+00,
                0.7000513740800000E+00,
                0.1315856097280000E+00
            };

            int[] n_vec =
            {
                -1,
                0, 1, 2,
                3, 4, 5,
                6, 7, 8,
                9, 10, 11,
                12
            };

            double[] x_vec =
            {
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                n = 0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                n = n_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void cheby_t01_poly_values(ref int n_data, ref int n, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CHEBY_T01_POLY_VALUES: values of shifted Chebyshev polynomials T01(n,x).
            //
            //  Discussion:
            //
            //    T01(n,x) = T(n,2*x-1)
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    19 July 2015
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
            //    Output, ref int N, the order of the function.
            //
            //    Output, ref double X, the point where the function is evaluated.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 25;

            double[] fx_vec =
            {
                0.0000000000000000,
                1.0000000000000000,
                0.7000000000000000,
                -0.0200000000000000,
                -0.7280000000000000,
                -0.9992000000000000,
                -0.6708800000000000,
                0.0599680000000000,
                0.7548352000000000,
                0.9968012800000000,
                0.6406865920000000,
                -0.0998400512000000,
                -0.7804626636800000,
                -0.9928076779520000,
                -1.0000000000000000,
                0.2063872000000000,
                -0.9784704000000000,
                0.2580224000000000,
                0.9870208000000000,
                0.0000000000000000,
                -0.9870208000000000,
                -0.2580224000000000,
                0.9784704000000000,
                -0.2063872000000000,
                1.0000000000000000
            };

            int[] n_vec =
            {
                -1,
                0, 1, 2,
                3, 4, 5,
                6, 7, 8,
                9, 10, 11,
                12, 7, 7,
                7, 7, 7,
                7, 7, 7,
                7, 7, 7
            };

            double[] x_vec =
            {
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.00,
                0.10,
                0.20,
                0.30,
                0.40,
                0.50,
                0.60,
                0.70,
                0.80,
                0.90,
                1.00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                n = 0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                n = n_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void cheby_u_poly_values(ref int n_data, ref int n, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CHEBY_U_POLY_VALUES returns values of Chebyshev polynomials U(n,x).
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      ChebyshevU[n,x]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    10 July 2015
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
            //    Output, ref int N, the order of the function.
            //
            //    Output, ref double X, the point where the function is evaluated.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 14;

            double[] fx_vec =
            {
                0.0000000000000000E+00,
                0.1000000000000000E+01,
                0.1600000000000000E+01,
                0.1560000000000000E+01,
                0.8960000000000000E+00,
                -0.1264000000000000E+00,
                -0.1098240000000000E+01,
                -0.1630784000000000E+01,
                -0.1511014400000000E+01,
                -0.7868390400000000E+00,
                0.2520719360000000E+00,
                0.1190154137600000E+01,
                0.1652174684160000E+01,
                0.1453325357056000E+01
            };

            int[] n_vec =
            {
                -1,
                0, 1, 2,
                3, 4, 5,
                6, 7, 8,
                9, 10, 11,
                12
            };

            double[] x_vec =
            {
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                n = 0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                n = n_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void cheby_u01_poly_values(ref int n_data, ref int n, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CHEBY_U01_POLY_VALUES: values of shifted Chebyshev polynomials U01(n,x).
            //
            //  Discussion:
            //
            //    U01(n,x) = U(n,2*x-1)
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    19 July 2015
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
            //    Output, ref int N, the order of the function.
            //
            //    Output, ref double X, the point where the function is evaluated.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 25;

            double[] fx_vec =
            {
                0.000000000000000,
                1.000000000000000,
                1.400000000000000,
                0.9600000000000000,
                -0.05600000000000000,
                -1.038400000000000,
                -1.397760000000000,
                -0.9184640000000000,
                0.1119104000000000,
                1.075138560000000,
                1.393283584000000,
                0.8754584576000000,
                -0.1676417433600000,
                -1.110156898304000,
                -8.000000000000000,
                1.511014400000000,
                -1.133260800000000,
                -0.1636352000000000,
                1.019801600000000,
                0.000000000000000,
                -1.019801600000000,
                0.1636352000000000,
                1.133260800000000,
                -1.511014400000000,
                8.000000000000000
            };

            int[] n_vec =
            {
                -1,
                0, 1, 2,
                3, 4, 5,
                6, 7, 8,
                9, 10, 11,
                12, 7, 7,
                7, 7, 7,
                7, 7, 7,
                7, 7, 7
            };

            double[] x_vec =
            {
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.00,
                0.10,
                0.20,
                0.30,
                0.40,
                0.50,
                0.60,
                0.70,
                0.80,
                0.90,
                1.00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                n = 0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                n = n_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void cheby_v_poly_values(ref int n_data, ref int n, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CHEBY_V_POLY_VALUES returns values of Chebyshev polynomials V(n,x).
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      u = Sqrt[(x+1)/2],
            //      ChebyshevT[2*n+1,u] / u
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    10 July 2015
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
            //    Output, ref int N, the order of the function.
            //
            //    Output, ref double X, the point where the function is evaluated.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 14;

            double[] fx_vec =
            {
                0.0000000000000000E+00,
                1.0000000000000000E+00,
                0.6000000000000000E+00,
                -0.0400000000000000E+00,
                -0.6640000000000000E+00,
                -1.0224000000000000E+00,
                -0.9718400000000000E+00,
                -0.5325440000000000E+00,
                0.1197696000000000E+00,
                0.7241753600000000E+00,
                1.0389109760000000E+00,
                0.9380822016000000E+00,
                0.4620205465600000E+00,
                -0.1988493271040000E+00
            };

            int[] n_vec =
            {
                -1,
                0, 1, 2,
                3, 4, 5,
                6, 7, 8,
                9, 10, 11,
                12
            };

            double[] x_vec =
            {
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                n = 0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                n = n_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void cheby_v01_poly_values(ref int n_data, ref int n, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CHEBY_V01_POLY_VALUES: values of shifted Chebyshev polynomials V01(n,x).
            //
            //  Discussion:
            //
            //    V01(n,x) = V(n,2*x-1)
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    19 July 2015
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
            //    Output, ref int N, the order of the function.
            //
            //    Output, ref double X, the point where the function is evaluated.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 25;

            double[] fx_vec =
            {
                0.0000000000000000,
                1.0000000000000000,
                0.4000000000000000,
                -0.4400000000000000,
                -1.0160000000000000,
                -0.9824000000000000,
                -0.3593600000000000,
                0.4792960000000000,
                1.0303744000000000,
                0.9632281600000000,
                0.3181450240000000,
                -0.5178251264000000,
                -1.0431002009600000,
                -0.9425151549440000,
                -15.000000000000000,
                3.1417984000000000,
                -1.3912448000000000,
                -1.2177792000000000,
                1.1837056000000000,
                1.0000000000000000,
                -0.8558976000000000,
                -0.8905088000000000,
                0.8752768000000000,
                0.1197696000000000,
                1.0000000000000000
            };

            int[] n_vec =
            {
                -1,
                0, 1, 2,
                3, 4, 5,
                6, 7, 8,
                9, 10, 11,
                12, 7, 7,
                7, 7, 7,
                7, 7, 7,
                7, 7, 7
            };

            double[] x_vec =
            {
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.00,
                0.10,
                0.20,
                0.30,
                0.40,
                0.50,
                0.60,
                0.70,
                0.80,
                0.90,
                1.00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                n = 0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                n = n_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void cheby_w_poly_values(ref int n_data, ref int n, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CHEBY_W_POLY_VALUES returns values of Chebyshev polynomials W(n,x).
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      u = Sqrt[(x+1)/2],
            //      ChebyshevU[2*n,u]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    10 July 2015
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
            //    Output, ref int N, the order of the function.
            //
            //    Output, ref double X, the point where the function is evaluated.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 14;

            double[] fx_vec =
            {
                0.000000000000000E+00,
                1.000000000000000E+00,
                2.600000000000000E+00,
                3.160000000000000E+00,
                2.456000000000000E+00,
                0.769600000000000E+00,
                -1.224640000000000E+00,
                -2.729024000000000E+00,
                -3.141798400000000E+00,
                -2.297853440000000E+00,
                -0.534767104000000E+00,
                1.442226073600000E+00,
                2.842328821760000E+00,
                3.105500041216000E+00
            };

            int[] n_vec =
            {
                -1,
                0, 1, 2,
                3, 4, 5,
                6, 7, 8,
                9, 10, 11,
                12
            };

            double[] x_vec =
            {
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                n = 0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                n = n_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void cheby_w01_poly_values(ref int n_data, ref int n, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CHEBY_W01_POLY_VALUES: values of shifted Chebyshev polynomials W01(n,x).
            //
            //  Discussion:
            //
            //    W01(n,x) = W(n,2*x-1)
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    19 July 2015
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
            //    Output, ref int N, the order of the function.
            //
            //    Output, ref double X, the point where the function is evaluated.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 25;

            double[] fx_vec =
            {
                0.000000000000000,
                1.000000000000000,
                2.400000000000000,
                2.360000000000000,
                0.904000000000000,
                -1.094400000000000,
                -2.436160000000000,
                -2.316224000000000,
                -0.806553600000000,
                1.187048960000000,
                2.468422144000000,
                2.268742041600000,
                0.707816714240000,
                -1.277798641664000,
                -1.000000000000000,
                -0.119769600000000,
                -0.875276800000000,
                0.890508800000000,
                0.855897600000000,
                -1.000000000000000,
                -1.183705600000000,
                1.217779200000000,
                1.391244800000000,
                -3.141798400000000,
                15.00000000000000
            };

            int[] n_vec =
            {
                -1,
                0, 1, 2,
                3, 4, 5,
                6, 7, 8,
                9, 10, 11,
                12, 7, 7,
                7, 7, 7,
                7, 7, 7,
                7, 7, 7
            };

            double[] x_vec =
            {
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.00,
                0.10,
                0.20,
                0.30,
                0.40,
                0.50,
                0.60,
                0.70,
                0.80,
                0.90,
                1.00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                n = 0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                n = n_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

    }
}