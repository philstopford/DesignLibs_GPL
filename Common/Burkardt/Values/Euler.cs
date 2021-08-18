namespace Burkardt.Values
{
    public static class Euler
    {
        public static void euler_number_values(ref int n_data, ref int n, ref int c)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    EULER_NUMBER_VALUES returns some values of the Euler numbers.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      EulerE[n]
            //
            //    These numbers rapidly get too big to store in an ordinary integer!
            //
            //    The terms of odd index are 0.
            //
            //    E(N) = -C(N,N-2) * E(N-2) - C(N,N-4) * E(N-4) - ... - C(N,0) * E(0).
            //
            //  First terms:
            //
            //    E0  = 1
            //    E1  = 0;
            //    E2  = -1
            //    E3  = 0;
            //    E4  = 5
            //    E5  = 0;
            //    E6  = -61
            //    E7  = 0;
            //    E8  = 1385
            //    E9  = 0;
            //    E10 = -50521
            //    E11 = 0;
            //    E12 = 2702765
            //    E13 = 0;
            //    E14 = -199360981
            //    E15 = 0;
            //    E16 = 19391512145
            //    E17 = 0;
            //    E18 = -2404879675441
            //    E19 = 0;
            //    E20 = 370371188237525
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 February 2015
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
            //    Output, ref int N, the order of the Euler number.
            //
            //    Output, ref int C, the value of the Euler number.
            //
        {
            int N_MAX = 8;

            int[] c_vec =
            {
                1, 0, -1, 5, -61, 1385, -50521, 2702765
            };

            int[] n_vec =
            {
                0, 1, 2, 4, 6, 8, 10, 12
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
                c = 0;
            }
            else
            {
                n = n_vec[n_data - 1];
                c = c_vec[n_data - 1];
            }
        }

        public static void euler_poly_values(ref int n_data, ref int n, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    EULER_POLY_VALUES returns some values of the Euler polynomials.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      EulerE[n,x]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 August 2004
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
            //    Output, ref int N, the order of the Euler polynomial.
            //
            //    Output, ref double X, the argument of the Euler polynomial.
            //
            //    Output, ref double FX, the value of the Euler polynomial.
            //
        {
            int N_MAX = 27;

            double[] fx_vec =
            {
                0.100000000000E+01,
                -0.300000000000E+00,
                -0.160000000000E+00,
                0.198000000000E+00,
                0.185600000000E+00,
                -0.403680000000E+00,
                -0.560896000000E+00,
                0.171878880000E+01,
                0.318043136000E+01,
                -0.125394670080E+02,
                -0.289999384576E+02,
                -0.625000000000E-01,
                -0.174240000000E+00,
                -0.297680000000E+00,
                -0.404320000000E+00,
                -0.475260000000E+00,
                -0.500000000000E+00,
                -0.475240000000E+00,
                -0.403680000000E+00,
                -0.292820000000E+00,
                -0.153760000000E+00,
                0.000000000000E+00,
                0.153760000000E+00,
                0.292820000000E+00,
                0.403680000000E+00,
                0.475240000000E+00,
                0.500000000000E+00
            };

            int[] n_vec =
            {
                0,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8,
                9,
                10,
                5,
                5,
                5,
                5,
                5,
                5,
                5,
                5,
                5,
                5,
                5,
                5,
                5,
                5,
                5,
                5
            };

            double[] x_vec =
            {
                0.2E+00,
                0.2E+00,
                0.2E+00,
                0.2E+00,
                0.2E+00,
                0.2E+00,
                0.2E+00,
                0.2E+00,
                0.2E+00,
                0.2E+00,
                0.2E+00,
                -0.5E+00,
                -0.4E+00,
                -0.3E+00,
                -0.2E+00,
                -0.1E+00,
                0.0E+00,
                0.1E+00,
                0.2E+00,
                0.3E+00,
                0.4E+00,
                0.5E+00,
                0.6E+00,
                0.7E+00,
                0.8E+00,
                0.9E+00,
                1.0E+00
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