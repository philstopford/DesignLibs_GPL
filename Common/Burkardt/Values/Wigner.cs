namespace Burkardt.Values
{
    public static class Wigner
    {
        public static void nine_j_values(ref int n_data, ref double j1, ref double j2, ref double j3,
                ref double j4, ref double j5, ref double j6, ref double j7, ref double j8, ref double j9,
                ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NINE_J_VALUES returns some values of the Wigner 9J function.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    09 February 2007
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
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref double J1, &J2, &J3, &J4, &J5, &J6, &J7, &J8, &J9,
            //    the arguments of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 9;

            double[] fx_vec =
            {
                0.0004270039294528318,
                -0.001228915451058514,
                -0.0001944260688400887,
                0.003338419923885592,
                -0.0007958936865080434,
                -0.004338208690251972,
                0.05379143536399187,
                0.006211299937499411,
                0.03042903097250921
            };
            double[] j1_vec =
            {
                1.0,
                1.5,
                2.0,
                1.0,
                1.5,
                2.0,
                0.5,
                1.0,
                1.5
            };
            double[] j2_vec =
            {
                8.0,
                8.0,
                8.0,
                3.0,
                3.0,
                3.0,
                0.5,
                0.5,
                0.5
            };
            double[] j3_vec =
            {
                7.0,
                7.0,
                7.0,
                2.0,
                2.0,
                2.0,
                1.0,
                1.0,
                1.0
            };
            double[] j4_vec =
            {
                6.5,
                6.5,
                6.5,
                4.0,
                4.0,
                4.0,
                2.0,
                2.0,
                2.0
            };
            double[] j5_vec =
            {
                7.5,
                7.5,
                7.5,
                1.5,
                1.5,
                1.5,
                1.0,
                1.0,
                1.0
            };
            double[] j6_vec =
            {
                7.5,
                7.5,
                7.5,
                3.0,
                3.0,
                3.0,
                1.5,
                1.5,
                1.5
            };
            double[] j7_vec =
            {
                6.0,
                6.0,
                6.0,
                3.5,
                3.5,
                3.5,
                1.5,
                1.5,
                1.5
            };
            double[] j8_vec =
            {
                10.0,
                10.0,
                10.0,
                2.0,
                2.0,
                2.0,
                0.5,
                0.5,
                0.5
            };
            double[] j9_vec =
            {
                6.0,
                6.0,
                6.0,
                2.0,
                2.0,
                2.0,
                1.5,
                1.5,
                1.5
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                j1 = 0.0;
                j2 = 0.0;
                j3 = 0.0;
                j4 = 0.0;
                j5 = 0.0;
                j6 = 0.0;
                j7 = 0.0;
                j8 = 0.0;
                j9 = 0.0;
                fx = 0.0;
            }
            else
            {
                j1 = j1_vec[n_data - 1];
                j2 = j2_vec[n_data - 1];
                j3 = j3_vec[n_data - 1];
                j4 = j4_vec[n_data - 1];
                j5 = j5_vec[n_data - 1];
                j6 = j6_vec[n_data - 1];
                j7 = j7_vec[n_data - 1];
                j8 = j8_vec[n_data - 1];
                j9 = j9_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void six_j_values(ref int n_data, ref double j1, ref double j2, ref double j3,
                ref double j4, ref double j5, ref double j6, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIX_J_VALUES returns some values of the Wigner 6J function.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      SixJSymbol[{j1,j2,j3},{j4,j5,j6}]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 February 2007
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
            //    Output, ref double J1, &J2, &J3, &J4, &J5, &J6, the arguments
            //    of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 15;

            double[] fx_vec =
            {
                0.03490905138373300,
                -0.03743025039659792,
                0.01890866390959560,
                0.007342448254928643,
                -0.02358935185081794,
                0.01913476955215437,
                0.001288017397724172,
                -0.01930018366290527,
                0.01677305949382889,
                0.005501147274850949,
                -0.02135439790896831,
                0.003460364451435387,
                0.02520950054795585,
                0.01483990561221713,
                0.002708577680633186
            };
            double[] j1_vec =
            {
                1.0,
                2.0,
                3.0,
                4.0,
                5.0,
                6.0,
                7.0,
                8.0,
                9.0,
                10.0,
                11.0,
                12.0,
                13.0,
                14.0,
                15.0
            };
            double[] j2_vec =
            {
                8.0,
                8.0,
                8.0,
                8.0,
                8.0,
                8.0,
                8.0,
                8.0,
                8.0,
                8.0,
                8.0,
                8.0,
                8.0,
                8.0,
                8.0
            };
            double[] j3_vec =
            {
                7.0,
                7.0,
                7.0,
                7.0,
                7.0,
                7.0,
                7.0,
                7.0,
                7.0,
                7.0,
                7.0,
                7.0,
                7.0,
                7.0,
                7.0
            };
            double[] j4_vec =
            {
                6.5,
                6.5,
                6.5,
                6.5,
                6.5,
                6.5,
                6.5,
                6.5,
                6.5,
                6.5,
                6.5,
                6.5,
                6.5,
                6.5,
                6.5
            };
            double[] j5_vec =
            {
                7.5,
                7.5,
                7.5,
                7.5,
                7.5,
                7.5,
                7.5,
                7.5,
                7.5,
                7.5,
                7.5,
                7.5,
                7.5,
                7.5,
                7.5
            };
            double[] j6_vec =
            {
                7.5,
                7.5,
                7.5,
                7.5,
                7.5,
                7.5,
                7.5,
                7.5,
                7.5,
                7.5,
                7.5,
                7.5,
                7.5,
                7.5,
                7.5
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                j1 = 0.0;
                j2 = 0.0;
                j3 = 0.0;
                j4 = 0.0;
                j5 = 0.0;
                j6 = 0.0;
                fx = 0.0;
            }
            else
            {
                j1 = j1_vec[n_data - 1];
                j2 = j2_vec[n_data - 1];
                j3 = j3_vec[n_data - 1];
                j4 = j4_vec[n_data - 1];
                j5 = j5_vec[n_data - 1];
                j6 = j6_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void three_j_values(ref int n_data, ref double j1, ref double j2, ref double j3,
                ref double m1, ref double m2, ref double m3, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    THREE_J_VALUES returns some values of the Wigner 3J function.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      ThreeJSymbol[{j1,m1},{j2,m2},{j3,m3}]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 February 2007
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
            //    Output, ref double J1, &J2, &J3, &M1, &M2, &M3, the arguments
            //    of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 8;

            double[] fx_vec =
            {
                0.2788866755113585,
                -0.09534625892455923,
                -0.06741998624632421,
                0.1533110351679666,
                -0.1564465546936860,
                0.1099450412156551,
                -0.05536235693131719,
                0.01799835451137786
            };
            double[] j1_vec =
            {
                1.0,
                2.0,
                3.0,
                4.0,
                5.0,
                6.0,
                7.0,
                8.0
            };
            double[] j2_vec =
            {
                4.5,
                4.5,
                4.5,
                4.5,
                4.5,
                4.5,
                4.5,
                4.5
            };
            double[] j3_vec =
            {
                3.5,
                3.5,
                3.5,
                3.5,
                3.5,
                3.5,
                3.5,
                3.5
            };
            double[] m1_vec =
            {
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0
            };
            double[] m2_vec =
            {
                -3.5,
                -3.5,
                -3.5,
                -3.5,
                -3.5,
                -3.5,
                -3.5,
                -3.5
            };
            double[] m3_vec =
            {
                2.5,
                2.5,
                2.5,
                2.5,
                2.5,
                2.5,
                2.5,
                2.5
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                j1 = 0.0;
                j2 = 0.0;
                j3 = 0.0;
                m1 = 0.0;
                m2 = 0.0;
                m3 = 0.0;
                fx = 0.0;
            }
            else
            {
                j1 = j1_vec[n_data - 1];
                j2 = j2_vec[n_data - 1];
                j3 = j3_vec[n_data - 1];
                m1 = m1_vec[n_data - 1];
                m2 = m2_vec[n_data - 1];
                m3 = m3_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

    }
}