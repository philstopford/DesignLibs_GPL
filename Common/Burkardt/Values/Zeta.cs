namespace Burkardt.Values
{
    public class Zeta
    {
        public static void zeta_values(ref int n_data, ref int n, ref double zeta)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ZETA_VALUES returns some values of the Riemann Zeta function.
            //
            //  Discussion:
            //
            //    ZETA(N) = sum ( 1 <= I < +oo ) 1 / I^N
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Zeta[n]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 August 2004
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
            //    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, int &N, the argument of the Zeta function.
            //
            //    Output, double &ZETA, the value of the Zeta function.
            //
        {
            int N_MAX = 15;

            int[] n_vec =
            {
                2,
                3,
                4,
                5,
                6,
                7,
                8,
                9,
                10,
                11,
                12,
                16,
                20,
                30,
                40
            };

            double[] zeta_vec =
            {
                0.164493406684822643647E+01,
                0.120205690315959428540E+01,
                0.108232323371113819152E+01,
                0.103692775514336992633E+01,
                0.101734306198444913971E+01,
                0.100834927738192282684E+01,
                0.100407735619794433939E+01,
                0.100200839292608221442E+01,
                0.100099457512781808534E+01,
                0.100049418860411946456E+01,
                0.100024608655330804830E+01,
                0.100001528225940865187E+01,
                0.100000095396203387280E+01,
                0.100000000093132743242E+01,
                0.100000000000090949478E+01
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
                zeta = 0.0;
            }
            else
            {
                n = n_vec[n_data - 1];
                zeta = zeta_vec[n_data - 1];
            }
        }

        public static void zeta_m1_values(ref int n_data, ref double p, ref double zeta_m1)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ZETA_M1_VALUES returns some values of the Riemann Zeta Minus One function.
            //
            //  Discussion:
            //
            //    ZETA_M1(N) = ZETA(N) - 1
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    15 January 2017
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
            //    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, double &P, the argument.
            //
            //    Output, double &ZETA_M1, the value.
            //
        {
            int N_MAX = 17;

            double[] p_vec =
            {
                2.0,
                2.5,
                3.0,
                3.5,
                4.0,
                5.0,
                6.0,
                7.0,
                8.0,
                9.0,
                10.0,
                11.0,
                12.0,
                16.0,
                20.0,
                30.0,
                40.0
            };

            double[] zeta_m1_vec =
            {
                0.64493406684822643647E+00,
                0.3414872573E+00,
                0.20205690315959428540E+00,
                0.1267338673E+00,
                0.8232323371113819152E-01,
                0.3692775514336992633E-01,
                0.1734306198444913971E-01,
                0.834927738192282684E-02,
                0.407735619794433939E-02,
                0.200839292608221442E-02,
                0.99457512781808534E-03,
                0.49418860411946456E-03,
                0.24608655330804830E-03,
                0.1528225940865187E-04,
                0.95396203387280E-06,
                0.93132743242E-10,
                0.90949478E-12
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                p = 0.0;
                zeta_m1 = 0.0;
            }
            else
            {
                p = p_vec[n_data - 1];
                zeta_m1 = zeta_m1_vec[n_data - 1];
            }
        }
    }
}