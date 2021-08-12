namespace Burkardt.TestValues
{
    public static class Synch
    {
        public static void synch1_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SYNCH1_VALUES returns some values of the synchrotron radiation function.
            //
            //  Discussion:
            //
            //    The function is defined by:
            //
            //      SYNCH1(x) = x * integral ( x <= t < +oo ) K(5/3)(t) dt
            //
            //    where K(5/3) is a modified Bessel function of order 5/3.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 September 2004
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
            //    Allan McLeod,
            //    Algorithm 757:
            //    MISCFUN: A software package to compute uncommon special functions,
            //    ACM Transactions on Mathematical Software,
            //    Volume 22, Number 3, September 1996, pages 288-301.
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
            int N_MAX = 20;

            double[] fx_vec =
            {
                0.26514864547487397044E+00,
                0.62050129979079045645E+00,
                0.85112572132368011206E+00,
                0.87081914687546885094E+00,
                0.65142281535536396975E+00,
                0.45064040920322354579E+00,
                0.30163590285073940285E+00,
                0.19814490804441305867E+00,
                0.12856571000906381300E+00,
                0.52827396697866818297E-01,
                0.42139298471720305542E-01,
                0.21248129774981984268E-01,
                0.13400258907505536491E-01,
                0.84260797314108699935E-02,
                0.12884516186754671469E-02,
                0.19223826430086897418E-03,
                0.28221070834007689394E-04,
                0.15548757973038189372E-05,
                0.11968634456097453636E-07,
                0.89564246772237127742E-10
            };

            double[] x_vec =
            {
                0.0019531250E+00,
                0.0312500000E+00,
                0.1250000000E+00,
                0.5000000000E+00,
                1.0000000000E+00,
                1.5000000000E+00,
                2.0000000000E+00,
                2.5000000000E+00,
                3.0000000000E+00,
                4.0000000000E+00,
                4.2500000000E+00,
                5.0000000000E+00,
                5.5000000000E+00,
                6.0000000000E+00,
                8.0000000000E+00,
                10.0000000000E+00,
                12.0000000000E+00,
                15.0000000000E+00,
                20.0000000000E+00,
                25.0000000000E+00
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

        public static void synch2_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SYNCH2_VALUES returns some values of the synchrotron radiation function.
            //
            //  Discussion:
            //
            //    The function is defined by:
            //
            //      SYNCH2(x) = x * K(2/3)(x)
            //
            //    where K(2/3) is a modified Bessel function of order 2/3.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 September 2004
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
            //    Allan McLeod,
            //    Algorithm 757:
            //    MISCFUN: A software package to compute uncommon special functions,
            //    ACM Transactions on Mathematical Software,
            //    Volume 22, Number 3, September 1996, pages 288-301.
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
            int N_MAX = 20;

            double[] fx_vec =
            {
                0.13430727275667378338E+00,
                0.33485265272424176976E+00,
                0.50404224110911078651E+00,
                0.60296523236016785113E+00,
                0.49447506210420826699E+00,
                0.36036067860473360389E+00,
                0.24967785497625662113E+00,
                0.16813830542905833533E+00,
                0.11117122348556549832E+00,
                0.46923205826101330711E-01,
                0.37624545861980001482E-01,
                0.19222123172484106436E-01,
                0.12209535343654701398E-01,
                0.77249644268525771866E-02,
                0.12029044213679269639E-02,
                0.18161187569530204281E-03,
                0.26884338006629353506E-04,
                0.14942212731345828759E-05,
                0.11607696854385161390E-07,
                0.87362343746221526073E-10
            };

            double[] x_vec =
            {
                0.0019531250E+00,
                0.0312500000E+00,
                0.1250000000E+00,
                0.5000000000E+00,
                1.0000000000E+00,
                1.5000000000E+00,
                2.0000000000E+00,
                2.5000000000E+00,
                3.0000000000E+00,
                4.0000000000E+00,
                4.2500000000E+00,
                5.0000000000E+00,
                5.5000000000E+00,
                6.0000000000E+00,
                8.0000000000E+00,
                10.0000000000E+00,
                12.0000000000E+00,
                15.0000000000E+00,
                20.0000000000E+00,
                25.0000000000E+00
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