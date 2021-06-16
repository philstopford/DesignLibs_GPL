﻿namespace Burkardt.RankingNS
{
    public static partial class Ranking
    {
        public static void gamma_log_values(ref int n_data, ref double x, ref double fx )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GAMMA_LOG_VALUES returns some values of the Log Gamma function.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Log[Gamma[x]]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 August 2004
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
        //    Output, double &X, the argument of the function.
        //
        //    Output, double &FX, the value of the function.
        //
        {
            int N_MAX = 20;

            double[] fx_vec =
            {
                0.1524063822430784E+01,
                0.7966778177017837E+00,
                0.3982338580692348E+00,
                0.1520596783998375E+00,
                0.0000000000000000E+00,
                -0.4987244125983972E-01,
                -0.8537409000331584E-01,
                -0.1081748095078604E+00,
                -0.1196129141723712E+00,
                -0.1207822376352452E+00,
                -0.1125917656967557E+00,
                -0.9580769740706586E-01,
                -0.7108387291437216E-01,
                -0.3898427592308333E-01,
                0.00000000000000000E+00,
                0.69314718055994530E+00,
                0.17917594692280550E+01,
                0.12801827480081469E+02,
                0.39339884187199494E+02,
                0.71257038967168009E+02
            }
            ;

            double[] x_vec =
            {
                0.20E+00,
                0.40E+00,
                0.60E+00,
                0.80E+00,
                1.00E+00,
                1.10E+00,
                1.20E+00,
                1.30E+00,
                1.40E+00,
                1.50E+00,
                1.60E+00,
                1.70E+00,
                1.80E+00,
                1.90E+00,
                2.00E+00,
                3.00E+00,
                4.00E+00,
                10.00E+00,
                20.00E+00,
                30.00E+00
            }
            ;

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