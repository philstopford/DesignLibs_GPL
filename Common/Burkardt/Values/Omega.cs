namespace Burkardt.Values
{
    public static class Omega
    {
        public static void omega_values(ref int n_data, ref int n, ref int c)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    OMEGA_VALUES returns some values of the OMEGA function.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by
            //
            //      Length [ FactorInteger [ n ] ]
            //
            //  First values:
            //
            //     N   OMEGA(N)
            //
            //     1    0
            //     2    1
            //     3    1
            //     4    1
            //     5    1
            //     6    2
            //     7    1
            //     8    1
            //     9    1
            //    10    2
            //    11    1
            //    12    2
            //    13    1
            //    14    2
            //    15    2
            //    16    1
            //    17    1
            //    18    2
            //    19    1
            //    20    2
            //
            //  Formula:
            //
            //    If N = 1, then
            //
            //      OMEGA(N) = 0
            //
            //    else if the prime factorization of N is
            //
            //      N = P1^E1 * P2^E2 * ... * PM^EM,
            //
            //    then
            //
            //      OMEGA(N) = M
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 April 2013
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
            //    Output, ref int N, the argument of the OMEGA function.
            //
            //    Output, ref int C, the value of the OMEGA function.
            //
        {
            int N_MAX = 23;

            int[] c_vec =
            {
                0, 1, 1, 1, 1,
                2, 1, 1, 1, 2,
                3, 1, 4, 4, 3,
                1, 5, 2, 2, 1,
                6, 7, 8
            };

            int[] n_vec =
            {
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
                30,
                101,
                210,
                1320,
                1764,
                2003,
                2310,
                2827,
                8717,
                12553,
                30030,
                510510,
                9699690
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

    }
}