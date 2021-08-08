namespace TestValues
{
    public static class Tau
    {
        public static void tau_values(ref int n_data, ref int n, ref int c)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TAU_VALUES returns some values of the Tau function.
            //
            //  Discussion:
            //
            //    TAU(N) is the number of divisors of N, including 1 and N.
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      DivisorSigma[1,n]
            //
            //  First values:
            //
            //     N   TAU(N)
            //
            //     1    1
            //     2    2
            //     3    2
            //     4    3
            //     5    2
            //     6    4
            //     7    2
            //     8    4
            //     9    3
            //    10    4
            //    11    2
            //    12    6
            //    13    2
            //    14    4
            //    15    4
            //    16    5
            //    17    2
            //    18    6
            //    19    2
            //    20    6
            //
            //  Formula:
            //
            //    If the prime factorization of N is
            //
            //      N = P1**E1 * P2**E2 * ... * PM**EM,
            //
            //    then
            //
            //      TAU(N) = ( E1 + 1 ) * ( E2 + 1 ) * ... * ( EM + 1 ).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    10 February 2003
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
            //    Output, ref int N, the argument of the Tau function.
            //
            //    Output, ref int C, the value of the Tau function.
            //
        {
            int N_MAX = 20;

            int[] c_vec =
            {
                1, 2, 2, 3, 2, 4, 2, 4, 3, 4,
                2, 12, 12, 4, 18, 24, 2, 8, 14, 28
            };

            int[] n_vec =
            {
                1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                23, 72, 126, 226, 300, 480, 521, 610, 832, 960
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