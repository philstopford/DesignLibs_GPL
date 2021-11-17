namespace Burkardt.Values;

public static class Moebius
{

    public static void moebius_values(ref int n_data, ref int n, ref int c)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MOEBIUS_VALUES returns some values of the Moebius function.
        //
        //  Discussion:
        //
        //    MU(N) is defined as follows:
        //
        //      MU(N) = 1 if N = 1;
        //              0 if N is divisible by the square of a prime;
        //              (-1)^K, if N is the product of K distinct primes.
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      MoebiusMu[n]
        //
        //  First values:
        //
        //     N  MU(N)
        //
        //     1    1
        //     2   -1
        //     3   -1
        //     4    0
        //     5   -1
        //     6    1
        //     7   -1
        //     8    0
        //     9    0
        //    10    1
        //    11   -1
        //    12    0
        //    13   -1
        //    14    1
        //    15    1
        //    16    0
        //    17   -1
        //    18    0
        //    19   -1
        //    20    0
        //
        //  Note:
        //
        //    As special cases, MU(N) is -1 if N is a prime, and MU(N) is 0
        //    if N is a square, cube, etc.
        //
        //  Formula:
        //
        //    The Moebius function is related to Euler's totient function:
        //
        //      PHI(N) = Sum ( D divides N ) MU(D) * ( N / D ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 February 2003
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
        //    Output, ref int N, the argument of the Moebius function.
        //
        //    Output, ref int C, the value of the Moebius function.
        //
    {
        const int N_MAX = 20;

        int[] c_vec =
        {
            1, -1, -1, 0, -1, 1, -1, 0, 0, 1,
            -1, 0, -1, 1, 1, 0, -1, 0, -1, 0
        };

        int[] n_vec =
        {
            1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
            11, 12, 13, 14, 15, 16, 17, 18, 19, 20
        };

        n_data = n_data switch
        {
            < 0 => 0,
            _ => n_data
        };

        n_data += 1;

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