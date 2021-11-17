namespace Burkardt.Values;

public static class Sigma
{
    public static void sigma_values(ref int n_data, ref int n, ref int c)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SIGMA_VALUES returns some values of the Sigma function.
        //
        //  Discussion:
        //
        //    SIGMA(N) is the sum of the distinct divisors of N, including 1 and N.
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      DivisorSigma[1,n]
        //
        //  First values:
        //
        //     N  SIGMA(N)
        //
        //     1    1
        //     2    3
        //     3    4
        //     4    7
        //     5    6
        //     6   12
        //     7    8
        //     8   15
        //     9   13
        //    10   18
        //    11   12
        //    12   28
        //    13   14
        //    14   24
        //    15   24
        //    16   31
        //    17   18
        //    18   39
        //    19   20
        //    20   42
        //
        //  Formula:
        //
        //    SIGMA(U*V) = SIGMA(U) * SIGMA(V) if U and V are relatively prime.
        //
        //    SIGMA(P**K) = ( P**(K+1) - 1 ) / ( P - 1 ) if P is prime.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 February 2003
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
        //    Output, ref int N, the argument of the Sigma function.
        //
        //    Output, ref int C, the value of the Sigma function.
        //
    {
        const int N_MAX = 20;

        int[] c_vec =
        {
            1, 3, 4, 7, 6, 12, 8, 15, 13, 18,
            72, 128, 255, 176, 576, 1170, 618, 984, 2232, 2340
        };

        int[] n_vec =
        {
            1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
            30, 127, 128, 129, 210, 360, 617, 815, 816, 1000
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