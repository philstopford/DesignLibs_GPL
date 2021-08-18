namespace Burkardt.Values
{
    public static class Phi
    {

        public static void phi_values(ref int n_data, ref int n, ref int c)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PHI_VALUES returns some values of the PHI function.
            //
            //  Discussion:
            //
            //    PHI(N) is the number of integers between 1 and N which are
            //    relatively prime to N.  I and J are relatively prime if they
            //    have no common factors.  The function PHI(N) is known as
            //    "Euler's totient function".
            //
            //    By convention, 1 and N are relatively prime.
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      EulerPhi[n]
            //
            //  First values:
            //
            //     N  PHI(N)
            //
            //     1    1
            //     2    1
            //     3    2
            //     4    2
            //     5    4
            //     6    2
            //     7    6
            //     8    4
            //     9    6
            //    10    4
            //    11   10
            //    12    4
            //    13   12
            //    14    6
            //    15    8
            //    16    8
            //    17   16
            //    18    6
            //    19   18
            //    20    8
            //
            //  Formula:
            //
            //    PHI(U*V) = PHI(U) * PHI(V) if U and V are relatively prime.
            //
            //    PHI(P^K) = P^(K-1) * ( P - 1 ) if P is prime.
            //
            //    PHI(N) = N * Product ( P divides N ) ( 1 - 1 / P )
            //
            //    N = Sum ( D divides N ) PHI(D).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    09 February 2003
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
            //    Output, ref int N, the argument of the PHI function.
            //
            //    Output, ref int C, the value of the PHI function.
            //
        {
            int N_MAX = 20;

            int[] c_vec =
            {
                1, 1, 2, 2, 4, 2, 6, 4, 6, 4,
                8, 8, 16, 20, 16, 40, 148, 200, 200, 648
            };

            int[] n_vec =
            {
                1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                20, 30, 40, 50, 60, 100, 149, 500, 750, 999
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