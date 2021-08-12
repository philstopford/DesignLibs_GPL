namespace Burkardt
{
    public static class Schroeder
    {
        public static void schroeder(int n, ref int[] s)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SCHROEDER generates the Schroeder numbers.
            //
            //  Discussion:
            //
            //    The Schroeder number S(N) counts the number of ways to insert
            //    parentheses into an expression of N items, with two or more items within
            //    a parenthesis.
            //
            //    Note that the Catalan number C(N) counts the number of ways
            //    to legally arrange a set of N left and N right parentheses.
            //
            //    The formula is:
            //
            //      S(N) = ( P(N)(3.0) - 3 P(N-1)(3.0) ) / ( 4 * ( N - 1 ) )
            //
            //    where P(N)(X) is the N-th Legendre polynomial.
            //
            //  Example:
            //
            //    N = 4
            //
            //    1234
            //    12(34)
            //    1(234)
            //    1(2(34))
            //    1(23)4
            //    1((23)4)
            //    (123)4
            //    (12)34
            //    (12)(34)
            //    (1(23))4
            //    ((12)3)4
            //
            //  First Values:
            //
            //           1
            //           1
            //           3
            //          11
            //          45
            //         197
            //         903
            //        4279
            //       20793
            //      103049
            //      518859
            //     2646723
            //    13648869
            //    71039373
            //
            //  Recursion:
            //
            //    S(1) = 1
            //    S(2) = 1
            //    S(N) = ( ( 6 * N - 9 ) * S(N-1) - ( N - 3 ) * S(N-2) ) / N
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 May 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    R P Stanley,
            //    Hipparchus, Plutarch, Schroeder, and Hough,
            //    American Mathematical Monthly,
            //    Volume 104, Number 4, 1997, pages 344-350.
            //
            //    Laurent Habsieger, Maxim Kazarian, Sergei Lando,
            //    On the Second Number of Plutarch,
            //    American Mathematical Monthly, May 1998, page 446.
            //
            //  Parameters:
            //
            //    Input, int N, the number of Schroeder numbers desired.
            //
            //    Output, int S[N], the Schroeder numbers.
            //
        {
            int i;

            if (n <= 0)
            {
                return;
            }

            s[0] = 1;

            if (n <= 1)
            {
                return;
            }

            s[1] = 1;

            if (n <= 2)
            {
                return;
            }

            for (i = 3; i <= n; i++)
            {
                s[i - 1] = ((6 * i - 9) * s[i - 2] - (i - 3) * s[i - 3]) / i;
            }
        }
    }
}