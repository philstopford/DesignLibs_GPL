namespace Burkardt.Sequence;

public static class Collatz
{
    public static int collatz_count(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COLLATZ_COUNT counts the number of terms in a Collatz sequence.
        //
        //  Discussion:
        //
        //    The rules for generation of the Collatz sequence are recursive.
        //    If T is the current entry of the sequence, (T is
        //    assumed to be a positive integer), then the next
        //    entry, U is determined as follows:
        //
        //      if T is 1 (or less)
        //        terminate the sequence;
        //      else if T is even
        //        U = T/2.
        //      else (if T is odd and not 1)
        //        U = 3*T+1;
        //
        //     N  Sequence                                                Length
        //
        //     1                                                               1
        //     2   1                                                           2
        //     3  10,  5, 16,  8,  4,  2,  1                                   8
        //     4   2   1                                                       3
        //     5  16,  8,  4,  2,  1                                           6
        //     6   3, 10,  5, 16,  8,  4,  2,  1                               9
        //     7  22, 11, 34, 17, 52, 26, 13, 40, 20, 10, 5, 16, 8, 4, 2, 1   17
        //     8   4,  2,  1                                                   4
        //     9  28, 14,  7, ...                                             20
        //    10   5, 16,  8,  4,  2,  1                                       7
        //    11  34, 17, 52, 26, 13, 40, 20, 10,  5, 16, 8, 4, 2, 1          15
        //    12   6,  3, 10,  5, 16,  8,  4,  2,  1                          10
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 March 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Eric Weisstein,
        //    CRC Concise Encyclopedia of Mathematics,
        //    CRC Press, 2002,
        //    Second edition,
        //    ISBN: 1584883472,
        //    LC: QA5.W45.
        //
        //  Parameters:
        //
        //    Input, int N, the first element of the sequence.
        //
        //    Output, int COLLATZ_COUNT, the number of elements in
        //    the Collatz sequence that begins with N.
        //
    {
        int count = 1;

        for (;;)
        {
            if (n <= 1)
            {
                break;
            }

            switch (n % 2)
            {
                case 0:
                    n /= 2;
                    break;
                default:
                    n = 3 * n + 1;
                    break;
            }

            count += 1;
        }

        return count;
    }

    public static void collatz_count_max(int n, ref int i_max, ref int j_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COLLATZ_COUNT_MAX seeks the maximum Collatz count for 1 through N.
        //
        //  Discussion:
        //
        //    For each integer I, we compute a sequence of values that 
        //    terminate when we reach 1.  The number of steps required to
        //    reach 1 is the "rank" of I, and we are searching the numbers
        //    from 1 to N for the number with maximum rank.
        //
        //    For a given I, the sequence is produced by:
        //
        //    1) J = 1, X(J) = I;
        //    2) If X(J) = 1, stop.
        //    3) J = J + 1; 
        //       if X(J-1) was even, X(J) = X(J-1)/2;
        //       else                X(J) = 3 * X(J-1) + 1;
        //    4) Go to 3
        //
        //  Example:
        //
        //            N      I_MAX J_MAX
        //
        //           10          9    20
        //          100         97   119
        //        1,000        871   179
        //       10,000      6,171   262
        //      100,000     77,031   351
        //    1,000,000    837,799   525
        //   10,000,000  8,400,511   686
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 April 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the maximum integer to check.
        //
        //    Output, int *I_MAX, *J_MAX, an integer I with the maximum rank,
        //    and the value of the maximum rank.
        //
    {
        int i;

        i_max = -1;
        j_max = -1;

        for (i = 1; i <= n; i++)
        {
            int j = 1;
            int x = i;

            while (x != 1)
            {
                j += 1;
                switch (x % 2)
                {
                    case 0:
                        x /= 2;
                        break;
                    default:
                        x = 3 * x + 1;
                        break;
                }
            }

            if (j_max >= j)
            {
                continue;
            }

            i_max = i;
            j_max = j;
        }
    }
}