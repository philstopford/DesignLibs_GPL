namespace Burkardt.Sequence
{
    public class Delannoy
    {
        public static int[] delannoy(int m, int n)

            //****************************************************************************0
            //
            //  Purpose:
            //
            //    DELANNOY returns the Delannoy numbers up to orders (M,N).
            //
            //  Discussion:
            //
            //    The Delannoy number A(M,N) counts the number of distinct paths
            //    from (0,0) to (M,N) in which the only steps used are
            //    (1,1), (1,0) and (0,1).
            //
            //  First values:
            //
            //      \N 0  1   2   3    4     5     6      7      8
            //     M-+--------------------------------------------
            //     0 | 1  1   1   1    1     1     1      1      1
            //     1 | 1  3   5   7    9    11    13     15     17
            //     2 | 1  5  13  25   41    61    85    113    145
            //     3 | 1  7  25  63  129   231   377    575    833
            //     4 | 1  9  41 129  321   681  1289   2241   3649
            //     5 | 1 11  61 231  681  1683  3653   7183  13073
            //     6 | 1 13  85 377 1289  3653  8989  19825  40081
            //     7 | 1 15 113 575 2241  7183 19825  48639 108545
            //     8 | 1 17 145 833 3649 13073 40081 108545 265729
            //
            //  Recursion:
            //
            //    A(0,0) = 1
            //    A(M,0) = 1
            //    A(0,N) = 1
            //    A(M,N) = A(M-1,N) + A(M,N-1) + A(M-1,N-1)
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 June 2018
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
            //    LC: QA5.W45
            //
            //  Parameters:
            //
            //    Input, int M, N, define the highest order number to 
            //    compute.
            //
            //    Output, int DELANNOY[(M+1)*(N+1)], the Delannoy numbers.
            //
        {
            int[] a;
            int i;
            int j;

            if (m < 0)
            {
                return null;
            }

            if (n < 0)
            {
                return null;
            }

            a = new int[(m + 1) * (n + 1)];

            a[0 + 0 * (m + 1)] = 1;

            for (i = 1; i <= m; i++)
            {
                a[i + 0 * (m + 1)] = 1;
            }

            for (j = 1; j <= n; j++)
            {
                a[0 + j * (m + 1)] = 1;
            }

            for (i = 1; i <= m; i++)
            {
                for (j = 1; j <= n; j++)
                {
                    a[i + j * (m + 1)] =
                        a[i - 1 + j * (m + 1)] + a[i + (j - 1) * (m + 1)] + a[i - 1 + (j - 1) * (m + 1)];
                }
            }

            return a;
        }
    }
}