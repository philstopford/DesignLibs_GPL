using Burkardt.FullertonFnLib;
using Burkardt.Types;

namespace Burkardt
{
    public static class Comb
    {
        public static int[] comb ( int n, int p, int l )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    COMB selects a subset of order P from a set of order N.
            //
            //  Discussion:
            //
            //    This subroutine finds the combination set of N things taken
            //    P at a time for a given lexicographic index.
            //
            //  Modified:
            //
            //    01 April 2016
            //
            //  Author:
            //
            //    Bill Buckles, Matthew Lybanon
            //
            //  Reference:
            //
            //    Bill Buckles, Matthew Lybanon,
            //    Algorithm 515: Generation of a Vector from the Lexicographical Index,
            //    ACM Transactions on Mathematical Software,
            //    Volume 3, Number 2, June 1977, pages 180-182.
            //
            //  Parameters:
            //
            //    Input, int N, the number of things in the set.
            //
            //    Input, int P, the number of things in each combination.
            //    0 < P < N.
            //
            //    Input, int L, the lexicographic index of the 
            //    desired combination.  1 <= L <= choose(N,P).
            //
            //    Output, int COMB[P], the combination set.
            //
        {
            int[] c;
            int i;
            int k;
            int p1;
            int r;

            c = new int[p];
            //
            //  Special case: P = 1
            //
            if ( p == 1 )
            {
                c[0] = l;
                return c;
            }
            //
            //  Initialize lower bound index.
            //
            k = 0;
            //
            //  Select elements in ascending order.
            //
            p1 = p - 1;
            c[0] = 0;

            for ( i = 1; i <= p1; i++ )
            {
                //
                //  Update lower bound as the previously selected element.
                //
                if ( 1 < i )
                {
                    c[i-1] = c[i-2];
                }
                //
                //  Check validity of each entry.
                //
                for ( ; ; )
                {
                    c[i-1] = c[i-1] + 1;
                    r = FullertonLib.i4_binom ( n - c[i-1], p - i );
                    k = k + r;

                    if ( l <= k )
                    {
                        break;
                    }
                }
                k = k - r;
            }

            c[p-1] = c[p1-1] + l - k;

            return c;
        }

        public static void comb_next(int n, int k, int[] a, ref bool done)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    COMB_NEXT computes combinations of K things out of N.
            //
            //  Discussion:
            //
            //    The combinations are computed one at a time, in lexicographical order.
            //
            //    10 April 1009: Thanks to "edA-qa mort-ora-y" for supplying a 
            //    correction to this code!
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 November 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Charles Mifsud,
            //    Combination in Lexicographic Order,
            //    ACM algorithm 154,
            //    Communications of the ACM,
            //    March 1963.
            //
            //  Parameters:
            //
            //    Input, int N, the total number of things.
            //
            //    Input, int K, the number of things in each combination.
            //
            //    Input/output, int A[K], contains the list of elements in
            //    the current combination.
            //
            //    Input/output, bool &DONE.  Set DONE to TRUE before the first call,
            //    and then use the output value from the previous call on subsequent
            //    calls.  The output value will be FALSE as long as there are more
            //    combinations to compute, and TRUE when the list is exhausted.
            //
        {
            int i;
            int j;

            if (done)
            {
                if (k <= 0)
                {
                    return;
                }

                typeMethods.i4vec_indicator1(k, ref a);

                done = false;
            }
            else
            {
                if (a[k - 1] < n)
                {
                    a[k - 1] = a[k - 1] + 1;
                    return;
                }

                for (i = k; 2 <= i; i--)
                {
                    if (a[i - 2] < n - k + i - 1)
                    {
                        a[i - 2] = a[i - 2] + 1;

                        for (j = i; j <= k; j++)
                        {
                            a[j - 1] = a[i - 2] + j - (i - 1);
                        }

                        return;
                    }
                }

                done = true;
            }

        }

        public static void comb_row_next(int n, ref int[] row)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    COMB_ROW_NEXT computes the next row of Pascal's triangle.
            //
            //  Discussion:
            //
            //    Row N contains the N+1 combinatorial coefficients
            //
            //      C(N,0), C(N,1), C(N,2), ... C(N,N)
            //
            //    The sum of the elements of row N is equal to 2**N.
            //
            //    The formula is:
            //
            //      C(N,K) = N! / ( K! * (N-K)! )
            //
            //  First terms:
            //
            //     N K:0  1   2   3   4   5   6   7  8  9 10
            //
            //     0   1
            //     1   1  1
            //     2   1  2   1
            //     3   1  3   3   1
            //     4   1  4   6   4   1
            //     5   1  5  10  10   5   1
            //     6   1  6  15  20  15   6   1
            //     7   1  7  21  35  35  21   7   1
            //     8   1  8  28  56  70  56  28   8  1
            //     9   1  9  36  84 126 126  84  36  9  1
            //    10   1 10  45 120 210 252 210 120 45 10  1
            //
            //  Recursion:
            //
            //    C(N,K) = C(N-1,K-1)+C(N-1,K)
            //
            //  Special values:
            //
            //    C(N,0) = C(N,N) = 1
            //    C(N,1) = C(N,N-1) = N
            //    C(N,N-2) = sum ( 1 <= I <= N ) N
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 December 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the row of the triangle desired.
            //
            //    Input/output, int ROW[N+1].  On input, row N-1 is
            //    contained in entries 0 through N-1.  On output, row N is contained
            //    in entries 0 through N.
            //
        {
            int i;

            if (n < 0)
            {
                return;
            }

            row[0] = 1;
            for (i = n - 1; 1 <= i; i--)
            {
                row[i] = row[i] + row[i - 1];
            }

            row[n] = 1;
        }

        public static void comb_unrank(int m, int n, int rank, ref int[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    COMB_UNRANK returns the RANK-th combination of N things out of M.
            //
            //  Discussion:
            //
            //    The combinations are ordered lexically.
            //
            //    Lexical order can be illustrated for the general case of N and M as
            //    follows:
            //
            //    1:       1,     2,     3,     ..., N-2, N-1, N
            //    2:       1,     2,     3,     ..., N-2, N-1, N+1
            //    3:       1,     2,     3,     ..., N-2, N-1, N+2
            //    ...
            //    M-N+1:   1,     2,     3,     ..., N-2, N-1, M
            //    M-N+2:   1,     2,     3,     ..., N-2, N,   N+1
            //    M-N+3:   1,     2,     3,     ..., N-2, N,   N+2
            //    ...
            //    LAST-2:  M-N,   M-N+1, M-N+3, ..., M-2, M-1, M
            //    LAST-1:  M-N,   M-N+2, M-N+3, ..., M-2, M-1, M
            //    LAST:    M-N+1, M-N+2, M-N+3, ..., M-2, M-1, M
            //
            //    There are a total of M!/(N!*(M-N)!) combinations of M
            //    things taken N at a time.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 June 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    B P Buckles, M Lybanon,
            //    Algorithm 515,
            //    Generation of a Vector from the Lexicographical Index,
            //    ACM Transactions on Mathematical Software,
            //    Volume 3, Number 2, pages 180-182, June 1977.
            //
            //  Parameters:
            //
            //    Input, int M, the size of the set.
            //
            //    Input, int N, the number of things in the combination.
            //    N must be greater than 0, and no greater than M.
            //
            //    Input, int RANK, the lexicographical index of combination
            //    sought.  RANK must be at least 1, and no greater than M!/(N!*(M-N)!).
            //
            //    Output, int A[N], array containing the combination set.
            //
        {
            int i;
            int j;
            int k;
            //
            //  Initialize lower bound index at zero.
            //
            k = 0;
            //
            //  Loop to select elements in ascending order.
            //
            for (i = 1; i <= n - 1; i++)
            {
                //
                //  Set lower bound element number for next element value.
                //
                a[i - 1] = 0;

                if (1 < i)
                {
                    a[i - 1] = a[i - 2];
                }

                //
                //  Check each element value.
                //
                for (;;)
                {
                    a[i - 1] = a[i - 1] + 1;
                    j = typeMethods.i4_choose(m - a[i - 1], n - i);
                    k = k + j;

                    if (rank <= k)
                    {
                        break;
                    }
                }

                k = k - j;
            }

            a[n - 1] = a[n - 2] + rank - k;

            return;
        }

    }
}