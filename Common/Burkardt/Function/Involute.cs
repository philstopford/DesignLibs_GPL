namespace Burkardt.Function
{
    public static class Involute
    {
        public static void involute_enum ( int n, ref int[] s )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    INVOLUTE_ENUM enumerates the involutions of N objects.
            //
            //  Discussion:
            //
            //    An involution is a permutation consisting only of fixed points and
            //    pairwise transpositions.
            //
            //    An involution is its own inverse permutation.
            //
            //  Recursion:
            //
            //    S(0) = 1
            //    S(1) = 1
            //    S(N) = S(N-1) + (N-1) * S(N-2)
            //
            //  First values:
            //
            //     N         S(N)
            //     0           1
            //     1           1
            //     2           2
            //     3           4
            //     4          10
            //     5          26
            //     6          76
            //     7         232
            //     8         764
            //     9        2620
            //    10        9496
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    28 May 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of objects to be permuted.
            //
            //    Output, int S[N+1], the number of involutions of 0, 1, 2, ... N
            //    objects.
            //
        {
            int i;

            if ( n < 0 )
            {
                return;
            }

            s[0] = 1;

            if ( n <= 0 )
            {
                return;
            }

            s[1] = 1;

            for ( i = 2; i <= n; i++ )
            {
                s[i] = s[i-1] + ( i - 1 ) * s[i-2];
            }
        }
    }
}