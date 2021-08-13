﻿namespace Burkardt.Sequence
{
    public static class Euler
    {
        public static void eulerian ( int n, ref int[] e )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    EULERIAN computes the Eulerian number E(N,K).
            //
            //  Definition:
            //
            //    A run in a permutation is a sequence of consecutive ascending values.
            //
            //    E(N,K) is the number of permutations of N objects which contain
            //    exactly K runs.
            //
            //  Examples:
            //
            //     N = 7
            //
            //     1     0     0     0     0     0     0
            //     1     1     0     0     0     0     0
            //     1     4     1     0     0     0     0
            //     1    11    11     1     0     0     0
            //     1    26    66    26     1     0     0
            //     1    57   302   302    57     1     0
            //     1   120  1191  2416  1191   120     1
            //
            //  Recursion:
            //
            //    E(N,K) = K * E(N-1,K) + (N-K+1) * E(N-1,K-1).
            //
            //  Properties:
            //
            //    E(N,1) = E(N,N) = 1.
            //    E(N,K) = 0 if K <= 0 or N < K.
            //    sum ( 1 <= K <= N ) E(N,K) = N!.
            //    X^N = sum ( 0 <= K <= N ) COMB(X+K-1, N ) E(N,K)
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 May 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Dennis Stanton and Dennis White,
            //    Constructive Combinatorics,
            //    Springer Verlag, 1986
            //
            //  Parameters:
            //
            //    Input, int N, the number of rows desired.
            //
            //    Output, int E[N*N], the first N rows of Eulerian numbers.
            //
        {
            int i;
            int j;

            if ( n < 1 )
            {
                return;
            }
            //
            //  Construct rows 1, 2, ..., N of the Eulerian triangle.
            //
            e[1-1+(1-1)*n] = 1;
            for ( j = 2; j <= n; j++ )
            {
                e[1-1+(j-1)*n] = 0;
            }

            for ( i = 2; i <= n; i++ )
            {
                e[i-1+(1-1)*n] = 1;
                for ( j = 2; j <= n; j++ )
                {
                    e[i-1+(j-1)*n] = j * e[i-2+(j-1)*n] + ( i - j + 1 ) * e[i-2+(j-2)*n];
                }
            }
        }
    }
}