namespace Burkardt.Sequence;

public static class Motzkin
{
    public static void motzkin ( int n, ref int[] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MOTZKIN returns the Motzkin numbers up to order N.
        //
        //  Discussion:
        //
        //    The Motzkin number A(N) counts the number of distinct paths
        //    from (0,0) to (0,N) in which the only steps used are
        //    (1,1), (1,-1) and (1,0), and the path is never allowed to
        //    go below the X axis.
        //
        //  First values:
        //
        //     N  A(N)
        //
        //     0    1
        //     1    1
        //     2    2
        //     3    4
        //     4    9
        //     5   21
        //     6   51
        //     7  127
        //     8  323
        //     9  835
        //    10 2188
        //
        //  Recursion:
        //
        //    A(N) = A(N-1) + sum ( 0 <= K <= N-2 ) A(K) * A(N-2-K)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 August 2004
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
        //    Input, int N, the highest order Motzkin number to compute.
        //
        //    Output, int A[N+1], the Motzkin numbers.
        //
    {
        int i;
        int j;

        switch (n)
        {
            case < 0:
                return;
        }

        a[0] = 1;

        for ( i = 1; i <= n; i++ )
        {
            a[i] = a[i-1];
            for ( j = 0; j <= i-2; j++ )
            {
                a[i] += a[j] * a[i-2-j];
            }
        }
    }

}