using Burkardt.FullertonFnLib;

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
    }
}