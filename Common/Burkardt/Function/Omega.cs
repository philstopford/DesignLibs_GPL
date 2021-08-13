using System;
using Burkardt.Types;

namespace Burkardt.Function
{
    public static class Omega
    {
        public static int omega ( int n )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    OMEGA returns OMEGA(N), the number of distinct prime divisors of N.
            //
            //  Discussion:
            //
            //    If N = 1, then
            //
            //      OMEGA(N) = 1
            //
            //    else if the prime factorization of N is
            //
            //      N = P1**E1 * P2^E2 * ... * PM^EM,
            //
            //    then
            //
            //      OMEGA(N) = M
            //
            //  First values:
            //
            //     N   OMEGA(N)
            //
            //     1    1
            //     2    1
            //     3    1
            //     4    1
            //     5    1
            //     6    2
            //     7    1
            //     8    1
            //     9    1
            //    10    2
            //    11    1
            //    12    2
            //    13    1
            //    14    2
            //    15    2
            //    16    1
            //    17    1
            //    18    2
            //    19    1
            //    20    2
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    18 November 2000
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the value to be analyzed.  N must be 1 or
            //    greater.
            //
            //    Output, int OMEGA, the value of OMEGA(N).  But if N is 0 or
            //    less, OMEGA is returned as 0, a nonsense value.  If there is
            //    not enough room for factoring, OMEGA is returned as -1.
            //
        {
            int FACTOR_MAX = 20;

            int[] factor = new int[FACTOR_MAX];
            int nfactor = 0;
            int nleft = 0;
            int[] power = new int[FACTOR_MAX];

            if ( n <= 0 )
            {
                return 0;
            }

            if ( n == 1 )
            {
                return 1;
            }
            //
            //  Factor N.
            //
            typeMethods.i4_factor ( n, FACTOR_MAX, ref nfactor, ref factor, ref power, ref nleft );

            if ( nleft != 1 )
            {
                Console.WriteLine("");
                Console.WriteLine("OMEGA - Fatal error!");
                Console.WriteLine("  Not enough factorization space.");
                return ( 1 );
            }

            return nfactor;
        }

    }
}