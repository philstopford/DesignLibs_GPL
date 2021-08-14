using System;
using Burkardt.Types;

namespace Burkardt.Function
{
    public static class Sigma
    {
        public static int sigma ( int n )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIGMA returns the value of SIGMA(N), the divisor sum.
            //
            //  Discussion:
            //
            //    SIGMA(N) is the sum of the distinct divisors of N, including 1 and N.
            //
            //    The formula is:
            //
            //      SIGMA(U*V) = SIGMA(U) * SIGMA(V) if U and V are relatively prime.
            //
            //      SIGMA(P**K) = ( P^(K+1) - 1 ) / ( P - 1 ) if P is prime.
            //
            //  First values:
            //
            //     N  SIGMA(N)
            //
            //     1    1
            //     2    3
            //     3    4
            //     4    7
            //     5    6
            //     6   12
            //     7    8
            //     8   15
            //     9   13
            //    10   18
            //    11   12
            //    12   28
            //    13   14
            //    14   24
            //    15   24
            //    16   31
            //    17   18
            //    18   39
            //    19   20
            //    20   42
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    11 February 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the value to be analyzed.
            //
            //    Output, int SIGMA, the value of SIGMA(N).  If N is less than
            //    or equal to 0, SIGMA will be returned as 0.  If there is not
            //    enough room for factoring N, SIGMA is returned as -1.
            //
        {
            int FACTOR_MAX = 20;

            int[] factor = new int[FACTOR_MAX];
            int i;
            int nfactor = 0;
            int nleft = 0;
            int[] power = new int[FACTOR_MAX];
            int value;

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
                Console.WriteLine("SIGMA - Fatal error!");
                Console.WriteLine("  Not enough factorization space.");
                return ( 1 );
            }

            value = 1;
            for ( i = 0; i < nfactor; i++ )
            {
                value = ( value * 
                          ( ( int ) Math.Pow ( ( double ) factor[i], power[i] + 1 ) - 1 ) ) 
                        / ( factor[i] - 1 );
            }

            return value;
        }

    }
}