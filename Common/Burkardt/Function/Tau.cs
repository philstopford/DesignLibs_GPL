using System;
using Burkardt.Types;

namespace Burkardt.Function;

public static class Tau
{
    public static int tau ( int n )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TAU returns the value of TAU(N), the number of distinct divisors of N.
        //
        //  Discussion:
        //
        //    TAU(N) is the number of divisors of N, including 1 and N.
        //
        //    The formula is:
        //
        //      If the prime factorization of N is
        //
        //        N = P1^E1 * P2^E2 * ... * PM^EM,
        //
        //      then
        //
        //        TAU(N) = ( E1 + 1 ) * ( E2 + 1 ) * ... * ( EM + 1 ).
        //
        //  First values:
        //
        //     N   TAU(N)
        //
        //     1    1
        //     2    2
        //     3    2
        //     4    3
        //     5    2
        //     6    4
        //     7    2
        //     8    4
        //     9    3
        //    10    4
        //    11    2
        //    12    6
        //    13    2
        //    14    4
        //    15    4
        //    16    5
        //    17    2
        //    18    6
        //    19    2
        //    20    6
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 December 1998
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
        //    Output, int TAU, the value of TAU(N).  But if N is 0 or
        //    less, TAU is returned as 0, a nonsense value.  If there is
        //    not enough room for factoring, TAU is returned as -1.
        //
    {
        int FACTOR_MAX = 20;

        int[] factor = new int[FACTOR_MAX];
        int i;
        int nfactor = 0;
        int nleft = 0;
        int[] power = new int[FACTOR_MAX];
        int value;

        switch (n)
        {
            case <= 0:
                return 0;
            case 1:
                return 1;
        }

        //
        //  Factor N.
        //
        typeMethods.i4_factor ( n, FACTOR_MAX, ref nfactor, ref factor, ref power, ref nleft );

        if ( nleft != 1 )
        {
            Console.WriteLine("");
            Console.WriteLine("TAU - Fatal error!");
            Console.WriteLine("  Not enough factorization space.");
            return 1;
        }

        value = 1;
        for ( i = 0; i < nfactor; i++ )
        {
            value *= ( power[i] + 1 );
        }

        return value;
    }

}