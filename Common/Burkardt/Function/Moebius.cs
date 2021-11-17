using System;
using Burkardt.Types;

namespace Burkardt.Function;

public static class Moebius
{
    public static int moebius ( int n )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MOEBIUS returns the value of MU(N), the Moebius function of N.
        //
        //  Discussion:
        //
        //    MU(N) is defined as follows:
        //
        //      MU(N) = 1 if N = 1;
        //              0 if N is divisible by the square of a prime;
        //              (-1)^K, if N is the product of K distinct primes.
        //
        //    As special cases, MU(N) is -1 if N is a prime, and MU(N) is 0
        //    if N is a square, cube, etc.
        //
        //    The Moebius function is related to Euler's totient function:
        //
        //      PHI(N) = Sum ( D divides N ) MU(D) * ( N / D ).
        //
        //  First values:
        //
        //     N  MU(N)
        //
        //     1    1
        //     2   -1
        //     3   -1
        //     4    0
        //     5   -1
        //     6    1
        //     7   -1
        //     8    0
        //     9    0
        //    10    1
        //    11   -1
        //    12    0
        //    13   -1
        //    14    1
        //    15    1
        //    16    0
        //    17   -1
        //    18    0
        //    19   -1
        //    20    0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 March 1999
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the value to be analyzed.
        //
        //    Output, int MOEBIUS, the value of MU(N).
        //    If N is less than or equal to 0, or there was not enough internal 
        //    space for factoring, MOEBIUS is returned as -1.
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
                return -1;
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
            Console.WriteLine("MOEBIUS - Fatal error!");
            Console.WriteLine( "  Not enough factorization space.");
            return 1;
        }

        value = 1;

        for ( i = 0; i < nfactor; i++ )
        {
            value = -value;

            switch (power[i])
            {
                case > 1:
                    return 0;
            }
        }

        return value;
    }

}