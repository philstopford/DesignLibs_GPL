using System;
using Burkardt.Uniform;

namespace Burkardt.RandomNS
{
    public static class Ksub
    {
        public static void ksub_random2 ( int n, int k, ref int seed, ref int[] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    KSUB_RANDOM2 selects a random subset of size K from a set of size N.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 May 2003
        //
        //  Author:
        //
        //    FORTRAN77 original version by Albert Nijenhuis, Herbert Wilf.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    A Nijenhuis and H Wilf,
        //    Combinatorial Algorithms,
        //    Academic Press, 1978, second edition,
        //    ISBN 0-12-519260-6.
        //
        //  Parameters:
        //
        //    Input, int N, the size of the set from which subsets are drawn.
        //
        //    Input, int K, number of elements in desired subsets.  K must
        //    be between 0 and N.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, int A[K].  A(I) is the I-th element of the
        //    output set.  The elements of A are in order.
        //
        {
            int available;
            int candidate;
            int have;
            int need;
            double r;

            if ( k < 0 || n < k )
            {
                Console.WriteLine("");
                Console.WriteLine("KSUB_RANDOM2 - Fatal error!");
                Console.WriteLine("  N = " + n + "");
                Console.WriteLine("  K = " + k + "");
                Console.WriteLine("  but 0 <= K <= N is required!");
                return;
            }

            if ( k == 0 )
            {
                return;
            }

            need = k;
            have = 0;
            available = n;
            candidate = 0;

            for ( ; ; )
            {
                candidate = candidate + 1;

                r = UniformRNG.r8_uniform_01 ( ref seed );

                if ( r * ( double ) available <= ( double ) need )
                {
                    need = need - 1;
                    a[have] = candidate;
                    have = have + 1;

                    if ( need <= 0 )
                    {
                        break;
                    }

                }

                available = available - 1;

            }

            return;
        }
    }
}