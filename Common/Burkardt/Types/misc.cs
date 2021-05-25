using System;
using Burkardt.Uniform;
using entropyRNG;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static int get_seed()
        {
            return RNG.nextint(1, Int32.MaxValue);
        }

        public static int[] perm_uniform ( int n, int base_, ref int seed )
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM_UNIFORM selects a random permutation of N objects.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    31 October 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Albert Nijenhuis, Herbert Wilf,
            //    Combinatorial Algorithms,
            //    Academic Press, 1978, second edition,
            //    ISBN 0-12-519260-6.
            //
            //  Parameters:
            //
            //    Input, int N, the number of objects to be permuted.
            //
            //    Input, int BASE, is 0 for a 0-based permutation and 1 for 
            //    a 1-based permutation.
            //
            //    Input/output, int *SEED, a seed for the random number generator.
            //
            //    Output, int PERM_UNIFORM[N], a permutation of (BASE, BASE+1, ..., BASE+N-1).
            //
        {
            int[] p = new int[n];
 
            for (int i = 0; i < n; i++ )
            {
                p[i] = i + base_;
            }

            for (int i = 0; i < n; i++ )
            {
                int j = UniformRNG.i4_uniform( i, n - 1, ref seed );
                int k    = p[i];
                p[i] = p[j];
                p[j] = k;
            }
 
            return p;
        }

    }
}