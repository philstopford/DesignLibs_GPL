using Burkardt.Uniform;

namespace Burkardt.RandomNS
{
    public static class Permutation
    {
        public static void random_permutation ( int n, ref double[] x, ref int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RANDOM_PERMUTATION applies a random permutation to an array.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 February 2016
        //
        //  Author:
        //
        //    Original C version by Warren Smith.
        //    This C++ version by John Burkardt.
        //
        //  Parameters:
        //
        //    Input, unsigned int N, indicates the size of X.
        //
        //    Input/output, double X[N+2].  On output, entries X[1] through
        //    X[N] have been randomly permuted.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        {
            int i;
            int j;
            double t;

            for ( i = 1; i < n; i++ )
            {
                j = UniformRNG.i4_uniform_ab ( i, n, ref seed );

                t = x[i];
                x[i] = x[j];
                x[j] = t;      
            }
        }
    }
}