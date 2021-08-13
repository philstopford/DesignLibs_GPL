using Burkardt.Uniform;

namespace Burkardt.Sequence
{
    public static class Vibonacci
    {
        public static void vibonacci ( int n, ref int seed, ref int[] v )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VIBONACCI computes the first N Vibonacci numbers.
        //
        //  Discussion:
        //
        //    The "Vibonacci numbers" are a generalization of the Fibonacci numbers:
        //      V(N+1) = +/- V(N) +/- V(N-1)
        //    where the signs are chosen randomly.
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
        //    Brian Hayes,
        //    The Vibonacci Numbers,
        //    American Scientist,
        //    July-August 1999, Volume 87, Number 4.
        //
        //    Divakar Viswanath,
        //    Random Fibonacci sequences and the number 1.13198824,
        //    Mathematics of Computation, 1998.
        //
        //  Parameters:
        //
        //    Input, int N, the highest number to compute.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, int V(N), the first N Vibonacci numbers.  By convention,
        //    V(1) and V(2) are taken to be 1.
        //
        {
            int i;
            int j;
            int s1;
            int s2;

            if ( n <= 0 )
            {
                return;
            }

            v[0] = 1;

            if ( n <= 1 )
            {
                return;
            }

            v[1] = 1;

            for ( i = 2; i < n; i++ )
            {
                j = UniformRNG.i4_uniform_ab ( 0, 1, ref seed );

                if ( j == 0 )
                {
                    s1 = -1;
                }
                else
                {
                    s1 = +1;
                }

                j = UniformRNG.i4_uniform_ab ( 0, 1, ref seed );

                if ( j == 0 )
                {
                    s2 = -1;
                }
                else
                {
                    s2 = +1;
                }

                v[i] = s1 * v[i-1] + s2 * v[i-2];

            }
        }

    }
}