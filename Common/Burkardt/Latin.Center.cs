using entropyRNG;
using Burkardt.Uniform;

namespace Burkardt.Latin.Center
{
    public static class LatinCenter
    {
        static int[] perm_uniform ( int n, int base_, ref int seed )
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

        
        public static double[] latin_center ( int dim_num, int point_num, ref int seed )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LATIN_CENTER returns center points in a Latin square.
        //
        //  Discussion:
        //
        //    In each spatial dimension, there will be exactly one
        //    point with the coordinate value 
        //
        //      ( 1, 3, 5, ..., 2*point_num-1 ) / ( 2 * point_num )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 March 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int POINT_NUM, the number of points.
        //
        //    Input/output, int *SEED, a seed for UNIFORM.
        //
        //    Output, double LATIN_CENTER[DIM_NUM*POINT_NUM], the points.
        //
        {
            int base_ = 0;
            double[] x = new double[dim_num * point_num];
            //
            //  For spatial dimension I, 
            //    pick a random permutation of 1 to POINT_NUM,
            //    force the corresponding I-th components of X to lie in the
            //    interval ( PERM[J]-1, PERM[J] ) / POINT_NUM.
            //
            int k = 0;
            for (int i = 0; i < dim_num; i++ )
            {
                int[] perm = perm_uniform ( point_num, base_, ref seed );

                for (int j = 0; j < point_num; j++ )
                {
                    double r = 0.5;
                    x[k] = ( perm[j] + r ) / point_num;
                    k = k + 1;
                }
            }

            return x;
        }

        
    }
}