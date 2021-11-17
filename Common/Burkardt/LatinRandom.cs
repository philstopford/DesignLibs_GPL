using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.Latin;

public static class Random
{
    public static double[] latin_random_new ( int dim_num, int point_num, ref int seed )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LATIN_RANDOM_NEW returns points in a Latin Random square.
        //
        //  Discussion:
        //
        //    In each spatial dimension, there will be exactly one
        //    point whose coordinate value lies between consecutive
        //    values in the list:
        //
        //      ( 0, 1, 2, ..., point_num ) / point_num
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 April 2003
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
        //    Input/output, int &SEED, a seed for UNIFORM.
        //
        //    Output, double LATIN_RANDOM_NEW[DIM_NUM,POINT_NUM], the points.
        //
    {
        double[] x = UniformRNG.r8mat_uniform_01 ( dim_num, point_num, ref seed );
        //
        //  For spatial dimension I, 
        //    pick a random permutation of 1 to POINT_NUM,
        //    force the corresponding I-th components of X to lie in the
        //    interval ( PERM[J]-1, PERM[J] ) / POINT_NUM.
        //
        for (int i = 0; i < dim_num; i++ )
        {
            int[] perm = typeMethods.perm_uniform ( point_num, 0, ref seed );

            for (int j = 0; j < point_num; j++ )
            {
                x[i+j*dim_num] = ( perm[j] + x[i+j*dim_num] ) 
                                 / point_num;
            }
        }
        return x;
    }

}