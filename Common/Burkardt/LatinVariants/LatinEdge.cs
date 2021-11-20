using Burkardt.Types;

namespace Burkardt.Latin;

public static partial class LatinVariants
{
    public static double[] latin_edge ( int dim_num, int point_num, ref int seed )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LATIN_EDGE returns edge points in a Latin square.
        //
        //  Discussion:
        //
        //    In each spatial dimension, there will be exactly one
        //    point with the coordinate value 
        //
        //      ( 0, 1, 2, ..., point_num-1 ) / ( point_num - 1 )
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
        //    Input, int POINT_NUM, the number of points, which should be
        //    at least 2!
        //
        //    Input/output, int *SEED, a seed for UNIFORM.
        //
        //    Output, double X[DIM_NUM*POINT_NUM], the points.
        //
    {
        const int base_ = 0;
        int i;
        int k;
        double[] x = new double[dim_num * point_num];
        switch (point_num)
        {
            //
            //  For spatial dimension I, 
            //    pick a random permutation of 1 to POINT_NUM,
            //    force the corresponding I-th components of X to lie in the
            //    interval ( PERM[J]-1, PERM[J] ) / POINT_NUM.
            //
            case 1:
            {
                k = 0;
                for ( i = 0; i < dim_num; i++ )
                {
                    x[k] = 0.5;
                    k += 1;
                }

                break;
            }
            default:
            {
                k = 0;
                for ( i = 0; i < dim_num; i++ )
                {
                    int[] perm = typeMethods.perm_uniform ( point_num, base_, ref seed );

                    int j;
                    for ( j = 0; j < point_num; j++ )
                    {
                        x[k] = ( double ) perm[j] / ( float ) ( point_num - 1 );
                        k += 1;
                    }
                }

                break;
            }
        }

        return x;
    }
}