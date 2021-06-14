using System;
using Burkardt.Types;

namespace Burkardt.RandomNS
{
    public static partial class BRandom
    {
        public static void direction_uniform_nd ( int dim_num, ref int seed, double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIRECTION_UNIFORM_ND generates a random direction vector in ND.
        //
        //  Discussion:
        //
        //    This is actually simply a random point on the unit sphere in ND.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 August 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the dimension of the space.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double W[DIM_NUM], a random direction vector, with unit norm.
        //
        {
            int i;
            double norm;
            //
            //  Sample the standard normal distribution.
            //
            typeMethods.r8vec_normal_01 ( dim_num, ref seed, w );
            //
            //  Compute the length of the vector.
            //
            norm = 0.0;
            for ( i = 0; i < dim_num; i++ )
            {
                norm = norm + w[i] * w[i];
            }
            norm = Math.Sqrt ( norm );
            //
            //  Normalize the vector.
            //
            for ( i = 0; i < dim_num; i++ )
            {
                w[i] = w[i] / norm;
            }
        }
    }
}