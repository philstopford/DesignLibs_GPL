using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.RandomNS
{
    public static partial class BRandom
    {
        public static double[] bad_in_simplex01 ( int dim_num, int point_num, ref int seed )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BAD_IN_SIMPLEX01 is a "bad" (nonuniform) sampling of the unit simplex.
            //
            //  Discussion:
            //
            //    The interior of the unit DIM_NUM-dimensional simplex is the set of
            //    points X(1:DIM_NUM) such that each X(I) is nonnegative, and
            //    sum(X(1:DIM_NUM)) <= 1.
            //
            //    Any point in the unit simplex CAN be chosen by this algorithm.
            //
            //    However, the points that are chosen tend to be clustered near
            //    the centroid.
            //
            //    This routine is supplied as an example of "bad" sampling.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    15 August 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the dimension of the space.
            //
            //    Input, int POINT_NUM, the number of points.
            //
            //    Input/output, int &SEED, a seed for the random
            //    number generator.
            //
            //    Output, double BAD_IN_SIMPLEX01[DIM_NUM*POINT_NUM], the points.
            //
        {
            double[] e;
            double e_sum;
            int i;
            int j;
            double[] x;

            x = new double[dim_num*point_num];

            for ( j = 0; j < point_num; j++ )
            {
                e = UniformRNG.r8vec_uniform_01_new ( dim_num + 1, ref seed );

                e_sum = typeMethods.r8vec_sum ( dim_num + 1, e );

                for ( i = 0; i < dim_num; i++ )
                {
                    x[i+j*dim_num] = e[i] / e_sum;
                }
            }

            return x;
        } 
    }
}