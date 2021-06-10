namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static void r8vec_binary_next ( int n, double[] bvec )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_BINARY_NEXT generates the next binary vector.
            //
            //  Discussion:
            //
            //    The vectors have the order
            //
            //      (0,0,...,0),
            //      (0,0,...,1),
            //      ...
            //      (1,1,...,1)
            //
            //    and the "next" vector after (1,1,...,1) is (0,0,...,0).  That is,
            //    we allow wrap around.
            //
            //  Example:
            //
            //    N = 3
            //
            //    Input      Output
            //    -----      ------
            //    0 0 0  =>  0 0 1
            //    0 0 1  =>  0 1 0
            //    0 1 0  =>  0 1 1
            //    0 1 1  =>  1 0 0
            //    1 0 0  =>  1 0 1
            //    1 0 1  =>  1 1 0
            //    1 1 0  =>  1 1 1
            //    1 1 1  =>  0 0 0
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    31 March 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the dimension of the vectors.
            //
            //    Input, double BVEC[N], the vector whose successor is desired.
            //
            //    Output, double BVEC[N], the successor to the input vector.
            //
        {
            int i;

            for ( i = n - 1; 0 <= i; i-- )
            {
                if ( bvec[i] == 0.0 )
                {
                    bvec[i] = 1.0;
                    return;
                }
                bvec[i] = 0.0;
            }

            return;
        }

    }
}