namespace Burkardt.Uniform
{
    public static partial class UniformRNG
    {
        public static void r83vec_uniform ( int n, double[] alo, double[] ahi, ref int seed, ref double[] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83VEC_UNIFORM returns a random R83VEC in a given range.
        //
        //  Discussion:
        //
        //    A is a two dimensional array of order N by 3, stored as a vector
        //    of rows: A(0,0), A(0,1), A(0,2), // A(1,0), A(1,1), A(1,2) // ...
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the vector.
        //
        //    Input, double ALO[3], AHI[3], the minimum and maximum values.
        //
        //    Input/output, int *SEED, a seed for the random number generator.
        //
        //    Output, double A[N*3], the vector of randomly chosen values.
        //
        {
            for ( int i = 0; i < n; i++ )
            {
                for ( int j = 0; j < 3; j++ )
                {
                    a[3*i+j] = r8_uniform ( alo[j], ahi[j], ref seed );
                }
            }
        }
    }
}