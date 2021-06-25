namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double[] r8vec_conjugate ( int n, double[] c )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_CONJUGATE reverses a vector and negates even-indexed entries.
            //
            //  Discussion:
            //
            //    There are many times in wavelet computations when such an operation
            //    is invoked.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 April 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the dimension of the vector.
            //
            //    Input, double C[N], the input vector.
            //
            //    Output, double D[N], the "conjugated" vector.
            //
        {
            double[] d;
            int i;

            d = new double[n];

            for ( i = 0; i < n; i++ )
            {
                d[i] = c[n-1-i];
            }
            for ( i = 1; i < n; i = i + 2 )
            {
                d[i] = - d[i];
            }

            return d;
        }
    }
}