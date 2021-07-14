using Burkardt;
using Burkardt.Types;

namespace Burkardt.Error
{
    public static class Inverse
    {
        public static double inverse_error ( int n, double[] a, double[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INVERSE_ERROR determines the error in an inverse matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 October 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[N*N], the matrix.
        //
        //    Input, double B[N*N], the inverse.
        //
        //    Output, double ERROR_FROBENIUS, the Frobenius norm
        //    of (A*B-I) + (B*A-I).
        //
        {
            double[] c;
            int j;
            double value;

            c = typeMethods.r8mat_mm_new ( n, n, n, a, b );

            for ( j = 0; j < n; j++ )
            {
                c[j+j*n] = c[j+j*n] - 1.0;
            }

            value = typeMethods.r8mat_norm_fro ( n, n, c );

            c = typeMethods.r8mat_mm_new ( n, n, n, b, a );

            for ( j = 0; j < n; j++ )
            {
                c[j+j*n] = c[j+j*n] - 1.0;
            }

            value = value + typeMethods.r8mat_norm_fro ( n, n, c );

            return value;
        }
    }
}