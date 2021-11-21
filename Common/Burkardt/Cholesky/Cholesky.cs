using Burkardt.Types;

namespace Burkardt.CholeskyNS;

public static class Cholesky
{
    public static double cholesky_upper_error ( int n, double[] a, double[] c )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHOLESKY_UPPER_ERROR determines the error in an upper Cholesky factor.
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
        //    Input, double C[N*N], the upper triangular Cholesky factor.
        //
        //    Output, double CHOLESKY_UPPER_ERROR, the Frobenius norm
        //    of the difference matrix A - C' * C.
        //
    {
        double[] ctc = typeMethods.r8mat_mtm_new ( n, n, n, c, c );

        double[] d = typeMethods.r8mat_sub_new ( n, n, a, ctc );
 
        double value = typeMethods.r8mat_norm_fro ( n, n, d );

        return value;
    }
}