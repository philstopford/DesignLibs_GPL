using Burkardt.Types;

namespace Burkardt.Error;

public static class LU
{
    public static double lu_error ( int n, double[] a, double[] l, double[] u )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LU_ERROR determines the error in an LU factorization.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 November 2013
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
        //    Input, double L[N*N], U[N*N], the LU factors.
        //
        //    Output, double LU_ERROR, the Frobenius norm
        //    of the difference matrix A - L * U.
        //
    {
        double[] lu = typeMethods.r8mat_mm_new ( n, n, n, l, u );

        double[] d = typeMethods.r8mat_sub_new ( n, n, a, lu );
 
        double value = typeMethods.r8mat_norm_fro ( n, n, d );

        return value;
    }
}