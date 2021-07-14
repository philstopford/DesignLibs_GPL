using Burkardt;
using Burkardt.Types;

namespace Burkardt.Error
{
    public static class Eigen
    {
        public static double eigen_error ( int n, int k, double[] a, double[] x, double[] lambda )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EIGEN_ERROR determines the error in a (right) eigensystem.
        //
        //  Discussion:
        //
        //    An R8MAT is a matrix of double values.
        //
        //    This routine computes the Frobenius norm of
        //
        //      A * X - X * LAMBDA
        //
        //    where
        //
        //      A is an N by N matrix,
        //      X is an N by K matrix (each of K columns is an eigenvector)
        //      LAMBDA is a K by K diagonal matrix of eigenvalues.
        //
        //    This routine assumes that A, X and LAMBDA are all real!
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
        //    Input, int K, the number of eigenvectors.
        //    K is usually 1 or N.
        //
        //    Input, double A[N*N], the matrix.
        //
        //    Input, double X[N*K]), the K eigenvectors.
        //
        //    Input, double LAMBDA[K], the K eigenvalues.
        //
        //    Output, double EIGEN_ERROR, the Frobenius norm
        //    of the difference matrix A * X - X * LAMBDA, which would be exactly zero
        //    if X and LAMBDA were exact eigenvectors and eigenvalues of A.
        //
        {
            double[] c;
            int i;
            int j;
            double value;

            c = typeMethods.r8mat_mm_new ( n, n, k, a, x );

            for ( j = 0; j < k; j++ )
            {
                for ( i = 0; i < n; i++ )
                {
                    c[i+n*j] = c[i+n*j] - lambda[j] * x[i+n*j];
                }
            }

            value = typeMethods.r8mat_norm_fro ( n, k, c );

            return value;
        }
    }
}