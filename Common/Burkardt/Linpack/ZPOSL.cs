using System.Numerics;
using Burkardt.BLAS;

namespace Burkardt.Linpack
{
    public static class ZPOSL
    {
        public static void zposl ( Complex[] a, int lda, int n, ref Complex[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZPOSL solves a complex hermitian positive definite system.
        //
        //  Discussion:
        //
        //    ZPOSL uses the factors computed by ZPOCO or ZPOFA.
        //
        //    A division by zero will occur if the input factor contains
        //    a zero on the diagonal.  Technically this indicates
        //    singularity but it is usually caused by improper subroutine
        //    arguments.  It will not occur if the subroutines are called
        //    correctly and INFO == 0.
        //
        //    To compute inverse(A) * C where C is a matrix with  p  columns
        //
        //      call zpoco(a,lda,n,rcond,z,info)
        //
        //      if (rcond is too small .or. info /= 0) then
        //        error
        //      end if
        //
        //      do j = 1, p
        //        call zposl(a,lda,n,c(1,j))
        //      end do
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 May 2006
        //
        //  Author:
        //
        //    C++ version by John Burkardt
        //
        //  Reference:
        //
        //    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
        //    LINPACK User's Guide,
        //    SIAM, (Society for Industrial and Applied Mathematics),
        //    3600 University City Science Center,
        //    Philadelphia, PA, 19104-2688.
        //
        //  Parameters:
        //
        //    Input, Complex A[LDA*N], the output from ZPOCO or ZPOFA.
        //
        //    Input, int LDA, the leading dimension of A.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input/output, Complex B[N].  On input, the right hand side.
        //    On output, the solution.
        //
        {
            int k;
            Complex t;
            //
            //  Solve hermitian(R) * Y = B.
            //
            for ( k = 1; k <= n; k++ )
            {
                t = BLAS1Z.zdotc ( k-1, a, 1, b, 1, xIndex:+0+(k-1)*lda );
                b[k-1] = ( b[k-1] - t ) / a[k-1+(k-1)*lda];
            }
            //
            //  Solve R * X = Y.
            //
            for ( k = n; 1 <= k; k-- )
            {
                b[k-1] = b[k-1] / a[k-1+(k-1)*lda];
                t = -b[k-1];
                BLAS1Z.zaxpy ( k-1, t, a, 1, ref b, 1 , xIndex:+0+(k-1)*lda);
            }
        }

    }
}