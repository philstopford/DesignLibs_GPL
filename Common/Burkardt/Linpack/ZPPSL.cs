using System.Numerics;
using Burkardt.BLAS;

namespace Burkardt.Linpack;

public static class ZPPSL
{
    public static void zppsl ( Complex[] ap, int n, ref Complex[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZPPSL solves a complex hermitian positive definite linear system.
        //
        //  Discussion:
        //
        //    The matrix is assumed to have been factored by ZPPCO or ZPPFA.
        //
        //    A division by zero will occur if the input factor contains
        //    a zero on the diagonal.  Technically this indicates
        //    singularity but it is usually caused by improper subroutine
        //    arguments.  It will not occur if the subroutines are called
        //    correctly and INFO == 0.
        //
        //    To compute inverse(A) * C where C is a matrix with P columns:
        //
        //      call zppco(ap,n,rcond,z,info)
        //
        //      if (rcond is too small .or. info /= 0) then
        //        error
        //      end if
        //
        //      do j = 1, p
        //        call zppsl(ap,n,c(1,j))
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
        //    Input, Complex AP[N*(N+1)/2], the output from ZPPCO or ZPPFA.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input/output, Complex B[N].  On input, the right hand side.
        //    On output, the solution.
        //
    {
        int k;
        Complex t;

        int kk = 0;
        for ( k = 1; k <= n; k++ )
        {
            t = BLAS1Z.zdotc ( k-1, ap, 1, b, 1, xIndex:+kk );
            kk += k;
            b[k-1] = ( b[k-1] - t ) / ap[kk-1];
        }

        for ( k = n; 1 <= k; k-- )
        {
            b[k-1] /= ap[kk-1];
            kk -= k;
            t = -b[k-1];
            BLAS1Z.zaxpy ( k-1, t, ap, 1, ref b, 1, xIndex:+kk );
        }
    }

}