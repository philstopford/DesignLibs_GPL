using Burkardt.BLAS;

namespace Burkardt.Linpack;

public static class DPPSL
{
    public static void dppsl ( double[] ap, int n, ref double[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DPPSL solves a real symmetric positive definite system factored by DPPCO or DPPFA.
        //
        //  Discussion:
        //
        //    To compute inverse(A) * C where C is a matrix with P columns
        //
        //      call dppco ( ap, n, rcond, z, info )
        //
        //      if ( rcond is too small .or. info /= 0 ) then
        //        exit
        //      end if
        //
        //      do j = 1, p
        //        call dppsl ( ap, n, c(1,j) )
        //      end do
        //
        //    A division by zero will occur if the input factor contains
        //    a zero on the diagonal.  Technically this indicates
        //    singularity but it is usually caused by improper subroutine
        //    arguments.  It will not occur if the subroutines are called
        //    correctly and INFO == 0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 May 2005
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Jack Dongarra, Cleve Moler, Jim Bunch, 
        //    Pete Stewart.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
        //    LINPACK User's Guide,
        //    SIAM, (Society for Industrial and Applied Mathematics),
        //    3600 University City Science Center,
        //    Philadelphia, PA, 19104-2688.
        //    ISBN 0-89871-172-X
        //
        //  Parameters:
        //
        //    Input, double AP[N*(N+1)/2], the output from DPPCO or DPPFA.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input/output, double B[N].  On input, the right hand side.
        //    On output, the solution.
        //
    {
        int k;
        double t;

        int kk = 0;

        for ( k = 1; k <= n; k++ )
        {
            t = BLAS1D.ddot ( k-1, ap, 1, b, 1, xIndex:+kk );
            kk += k;
            b[k-1] = ( b[k-1] - t ) / ap[kk-1];
        }

        for ( k = n; 1 <= k; k-- )
        {
            b[k-1] /= ap[kk-1];
            kk -= k;
            t = -b[k-1];
            BLAS1D.daxpy ( k-1, t, ap, 1, ref b, 1, xIndex:+kk );
        }
    }
}