using Burkardt.BLAS;

namespace Burkardt.Linpack
{
    public static class DPOSL
    {
        public static void dposl(double[] a, int lda, int n, double[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DPOSL solves a linear system factored by DPOCO or DPOFA.
        //
        //  Discussion:
        //
        //    To compute inverse(A) * C where C is a matrix with P columns:
        //
        //      call dpoco ( a, lda, n, rcond, z, info )
        //
        //      if ( rcond is not too small .and. info == 0 ) then
        //        do j = 1, p
        //          call dposl ( a, lda, n, c(1,j) )
        //        end do
        //      end if
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
        //    23 May 2005
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
        //    Input, double A[LDA*N], the output from DPOCO or DPOFA.
        //
        //    Input, int LDA, the leading dimension of the array A.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input/output, double B[N].  On input, the right hand side.
        //    On output, the solution.
        //
        {
            int k;
            double t;
            //
            //  Solve R' * Y = B.
            //
            for (k = 1; k <= n; k++)
            {
                t = BLAS1D.ddot(k - 1, a, 1, b, 1, xIndex: + 0 + (k - 1) * lda);
                b[k - 1] = (b[k - 1] - t) / a[k - 1 + (k - 1) * lda];
            }

            //
            //  Solve R * X = Y.
            //
            for (k = n; 1 <= k; k--)
            {
                b[k - 1] = b[k - 1] / a[k - 1 + (k - 1) * lda];
                t = -b[k - 1];
                BLAS1D.daxpy(k - 1, t, a, 1, ref b, 1, xIndex: + 0 + (k - 1) * lda);
            }
        }
    }
}