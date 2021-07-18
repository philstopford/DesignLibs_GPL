using System;
using Burkardt.BLAS;

namespace Burkardt.Linpack
{
    public static class DPBSL
    {
        public static void dpbsl(double[] abd, int lda, int n, int m, ref double[] b)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DPBSL solves a real SPD band system factored by DPBCO or DPBFA.
            //
            //  Discussion:
            //
            //    The matrix is assumed to be a symmetric positive definite (SPD)
            //    band matrix.
            //
            //    To compute inverse(A) * C  where C is a matrix with P columns:
            //
            //      call dpbco ( abd, lda, n, rcond, z, info )
            //
            //      if ( rcond is too small .or. info /= 0) go to ...
            //
            //      do j = 1, p
            //        call dpbsl ( abd, lda, n, c(1,j) )
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
            //    Input, double ABD[LDA*N], the output from DPBCO or DPBFA.
            //
            //    Input, int LDA, the leading dimension of the array ABD.
            //
            //    Input, int N, the order of the matrix.
            //
            //    Input, int M, the number of diagonals above the main diagonal.
            //
            //    Input/output, double B[N].  On input, the right hand side.
            //    On output, the solution.
            //
        {
            int k;
            int la;
            int lb;
            int lm;
            double t;
            //
            //  Solve R'*Y = B.
            //
            for (k = 1; k <= n; k++)
            {
                lm = Math.Min(k - 1, m);
                la = m + 1 - lm;
                lb = k - lm;
                t = BLAS1D.ddot(lm, abd, 1, b, 1, xIndex: +la - 1 + (k - 1) * lda, yIndex: +lb - 1);
                b[k - 1] = (b[k - 1] - t) / abd[m + (k - 1) * lda];
            }

            //
            //  Solve R*X = Y.
            //
            for (k = n; 1 <= k; k--)
            {
                lm = Math.Min(k - 1, m);
                la = m + 1 - lm;
                lb = k - lm;
                b[k - 1] = b[k - 1] / abd[m + (k - 1) * lda];
                t = -b[k - 1];
                BLAS1D.daxpy(lm, t, abd, 1, ref b, 1, xIndex: +la - 1 + (k - 1) * lda, yIndex: +lb - 1);
            }
        }
    }
}