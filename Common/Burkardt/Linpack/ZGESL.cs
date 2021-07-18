using System.Numerics;
using Burkardt.BLAS;

namespace Burkardt.Linpack
{
    public static class ZGESL
    {
        public static void zgesl(Complex[] a, int lda, int n, int[] ipvt,
                ref Complex[] b, int job)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ZGESL solves a complex system factored by ZGECO or ZGEFA.
            //
            //  Discussion:
            //
            //    A division by zero will occur if the input factor contains a
            //    zero on the diagonal.  Technically this indicates singularity
            //    but it is often caused by improper arguments or improper
            //    setting of LDA.  It will not occur if the subroutines are
            //    called correctly and if ZGECO has set 0.0 < RCOND
            //    or ZGEFA has set INFO == 0.
            //
            //    To compute inverse(A) * C where C is a matrix with P columns:
            //
            //      call zgeco(a,lda,n,ipvt,rcond,z)
            //
            //      if (rcond is not too small) then
            //        do j = 1, p
            //          call zgesl ( a, lda, n, ipvt, c(1,j), 0 )
            //        end do
            //      end if
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
            //    Input, Complex A[LDA*N], the factored matrix information,
            //    as output from ZGECO or ZGEFA.
            //
            //    Input, int LDA, the leading dimension of A.
            //
            //    Input, int N, the order of the matrix.
            //
            //    Input, int IPVT[N], the pivot vector from ZGECO or ZGEFA.
            //
            //    Input/output, Complex B[N].  On input, the right hand side.
            //    On output, the solution.
            //
            //    Input, int JOB.
            //    0, to solve A*X = B,
            //    nonzero, to solve hermitian(A)*X = B where hermitian(A) is the
            //    conjugate transpose.
            //
        {
            int k;
            int l;
            Complex t;
            //
            //  JOB = 0, solve A * X = B.
            //
            //  First solve L * Y = B.
            //
            if (job == 0)
            {
                for (k = 1; k <= n - 1; k++)
                {
                    l = ipvt[k - 1];
                    t = b[l - 1];
                    if (l != k)
                    {
                        b[l - 1] = b[k - 1];
                        b[k - 1] = t;
                    }

                    BLAS1Z.zaxpy(n - k, t, a, 1, ref b, 1, xIndex: +k + (k - 1) * lda, yIndex: +k);
                }

                //
                //  Now solve U * X = Y.
                //
                for (k = n; 1 <= k; k--)
                {
                    b[k - 1] = b[k - 1] / a[k - 1 + (k - 1) * lda];
                    t = -b[k - 1];
                    BLAS1Z.zaxpy(k - 1, t, a, 1, ref b, 1, xIndex: +0 + (k - 1) * lda);
                }
            }
            //
            //  JOB nonzero, solve hermitian(A) * X = B.
            //
            //  First solve hermitian(U) * Y = B.
            //
            else
            {
                for (k = 1; k <= n; k++)
                {
                    t = BLAS1Z.zdotc(k - 1, a, 1, b, 1, xIndex: +0 + (k - 1) * lda);
                    b[k - 1] = (b[k - 1] - t) / Complex.Conjugate(a[k - 1 + (k - 1) * lda]);
                }

                //
                //  Now solve hermitian(L) * X = Y.
                //
                for (k = n - 1; 1 <= k; k--)
                {
                    b[k - 1] = b[k - 1] + BLAS1Z.zdotc(n - k, a, 1, b, 1, xIndex: +k + (k - 1) * lda, yIndex: +k);
                    l = ipvt[k - 1];
                    if (l != k)
                    {
                        t = b[l - 1];
                        b[l - 1] = b[k - 1];
                        b[k - 1] = t;
                    }
                }
            }
        }

    }
}