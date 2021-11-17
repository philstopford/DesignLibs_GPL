using System;
using System.Numerics;
using Burkardt.BLAS;

namespace Burkardt.Linpack;

public static class ZGBSL
{
    public static void zgbsl(Complex[] abd, int lda, int n, int ml, int mu,
            int[] ipvt, ref Complex[] b, int job)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZGBSL solves a complex band system factored by ZGBCO or ZGBFA.
        //
        //  Discussion:
        //
        //    ZGBSL can solve A * X = B or hermitan ( A ) * X = B.
        //
        //    A division by zero will occur if the input factor contains a
        //    zero on the diagonal.  Technically this indicates singularity
        //    but it is often caused by improper arguments or improper
        //    setting of LDA.  It will not occur if the subroutines are
        //    called correctly and if ZGBCO has set 0.0 < RCOND
        //    or ZGBFA has set INFO = 0.
        //
        //    To compute inverse ( A ) * C where C is a matrix with P columns:
        //
        //      call zgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
        //
        //      if ( rcond is not too small ) then
        //        do j = 1, p
        //          call zgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
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
        //    Input, Complex ABD[LDA*N], the output from ZGBCO or ZGBFA.
        //
        //    Input, int LDA, the leading dimension of ABD.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int ML, the number of diagonals below the main diagonal.
        //
        //    Input, int MU, the number of diagonals above the main diagonal.
        //
        //    Input, int IPVT[N], the pivot vector from ZGBCO or ZGBFA.
        //
        //    Input/output, Complex B[N].  On input, the right hand side.
        //    On output, the solution.
        //
        //    Input, int JOB.
        //    0, to solve A*x = b,
        //    nonzero, to solve hermitian(A)*x = b, where hermitian(A) is the
        //    conjugate transpose.
        //
    {
        int k;
        int l;
        int la;
        int lb;
        int lm;
        int m;
        Complex t;

        m = mu + ml + 1;
        switch (job)
        {
            //
            //  JOB = 0, solve A * X = B.
            //
            case 0:
            {
                //
                //  First solve L * Y = B.
                //
                if (ml != 0)
                {
                    for (k = 1; k <= n - 1; k++)
                    {
                        lm = Math.Min(ml, n - k);
                        l = ipvt[k - 1];
                        t = b[l - 1];

                        if (l != k)
                        {
                            b[l - 1] = b[k - 1];
                            b[k - 1] = t;
                        }

                        BLAS1Z.zaxpy(lm, t, abd, 1, ref b, 1, xIndex: +m + (k - 1) * lda, yIndex: +k);
                    }
                }

                //
                //  Now solve U * X = Y.
                //
                for (k = n; 1 <= k; k--)
                {
                    b[k - 1] /= abd[m - 1 + (k - 1) * lda];
                    lm = Math.Min(k, m) - 1;
                    la = m - lm;
                    lb = k - lm;
                    t = -b[k - 1];
                    BLAS1Z.zaxpy(lm, t, abd, 1, ref b, 1, xIndex: +la - 1 + (k - 1) * lda, yIndex: +lb - 1);
                }

                break;
            }
            //
            default:
            {
                //
                //  First solve hermitian ( U ) * Y = B.
                //
                for (k = 1; k <= n; k++)
                {
                    lm = Math.Min(k, m) - 1;
                    la = m - lm;
                    lb = k - lm;
                    t = BLAS1Z.zdotc(lm, abd, 1, b, 1, xIndex: +la - 1 + (k - 1) * lda, yIndex: +lb - 1);
                    b[k - 1] = (b[k - 1] - t) / Complex.Conjugate(abd[m - 1 + (k - 1) * lda]);
                }

                //
                //  Now solve hermitian ( L ) * X = Y.
                //
                if (ml != 0)
                {
                    for (k = n - 1; 1 <= k; k--)
                    {
                        lm = Math.Min(ml, n - k);
                        b[k - 1] += BLAS1Z.zdotc(lm, abd, 1, b, 1, xIndex: +m + (k - 1) * lda, yIndex: +k);
                        l = ipvt[k - 1];

                        if (l != k)
                        {
                            t = b[l - 1];
                            b[l - 1] = b[k - 1];
                            b[k - 1] = t;
                        }
                    }
                }

                break;
            }
        }
    }

}