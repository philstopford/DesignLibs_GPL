using System;
using Burkardt.BLAS;
using Burkardt.Types;

namespace Burkardt.Linpack
{
    public static class DSICO
    {
        public static double dsico(ref double[] a, int lda, int n, ref int[] kpvt, ref double[] z)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DSICO factors a real symmetric matrix and estimates its condition.
            //
            //  Discussion:
            //
            //    If RCOND is not needed, DSIFA is slightly faster.
            //
            //    To solve A * X = B, follow DSICO by DSISL.
            //
            //    To compute inverse(A)*C, follow DSICO by DSISL.
            //
            //    To compute inverse(A), follow DSICO by DSIDI.
            //
            //    To compute determinant(A), follow DSICO by DSIDI.
            //
            //    To compute inertia(A), follow DSICO by DSIDI.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 June 2005
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
            //    Input/output, double A[LDA*N].  On input, the symmetric
            //    matrix to be factored.  Only the diagonal and upper triangle are used.
            //    On output, a block diagonal matrix and the multipliers which
            //    were used to obtain it.  The factorization can be written A = U*D*U'
            //    where U is a product of permutation and unit upper triangular
            //    matrices, U' is the transpose of U, and D is block diagonal
            //    with 1 by 1 and 2 by 2 blocks.
            //
            //    Input, int LDA, the leading dimension of the array A.
            //
            //    Input, int N, the order of the matrix.
            //
            //    Output, int KPVT[N], pivot indices.
            //
            //    Output, double Z[N], a work vector whose contents are usually
            //    unimportant.  If A is close to a singular matrix, then Z is an
            //    approximate null vector in the sense that
            //      norm(A*Z) = RCOND * norm(A) * norm(Z).
            //
            //    Output, double DSICO, an estimate of the reciprocal condition number RCOND
            //    of A.  For the system A*X = B, relative perturbations in A and B of size
            //    EPSILON may cause relative perturbations in X of size EPSILON/RCOND.
            //    If RCOND is so small that the logical expression
            //      1.0 + RCOND == 1.0D+00
            //    is true, then A may be singular to working precision.  In particular,
            //    RCOND is zero if exact singularity is detected or the estimate underflows.
            //
        {
            double ak;
            double akm1;
            double anorm;
            double bk;
            double bkm1;
            double denom;
            double ek;
            int i;
            int j;
            int k;
            int kp;
            int kps;
            int ks;
            double rcond;
            double s;
            double t;
            double ynorm;
            //
            //  Find the norm of A, using only entries in the upper half of the matrix.
            //
            for (j = 1; j <= n; j++)
            {
                z[j - 1] = BLAS1D.dasum(j, a, 1, index: +0 + (j - 1) * lda);
                for (i = 1; i <= j - 1; i++)
                {
                    z[i - 1] = z[i - 1] + Math.Abs(a[i - 1 + (j - 1) * lda]);
                }
            }

            anorm = 0.0;
            for (i = 1; i <= n; i++)
            {
                anorm = Math.Max(anorm, z[i - 1]);
            }

            //
            //  Factor.
            //
            DSIFA.dsifa(ref a, lda, n, ref kpvt);
            //
            //  RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))).
            //
            //  Estimate = norm(Z)/norm(Y) where A*Z = Y and A*Y = E.
            //
            //  The components of E are chosen to cause maximum local
            //  growth in the elements of W where U*D*W = E.
            //
            //  The vectors are frequently rescaled to avoid overflow.
            //
            //  Solve U * D * W = E.
            //
            ek = 1.0;
            for (i = 1; i <= n; i++)
            {
                z[i - 1] = 0.0;
            }

            k = n;

            while (k != 2)
            {
                if (kpvt[k - 1] < 0)
                {
                    ks = 2;
                }
                else
                {
                    ks = 1;
                }

                kp = Math.Abs(kpvt[k - 1]);
                kps = k + 1 - ks;

                if (kp != kps)
                {
                    t = z[kps - 1];
                    z[kps - 1] = z[kp - 1];
                    z[kp - 1] = t;
                }

                if (z[k - 1] != 0.0)
                {
                    ek = ek * typeMethods.r8_sign(z[k - 1]);
                }

                z[k - 1] = z[k - 1] + ek;
                BLAS1D.daxpy(k - ks, z[k - 2], a, 1, ref z, 1, xIndex: +0 + (k - 1) * lda);

                if (ks != 1)
                {
                    if (z[k - 2] != 0.0)
                    {
                        ek = ek * typeMethods.r8_sign(z[k - 2]);
                    }

                    z[k - 2] = z[k - 2] + ek;
                    BLAS1D.daxpy(k - ks, z[k - 2], a, 1, ref z, 1, xIndex: +0 + (k - 2) * lda);
                }

                if (ks != 2)
                {
                    if (Math.Abs(a[k - 1 + (k - 1) * lda]) < Math.Abs(z[k - 1]))
                    {
                        s = Math.Abs(a[k - 1 + (k - 1) * lda]) / Math.Abs(z[k - 1]);
                        for (i = 1; i <= n; i++)
                        {
                            z[i - 1] = s * z[i - 1];
                        }

                        ek = s * ek;
                    }

                    if (a[k - 1 + (k - 1) * lda] != 0.0)
                    {
                        z[k - 1] = z[k - 1] / a[k - 1 + (k - 1) * lda];
                    }
                    else
                    {
                        z[k - 1] = 1.0;
                    }
                }
                else
                {
                    ak = a[k - 1 + (k - 1) * lda] / a[k - 2 + (k - 1) * lda];
                    akm1 = a[k - 2 + (k - 2) * lda] / a[k - 2 + (k - 1) * lda];
                    bk = z[k - 1] / a[k - 2 + (k - 1) * lda];
                    bkm1 = z[k - 2] / a[k - 2 + (k - 1) * lda];
                    denom = ak * akm1 - 1.0;
                    z[k - 1] = (akm1 * bk - bkm1) / denom;
                    z[k - 2] = (ak * bkm1 - bk) / denom;
                }

                k = k - ks;
            }

            s = BLAS1D.dasum(n, z, 1);
            for (i = 1; i <= n; i++)
            {
                z[i - 1] = z[i - 1] / s;
            }

            //
            //  Solve U' * Y = W.
            //
            k = 1;

            while (k <= n)
            {
                if (kpvt[k - 1] < 0)
                {
                    ks = 2;
                }
                else
                {
                    ks = 1;
                }

                if (k != 1)
                {
                    z[k - 1] = z[k - 1] + BLAS1D.ddot(k - 1, a, 1, z, 1, xIndex: +0 + (k - 1) * lda);

                    if (ks == 2)
                    {
                        z[k] = z[k] + BLAS1D.ddot(k - 1, a, 1, z, 1, xIndex: +0 + k * lda);
                    }

                    kp = Math.Abs(kpvt[k - 1]);

                    if (kp != k)
                    {
                        t = z[k - 1];
                        z[k - 1] = z[kp - 1];
                        z[kp - 1] = t;
                    }
                }

                k = k + ks;
            }

            s = BLAS1D.dasum(n, z, 1);
            for (i = 1; i <= n; i++)
            {
                z[i - 1] = z[i - 1] / s;
            }

            ynorm = 1.0;
            //
            //  Solve U * D * V = Y.
            //
            k = n;

            while (k != 0)
            {
                if (kpvt[k - 1] < 0)
                {
                    ks = 2;
                }
                else
                {
                    ks = 1;
                }

                if (k != ks)
                {
                    kp = Math.Abs(kpvt[k - 1]);
                    kps = k + 1 - ks;

                    if (kp != kps)
                    {
                        t = z[kps - 1];
                        z[kps - 1] = z[kp - 1];
                        z[kp - 1] = t;
                    }

                    BLAS1D.daxpy(k - ks, z[k - 1], a, 1, ref z, 1, xIndex: +0 + (k - 1) * lda);

                    if (ks == 2)
                    {
                        BLAS1D.daxpy(k - ks, z[k - 2], a, 1, ref z, 1, xIndex: +0 + (k - 2) * lda);
                    }
                }

                if (ks != 2)
                {
                    if (Math.Abs(a[k - 1 + (k - 1) * lda]) < Math.Abs(z[k - 1]))
                    {
                        s = Math.Abs(a[k - 1 + (k - 1) * lda]) / Math.Abs(z[k - 1]);
                        for (i = 1; i <= n; i++)
                        {
                            z[i - 1] = s * z[i - 1];
                        }

                        ynorm = s * ynorm;
                    }

                    if (a[k - 1 + (k - 1) * lda] != 0.0)
                    {
                        z[k - 1] = z[k - 1] / a[k - 1 + (k - 1) * lda];
                    }
                    else
                    {
                        z[k - 1] = 1.0;
                    }
                }
                else
                {
                    ak = a[k - 1 + (k - 1) * lda] / a[k - 2 + (k - 1) * lda];
                    akm1 = a[k - 2 + (k - 2) * lda] / a[k - 2 + (k - 1) * lda];
                    bk = z[k - 1] / a[k - 2 + (k - 1) * lda];
                    bkm1 = z[k - 2] / a[k - 2 + (k - 1) * lda];
                    denom = ak * akm1 - 1.0;
                    z[k - 1] = (akm1 * bk - bkm1) / denom;
                    z[k - 2] = (ak * bkm1 - bk) / denom;
                }

                k = k - ks;
            }

            s = 1.0 / BLAS1D.dasum(n, z, 1);
            for (i = 1; i <= n; i++)
            {
                z[i - 1] = s * z[i - 1];
            }

            ynorm = s * ynorm;
            //
            //  Solve U' * Z = V.
            //
            k = 1;

            while (k <= n)
            {
                if (kpvt[k - 1] < 0)
                {
                    ks = 2;
                }
                else
                {
                    ks = 1;
                }

                if (k != 1)
                {
                    z[k - 1] = z[k - 1] + BLAS1D.ddot(k - 1, a, 1, z, 1, xIndex: +0 + (k - 1) * lda);
                    if (ks == 2)
                    {
                        z[k] = z[k] + BLAS1D.ddot(k - 1, a, 1, z, 1, xIndex: +0 + k * lda);
                    }

                    kp = Math.Abs(kpvt[k - 1]);

                    if (kp != k)
                    {
                        t = z[k - 1];
                        z[k - 1] = z[kp - 1];
                        z[kp - 1] = t;
                    }
                }

                k = k + ks;
            }

            //
            //  Make ZNORM = 1.0.
            //
            s = 1.0 / BLAS1D.dasum(n, z, 1);
            for (i = 1; i <= n; i++)
            {
                z[i - 1] = s * z[i - 1];
            }

            ynorm = s * ynorm;

            if (anorm != 0.0)
            {
                rcond = ynorm / anorm;
            }
            else
            {
                rcond = 0.0;
            }

            return rcond;
        }

    }
}