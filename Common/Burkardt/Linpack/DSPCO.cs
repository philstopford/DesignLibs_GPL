using System;
using Burkardt.BLAS;
using Burkardt.Types;

namespace Burkardt.Linpack
{
    public static class DSPCO
    {
        public static double dspco(ref double[] ap, int n, ref int[] kpvt, ref double[] z)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DSPCO factors a real symmetric matrix stored in packed form.
            //
            //  Discussion:
            //
            //    DSPCO uses elimination with symmetric pivoting and estimates
            //    the condition of the matrix.
            //
            //    If RCOND is not needed, DSPFA is slightly faster.
            //
            //    To solve A*X = B, follow DSPCO by DSPSL.
            //
            //    To compute inverse(A)*C, follow DSPCO by DSPSL.
            //
            //    To compute inverse(A), follow DSPCO by DSPDI.
            //
            //    To compute determinant(A), follow DSPCO by DSPDI.
            //
            //    To compute inertia(A), follow DSPCO by DSPDI.
            //
            //  Packed storage:
            //
            //    The following program segment will pack the upper triangle of a
            //    symmetric matrix.
            //
            //      k = 0
            //      do j = 1, n
            //        do i = 1, j
            //          k = k + 1
            //          ap[k-1] = a(i,j)
            //        }
            //      }
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
            //    Input/output, double AP[N*(N+1)/2].  On input, the packed form
            //    of a symmetric matrix A.  The columns of the upper triangle are stored
            //    sequentially in a one-dimensional array.  On output, a block diagonal
            //    matrix and the multipliers which were used to obtain it, stored in
            //    packed form.  The factorization can be written A = U*D*U'
            //    where U is a product of permutation and unit upper triangular
            //    matrices, U' is the transpose of U, and D is block diagonal
            //    with 1 by 1 and 2 by 2 blocks.
            //
            //    Input, int N, the order of the matrix.
            //
            //    Output, int KPVT[N], the pivot indices.
            //
            //    Output, double Z[N] a work vector whose contents are usually
            //    unimportant.  If A is close to a singular matrix, then Z is an
            //    approximate null vector in the sense that
            //      norm(A*Z) = RCOND * norm(A) * norm(Z).
            //
            //    Output, double DSPCO, an estimate of the reciprocal condition number RCOND
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
            int ij;
            int ik;
            int ikm1;
            int ikp1;
            int j;
            int j1;
            int k;
            int kk;
            int km1k;
            int km1km1;
            int kp;
            int kps;
            int ks;
            double rcond;
            double s;
            double t;
            double ynorm;
            //
            //  Find norm of A using only upper half.
            //
            j1 = 1;
            for (j = 1; j <= n; j++)
            {
                z[j - 1] = BLAS1D.dasum(j, ap, 1, index: +j1 - 1);
                ij = j1;
                j1 = j1 + j;
                for (i = 1; i <= j - 1; i++)
                {
                    z[i - 1] = z[i - 1] + Math.Abs(ap[ij - 1]);
                    ij = ij + 1;
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
            DSPFA.dspfa(ref ap, n, ref kpvt);
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
            ik = (n * (n - 1)) / 2;

            while (k != 0)
            {
                kk = ik + k;
                ikm1 = ik - (k - 1);

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
                BLAS1D.daxpy(k - ks, z[k - 1], ap, 1, ref z, 1, xIndex: +ik);

                if (ks != 1)
                {
                    if (z[k - 2] != 0.0)
                    {
                        ek = ek * typeMethods.r8_sign(z[k - 2]);
                    }

                    z[k - 2] = z[k - 2] + ek;
                    BLAS1D.daxpy(k - ks, z[k - 2], ap, 1, ref z, 1, xIndex: +ikm1);
                }

                if (ks != 2)
                {
                    if (Math.Abs(ap[kk - 1]) < Math.Abs(z[k - 1]))
                    {
                        s = Math.Abs(ap[kk - 1]) / Math.Abs(z[k - 1]);
                        for (i = 1; i <= n; i++)
                        {
                            z[i - 1] = s * z[i - 1];
                        }

                        ek = s * ek;
                    }

                    if (ap[kk - 1] != 0.0)
                    {
                        z[k - 1] = z[k - 1] / ap[kk - 1];
                    }
                    else
                    {
                        z[k - 1] = 1.0;
                    }
                }
                else
                {
                    km1k = ik + k - 1;
                    km1km1 = ikm1 + k - 1;
                    ak = ap[kk - 1] / ap[km1k - 1];
                    akm1 = ap[km1km1 - 1] / ap[km1k - 1];
                    bk = z[k - 1] / ap[km1k - 1];
                    bkm1 = z[k - 2] / ap[km1k - 1];
                    denom = ak * akm1 - 1.0;
                    z[k - 1] = (akm1 * bk - bkm1) / denom;
                    z[k - 2] = (ak * bkm1 - bk) / denom;
                }

                k = k - ks;
                ik = ik - k;
                if (ks == 2)
                {
                    ik = ik - (k + 1);
                }
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
            ik = 0;

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
                    z[k - 1] = z[k - 1] + BLAS1D.ddot(k - 1, ap, 1, z, 1, xIndex: +ik);
                    ikp1 = ik + k;

                    if (ks == 2)
                    {
                        z[k] = z[k] + BLAS1D.ddot(k - 1, ap, 1, z, 1, xIndex: +ikp1);
                    }

                    kp = Math.Abs(kpvt[k - 1]);

                    if (kp != k)
                    {
                        t = z[k - 1];
                        z[k - 1] = z[kp - 1];
                        z[kp - 1] = t;
                    }
                }

                ik = ik + k;
                if (ks == 2)
                {
                    ik = ik + (k + 1);
                }

                k = k + ks;
            }

            for (i = 1; i <= n; i++)
            {
                z[i - 1] = s * z[i - 1];
            }

            ynorm = 1.0;
            //
            //  Solve U * D * V = Y.
            //
            k = n;

            ik = (n * (n - 1)) / 2;

            while (0 < k)
            {
                kk = ik + k;
                ikm1 = ik - (k - 1);

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

                    BLAS1D.daxpy(k - ks, z[k - 1], ap, 1, ref z, 1, xIndex: +ik);

                    if (ks == 2)
                    {
                        BLAS1D.daxpy(k - ks, z[k - 2], ap, 1, ref z, 1, xIndex: +ikm1);
                    }
                }

                if (ks != 2)
                {
                    if (Math.Abs(ap[kk - 1]) < Math.Abs(z[k - 1]))
                    {
                        s = Math.Abs(ap[kk - 1]) / Math.Abs(z[k - 1]);
                        for (i = 1; i <= n; i++)
                        {
                            z[i - 1] = s * z[i - 1];
                        }

                        ynorm = s * ynorm;
                    }

                    if (ap[kk - 1] != 0.0)
                    {
                        z[k - 1] = z[k - 1] / ap[kk - 1];
                    }
                    else
                    {
                        z[k - 1] = 1.0;
                    }
                }
                else
                {
                    km1k = ik + k - 1;
                    km1km1 = ikm1 + k - 1;
                    ak = ap[kk - 1] / ap[km1k - 1];
                    akm1 = ap[km1km1 - 1] / ap[km1k - 1];
                    bk = z[k - 1] / ap[km1k - 1];
                    bkm1 = z[k - 2] / ap[km1k - 1];
                    denom = ak * akm1 - 1.0;
                    z[k - 1] = (akm1 * bk - bkm1) / denom;
                    z[k - 2] = (ak * bkm1 - bk) / denom;
                }

                k = k - ks;
                ik = ik - k;
                if (ks == 2)
                {
                    ik = ik - (k + 1);
                }
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
            ik = 0;

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
                    z[k - 1] = z[k - 1] + BLAS1D.ddot(k - 1, ap, 1, z, 1, xIndex: +ik);
                    ikp1 = ik + k;

                    if (ks == 2)
                    {
                        z[k] = z[k] + BLAS1D.ddot(k - 1, ap, 1, z, 1, xIndex: +ikp1);
                    }

                    kp = Math.Abs(kpvt[k - 1]);

                    if (kp != k)
                    {
                        t = z[k - 1];
                        z[k - 1] = z[kp - 1];
                        z[kp - 1] = t;
                    }
                }

                ik = ik + k;
                if (ks == 2)
                {
                    ik = ik + (k + 1);
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