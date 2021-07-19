using System;
using System.Numerics;
using Burkardt.BLAS;
using Burkardt.Types;

namespace Burkardt.Linpack
{
    public static class ZHICO
    {
        public static double zhico(ref Complex[] a, int lda, int n, ref int[] ipvt)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ZHICO factors a complex hermitian matrix and estimates its condition.
            //
            //  Discussion:
            //
            //    If RCOND is not needed, ZHIFA is slightly faster.
            //
            //    To solve A*X = B, follow ZHICO by ZHISL.
            //
            //    To compute inverse(A)*C, follow ZHICO by ZHISL.
            //
            //    To compute inverse(A), follow ZHICO by ZHIDI.
            //
            //    To compute determinant(A), follow ZHICO by ZHIDI.
            //
            //    To compute inertia(A), follow ZHICO by ZHIDI.
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
            //    Input/output, Complex A[LDA*N]; on input, the hermitian matrix
            //    to be factored.  On output, a block diagonal matrix and the multipliers
            //    which were used to obtain it.  The factorization can be written
            //    A = U*D*hermitian(U) where U is a product of permutation and unit
            //    upper triangular matrices, hermitian(U) is the Complex.Conjugateugate transpose
            //    of U, and D is block diagonal with 1 by 1 and 2 by 2 blocks.
            //    Only the diagonal and upper triangle are used.
            //
            //    Input, int LDA, the leading dimension of A.
            //
            //    Input, int N, the order of the matrix.
            //
            //    Output, int IPVT[N], the pivot indices.
            //
            //    Output, double ZHICO, an estimate of RCOND, the reciprocal condition of
            //    the matrix.  For the system A*X = B, relative perturbations in A and B
            //    of size EPSILON may cause relative perturbations in X of size
            //    (EPSILON/RCOND).  If RCOND is so small that the logical expression
            //      1.0 + RCOND == 1.0
            //    is true, then A may be singular to working precision.  In particular,
            //    RCOND is zero if exact singularity is detected or the estimate underflows.
            //
            //  Local Parameter:
            //
            //    Workspace, Complex Z[N], a work vector whose contents are usually
            //    unimportant.  If A is close to a singular matrix, then Z is an
            //    approximate null vector in the sense that
            //      norm(A*Z) = RCOND * norm(A) * norm(Z).
            //
        {
            Complex ak;
            Complex akm1;
            double anorm;
            Complex bk;
            Complex bkm1;
            Complex denom;
            Complex ek;
            int i;
            int j;
            int k;
            int kp;
            int kps;
            int ks;
            double rcond;
            double s;
            Complex t;
            double ynorm;
            Complex[] z;
            //
            //  Find norm of A using only upper half.
            //
            z = new Complex [n];

            for (j = 1; j <= n; j++)
            {
                z[j - 1] = new Complex(BLAS1Z.dzasum(j, a, 1, index: +0 + (j - 1) * lda), 0.0);
                for (i = 1; i <= j - 1; i++)
                {
                    z[i - 1] =
                        new Complex((z[i - 1].Real) + typeMethods.zabs1(a[i - 1 + (j - 1) * lda]), 0.0);
                }
            }

            anorm = 0.0;
            for (j = 0; j < n; j++)
            {
                anorm = Math.Max(anorm, (z[j].Real));
            }

            //
            //  Factor.
            //
            ZHIFA.zhifa(ref a, lda, n, ref ipvt);
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
            //  Solve U*D*W = E.
            //
            ek = new Complex(1.0, 0.0);
            for (i = 0; i < n; i++)
            {
                z[i] = new Complex(0.0, 0.0);
            }

            k = n;

            while (0 < k)
            {
                if (ipvt[k - 1] < 0)
                {
                    ks = 2;
                }
                else
                {
                    ks = 1;
                }

                kp = Math.Abs(ipvt[k - 1]);
                kps = k + 1 - ks;

                if (kp != kps)
                {
                    t = z[kps - 1];
                    z[kps - 1] = z[kp - 1];
                    z[kp - 1] = t;
                }

                if (typeMethods.zabs1(z[k - 1]) != 0.0)
                {
                    ek = typeMethods.zsign1(ek, z[k - 1]);
                }

                z[k - 1] = z[k - 1] + ek;
                BLAS1Z.zaxpy(k - ks, z[k - 1], a, 1, ref z, 1, xIndex: +0 + (k - 1) * lda);

                if (ks != 1)
                {
                    if (typeMethods.zabs1(z[k - 2]) != 0.0)
                    {
                        ek = typeMethods.zsign1(ek, z[k - 2]);
                    }

                    z[k - 2] = z[k - 2] + ek;
                    BLAS1Z.zaxpy(k - ks, z[k - 2], a, 1, ref z, 1, xIndex: +0 + (k - 2) * lda);
                }

                if (ks != 2)
                {
                    if (typeMethods.zabs1(a[k - 1 + (k - 1) * lda]) < typeMethods.zabs1(z[k - 1]))
                    {
                        s = typeMethods.zabs1(a[k - 1 + (k - 1) * lda]) / typeMethods.zabs1(z[k - 1]);
                        BLAS1Z.zdscal(n, s, ref z, 1);
                        ek = new Complex(s, 0.0) * ek;
                    }

                    if (typeMethods.zabs1(a[k - 1 + (k - 1) * lda]) != 0.0)
                    {
                        z[k - 1] = z[k - 1] / a[k - 1 + (k - 1) * lda];
                    }
                    else
                    {
                        z[k - 1] = new Complex(1.0, 0.0);
                    }
                }
                else
                {
                    ak = a[k - 1 + (k - 1) * lda] / Complex.Conjugate(a[k - 2 + (k - 1) * lda]);
                    akm1 = a[k - 2 + (k - 2) * lda] / a[k - 2 + (k - 1) * lda];
                    bk = z[k - 1] / Complex.Conjugate(a[k - 2 + (k - 1) * lda]);
                    bkm1 = z[k - 2] / a[k - 2 + (k - 1) * lda];
                    denom = ak * akm1 - new Complex(1.0, 0.0);
                    z[k - 1] = (akm1 * bk - bkm1) / denom;
                    z[k - 2] = (ak * bkm1 - bk) / denom;
                }

                k = k - ks;
            }

            s = 1.0 / BLAS1Z.dzasum(n, z, 1);
            BLAS1Z.zdscal(n, s, ref z, 1);
            //
            //  Solve hermitian(U) * Y = W.
            //
            k = 1;

            while (k <= n)
            {
                if (ipvt[k - 1] < 0)
                {
                    ks = 2;
                }
                else
                {
                    ks = 1;
                }

                if (k != 1)
                {
                    z[k - 1] = z[k - 1] + BLAS1Z.zdotc(k - 1, a, 1, z, 1, xIndex: +0 + (k - 1) * lda);

                    if (ks == 2)
                    {
                        z[k] = z[k] + BLAS1Z.zdotc(k - 1, a, 1, z, 1, xIndex: +0 + k * lda);
                    }

                    kp = Math.Abs(ipvt[k - 1]);

                    if (kp != k)
                    {
                        t = z[k - 1];
                        z[k - 1] = z[kp - 1];
                        z[kp - 1] = t;
                    }
                }

                k = k + ks;
            }

            s = 1.0 / BLAS1Z.dzasum(n, z, 1);
            BLAS1Z.zdscal(n, s, ref z, 1);
            ynorm = 1.0;
            //
            //  Solve U*D*V = Y.
            //
            k = n;

            while (0 < k)
            {
                if (ipvt[k - 1] < 0)
                {
                    ks = 2;
                }
                else
                {
                    ks = 1;
                }

                if (k != ks)
                {
                    kp = Math.Abs(ipvt[k - 1]);
                    kps = k + 1 - ks;

                    if (kp != kps)
                    {
                        t = z[kps - 1];
                        z[kps - 1] = z[kp - 1];
                        z[kp - 1] = t;
                    }

                    BLAS1Z.zaxpy(k - ks, z[k - 1], a, 1, ref z, 1, xIndex: +0 + (k - 1) * lda);

                    if (ks == 2)
                    {
                        BLAS1Z.zaxpy(k - ks, z[k - 2], a, 1, ref z, 1, xIndex: +0 + (k - 2) * lda);
                    }
                }

                if (ks != 2)
                {
                    if (typeMethods.zabs1(a[k - 1 + (k - 1) * lda]) < typeMethods.zabs1(z[k - 1]))
                    {
                        s = typeMethods.zabs1(a[k - 1 + (k - 1) * lda]) / typeMethods.zabs1(z[k - 1]);
                        BLAS1Z.zdscal(n, s, ref z, 1);
                        ynorm = s * ynorm;
                    }

                    if (typeMethods.zabs1(a[k - 1 + (k - 1) * lda]) != 0.0)
                    {
                        z[k - 1] = z[k - 1] / a[k - 1 + (k - 1) * lda];
                    }
                    else
                    {
                        z[k - 1] = new Complex(1.0, 0.0);
                    }
                }
                else
                {
                    ak = a[k - 1 + (k - 1) * lda] / Complex.Conjugate(a[k - 2 + (k - 1) * lda]);
                    akm1 = a[k - 2 + (k - 2) * lda] / a[k - 2 + (k - 1) * lda];
                    bk = z[k - 1] / Complex.Conjugate(a[k - 2 + (k - 1) * lda]);
                    bkm1 = z[k - 2] / a[k - 2 + (k - 1) * lda];
                    denom = ak * akm1 - new Complex(1.0, 0.0);
                    z[k - 1] = (akm1 * bk - bkm1) / denom;
                    z[k - 2] = (ak * bkm1 - bk) / denom;
                }

                k = k - ks;
            }

            s = 1.0 / BLAS1Z.dzasum(n, z, 1);
            BLAS1Z.zdscal(n, s, ref z, 1);
            ynorm = s * ynorm;
            //
            //  Solve hermitian(U) * Z = V.
            //
            k = 1;

            while (k <= n)
            {
                if (ipvt[k - 1] < 0)
                {
                    ks = 2;
                }
                else
                {
                    ks = 1;
                }

                if (k != 1)
                {
                    z[k - 1] = z[k - 1] + BLAS1Z.zdotc(k - 1, a, 1, z, 1, xIndex: +0 + (k - 1) * lda);

                    if (ks == 2)
                    {
                        z[k] = z[k] + BLAS1Z.zdotc(k - 1, a, 1, z, 1, xIndex: +0 + k * lda);
                    }

                    kp = Math.Abs(ipvt[k - 1]);

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
            //  Make ZNORM = 1.
            //
            s = 1.0 / BLAS1Z.dzasum(n, z, 1);
            BLAS1Z.zdscal(n, s, ref z, 1);
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