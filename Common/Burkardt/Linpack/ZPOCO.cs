using System;
using System.Numerics;
using Burkardt.BLAS;
using Burkardt.Types;

namespace Burkardt.Linpack
{
    public static class ZPOCO
    {
        public static double zpoco(ref Complex[] a, int lda, int n, ref int info)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ZPOCO factors a complex hermitian positive definite matrix.
            //
            //  Discussion:
            //
            //    The routine also estimates the condition of the matrix.
            //
            //    If RCOND is not needed, ZPOFA is slightly faster.
            //
            //    To solve A*X = B, follow ZPOCO by ZPOSL.
            //
            //    To compute inverse(A)*C, follow ZPOCO by ZPOSL.
            //
            //    To compute determinant(A), follow ZPOCO by ZPODI.
            //
            //    To compute inverse(A), follow ZPOCO by ZPODI.
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
            //    Input/output, Complex A[LDA*N]; on input, the hermitian matrix to be
            //    factored.  On output, an upper triangular matrix R so that
            //      A = hermitian(R)*R
            //    where hermitian(R) is the conjugate transpose.  The strict lower
            //    triangle is unaltered.  If INFO /= 0, the factorization is not complete.
            //
            //    Input, int LDA, the leading dimension of A.
            //
            //    Input, int N, the order of the matrix.
            //
            //    Output, int *INFO.
            //    0, for normal return.
            //    K, signals an error condition.  The leading minor of order K is not
            //    positive definite.
            //
            //    Output, double ZPOCO, an estimate of RCOND, the reciprocal condition of
            //    the matrix.  For the system A*X = B, relative perturbations in A and B
            //    of size EPSILON may cause relative perturbations in X of size
            //    (EPSILON/RCOND).  If RCOND is so small that the logical expression
            //      1.0 + RCOND == 1.0
            //    is true, then A may be singular to working precision.  In particular,
            //    RCOND is zero if exact singularity is detected or the estimate underflows.
            //
            //  Local Parameters:
            //
            //    Workspace, Complex Z[N], a work vector whose contents are usually
            //    unimportant.  If A is close to a singular matrix, then Z is an
            //    approximate null vector in the sense that
            //      norm(A*Z) = RCOND * norm(A) * norm(Z).
            //
        {
            double anorm;
            Complex ek;
            int i;
            int j;
            int k;
            int kp1;
            double rcond;
            double s;
            double sm;
            Complex t;
            Complex wk;
            Complex wkm;
            double ynorm;
            Complex[] z;
            //
            //  Find norm of A using only upper half.
            //
            z = new Complex [n];

            for (j = 1; j <= n; j++)
            {
                z[j - 1] = new Complex(BLAS1Z.dzasum(j, a, 1, index: +0 + (j - 1) * lda), 0.0);
                for (i = 1; i < j; i++)
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
            info = ZPOFA.zpofa(ref a, lda, n);

            if (info != 0)
            {
                rcond = 0.0;
                return rcond;
            }

            //
            //  RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))).
            //
            //  Estimate = norm(Z)/norm(Y) where A*Z = Y and A*Y = E.
            //
            //  The components of E are chosen to cause maximum local
            //  growth in the elements of W where hermitian(R)*W = E.
            //
            //  The vectors are frequently rescaled to avoid overflow.
            //
            //  Solve hermitian(R)*W = E.
            //
            ek = new Complex(1.0, 0.0);
            for (j = 0; j < n; j++)
            {
                z[j] = new Complex(0.0, 0.0);
            }

            for (k = 1; k <= n; k++)
            {
                if (typeMethods.zabs1(z[k - 1]) != 0.0)
                {
                    ek = typeMethods.zsign1(ek, -z[k - 1]);
                }

                if ((a[k - 1 + (k - 1) * lda].Real) < typeMethods.zabs1(ek - z[k - 1]))
                {
                    s = (a[k - 1 + (k - 1) * lda].Real) / typeMethods.zabs1(ek - z[k - 1]);
                    BLAS1Z.zdscal(n, s, ref z, 1);
                    ek = new Complex(s, 0.0) * ek;
                }

                wk = ek - z[k - 1];
                wkm = -ek - z[k - 1];
                s = typeMethods.zabs1(wk);
                sm = typeMethods.zabs1(wkm);
                wk = wk / a[k - 1 + (k - 1) * lda];
                wkm = wkm / a[k - 1 + (k - 1) * lda];
                kp1 = k + 1;

                if (kp1 <= n)
                {
                    for (j = kp1; j <= n; j++)
                    {
                        sm = sm + typeMethods.zabs1(z[j - 1] + wkm * Complex.Conjugate(a[k - 1 + (j - 1) * lda]));
                        z[j - 1] = z[j - 1] + wk * Complex.Conjugate(a[k - 1 + (j - 1) * lda]);
                        s = s + typeMethods.zabs1(z[j - 1]);
                    }

                    if (s < sm)
                    {
                        t = wkm - wk;
                        wk = wkm;
                        for (j = kp1; j <= n; j++)
                        {
                            z[j - 1] = z[j - 1] + t * Complex.Conjugate(a[k - 1 + (j - 1) * lda]);
                        }
                    }
                }

                z[k - 1] = wk;
            }

            s = 1.0 / BLAS1Z.dzasum(n, z, 1);
            BLAS1Z.zdscal(n, s, ref z, 1);
            //
            //  Solve R * Y = W.
            //
            for (k = n; 1 <= k; k--)
            {
                if ((a[k - 1 + (k - 1) * lda].Real) < typeMethods.zabs1(z[k - 1]))
                {
                    s = (a[k - 1 + (k - 1) * lda].Real) / typeMethods.zabs1(z[k - 1]);
                    BLAS1Z.zdscal(n, s, ref z, 1);
                }

                z[k - 1] = z[k - 1] / a[k - 1 + (k - 1) * lda];
                t = -z[k - 1];
                BLAS1Z.zaxpy(k - 1, t, a, 1, ref z, 1, xIndex: +0 + (k - 1) * lda);
            }

            s = 1.0 / BLAS1Z.dzasum(n, z, 1);
            BLAS1Z.zdscal(n, s, ref z, 1);
            ynorm = 1.0;
            //
            //  Solve hermitian(R) * V = Y.
            //
            for (k = 1; k <= n; k++)
            {
                z[k - 1] = z[k - 1] - BLAS1Z.zdotc(k - 1, a, 1, z, 1, xIndex: +0 + (k - 1) * lda);

                if ((a[k - 1 + (k - 1) * lda].Real) < typeMethods.zabs1(z[k - 1]))
                {
                    s = (a[k - 1 + (k - 1) * lda].Real) / typeMethods.zabs1(z[k - 1]);
                    BLAS1Z.zdscal(n, s, ref z, 1);
                    ynorm = s * ynorm;
                }

                z[k - 1] = z[k - 1] / a[k - 1 + (k - 1) * lda];
            }

            s = 1.0 / BLAS1Z.dzasum(n, z, 1);
            BLAS1Z.zdscal(n, s, ref z, 1);
            ynorm = s * ynorm;
            //
            //  Solve R * Z = V.
            //
            for (k = n; 1 <= k; k--)
            {
                if ((a[k - 1 + (k - 1) * lda].Real) < typeMethods.zabs1(z[k - 1]))
                {
                    s = (a[k - 1 + (k - 1) * lda].Real) / typeMethods.zabs1(z[k - 1]);
                    BLAS1Z.zdscal(n, s, ref z, 1);
                    ynorm = s * ynorm;
                }

                z[k - 1] = z[k - 1] / a[k - 1 + (k - 1) * lda];
                t = -z[k - 1];
                BLAS1Z.zaxpy(k - 1, t, a, 1, ref z, 1, xIndex: +0 + (k - 1) * lda);
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