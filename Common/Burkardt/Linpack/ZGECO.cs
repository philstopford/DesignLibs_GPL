using System;
using System.Numerics;
using Burkardt.BLAS;
using Burkardt.Types;

namespace Burkardt.Linpack;

public static class ZGECO
{
    public static double zgeco(ref Complex[] a, int lda, int n, ref int[] ipvt)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZGECO factors a complex matrix and estimates its condition.
        //
        //  Discussion:
        //
        //    If RCOND is not needed, ZGEFA is slightly faster.
        //
        //    To solve A*X = B, follow ZGECO by ZGESL.
        //
        //    To compute inverse(A)*C, follow ZGECO by ZGESL.
        //
        //    To compute determinant(A), follow ZGECO by ZGEDI.
        //
        //    To compute inverse(A), follow ZGECO by ZGEDI.
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
        //    Input/output, Complex A[LDA*N], on input, the matrix to be
        //    factored.  On output, an upper triangular matrix and the multipliers
        //    used to obtain it.  The factorization can be written A = L*U where
        //    L is a product of permutation and unit lower triangular matrices
        //    and U is upper triangular.
        //
        //    Input, int LDA, the leading dimension of A.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Output, int IPVT[N], the pivot indices.
        //
        //    Output, double SGECO, an estimate of the reciprocal condition of A.
        //    For the system A*X = B, relative perturbations in A and B of size
        //    EPSILON may cause relative perturbations in X of size (EPSILON/RCOND).
        //    If RCOND is so small that the logical expression
        //      1.0 + RCOND == 1.0
        //    is true, then A may be singular to working precision.  In particular,
        //    RCOND is zero if exact singularity is detected or the estimate
        //    underflows.
        //
        //  Local Parameters:
        //
        //    Workspace, Complex Z[N], a work vector whose contents are usually
        //    unimportant.  If A is close to a singular matrix, then Z is
        //    an approximate null vector in the sense that
        //      norm ( A * Z ) = RCOND * norm ( A ) * norm ( Z ).
        //
    {
        int i;
        int j;
        int k;
        int l;
        double rcond;
        double s;
        Complex t;

        Complex[] z = new Complex [n];
        //
        //  Compute the 1-norm of A.
        //
        double anorm = 0.0;
        for (j = 0; j < n; j++)
        {
            anorm = Math.Max(anorm, BLAS1Z.dzasum(n, a, 1, index: +0 + j * lda));
        }

        //
        //  Factor.
        //
        ZGEFA.zgefa(ref a, lda, n, ref ipvt);
        //
        //  RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))).
        //
        //  Estimate = norm(Z)/norm(Y) where A*Z = Y and hermitian(A)*Y = E.
        //
        //  Hermitian(A) is the Complex.Conjugateugate transpose of A.
        //
        //  The components of E are chosen to cause maximum local
        //  growth in the elements of W where hermitian(U)*W = E.
        //
        //  The vectors are frequently rescaled to avoid overflow.
        //
        //  Solve hermitian(U)*W = E.
        //
        Complex ek = new(1.0, 0.0);
        for (i = 0; i < n; i++)
        {
            z[i] = new Complex(0.0, 0.0);
        }

        for (k = 1; k <= n; k++)
        {
            if (typeMethods.zabs1(z[k - 1]) != 0.0)
            {
                ek = typeMethods.zsign1(ek, -z[k - 1]);
            }

            if (typeMethods.zabs1(a[k - 1 + (k - 1) * lda]) < typeMethods.zabs1(ek - z[k - 1]))
            {
                s = typeMethods.zabs1(a[k - 1 + (k - 1) * lda]) / typeMethods.zabs1(ek - z[k - 1]);
                BLAS1Z.zdscal(n, s, ref z, 1);
                ek = new Complex(s, 0.0) * ek;
            }

            Complex wk = ek - z[k - 1];
            Complex wkm = -ek - z[k - 1];
            s = typeMethods.zabs1(wk);
            double sm = typeMethods.zabs1(wkm);

            if (typeMethods.zabs1(a[k - 1 + (k - 1) * lda]) != 0.0)
            {
                wk /= Complex.Conjugate(a[k - 1 + (k - 1) * lda]);
                wkm /= Complex.Conjugate(a[k - 1 + (k - 1) * lda]);
            }
            else
            {
                wk = new Complex(1.0, 0.0);
                wkm = new Complex(1.0, 0.0);
            }

            for (j = k + 1; j <= n; j++)
            {
                sm += typeMethods.zabs1(z[j - 1] + wkm * Complex.Conjugate(a[k - 1 + (j - 1) * lda]));
                z[j - 1] += wk * Complex.Conjugate(a[k - 1 + (j - 1) * lda]);
                s += typeMethods.zabs1(z[j - 1]);
            }

            if (s < sm)
            {
                t = wkm - wk;
                wk = wkm;
                for (j = k + 1; j <= n; j++)
                {
                    z[j - 1] += t * Complex.Conjugate(a[k - 1 + (j - 1) * lda]);
                }
            }

            z[k - 1] = wk;
        }

        s = 1.0 / BLAS1Z.dzasum(n, z, 1);
        BLAS1Z.zdscal(n, s, ref z, 1);
        //
        //  Solve hermitian(L) * Y = W.
        //
        for (k = n; 1 <= k; k--)
        {
            if (k < n)
            {
                z[k - 1] += BLAS1Z.zdotc(n - k, a, 1, z, 1, xIndex: +k + (k - 1) * lda, yIndex: +k);
            }

            if (1.0 < typeMethods.zabs1(z[k - 1]))
            {
                s = 1.0 / typeMethods.zabs1(z[k - 1]);
                BLAS1Z.zdscal(n, s, ref z, 1);
            }

            l = ipvt[k - 1];

            t = z[l - 1];
            z[l - 1] = z[k - 1];
            z[k - 1] = t;
        }

        s = 1.0 / BLAS1Z.dzasum(n, z, 1);
        BLAS1Z.zdscal(n, s, ref z, 1);

        double ynorm = 1.0;
        //
        //  Solve L * V = Y.
        //
        for (k = 1; k <= n; k++)
        {
            l = ipvt[k - 1];

            t = z[l - 1];
            z[l - 1] = z[k - 1];
            z[k - 1] = t;

            if (k < n)
            {
                BLAS1Z.zaxpy(n - k, t, a, 1, ref z, 1, xIndex: +k + (k - 1) * lda, yIndex: +k);
            }

            if (!(1.0 < typeMethods.zabs1(z[k - 1])))
            {
                continue;
            }

            s = 1.0 / typeMethods.zabs1(z[k - 1]);
            BLAS1Z.zdscal(n, s, ref z, 1);
            ynorm = s * ynorm;
        }

        s = 1.0 / BLAS1Z.dzasum(n, z, 1);
        BLAS1Z.zdscal(n, s, ref z, 1);
        ynorm = s * ynorm;
        //
        //  Solve U * Z = V.
        //
        for (k = n; 1 <= k; k--)
        {
            if (typeMethods.zabs1(a[k - 1 + (k - 1) * lda]) < typeMethods.zabs1(z[k - 1]))
            {
                s = typeMethods.zabs1(a[k - 1 + (k - 1) * lda]) / typeMethods.zabs1(z[k - 1]);
                BLAS1Z.zdscal(n, s, ref z, 1);
                ynorm = s * ynorm;
            }

            if (typeMethods.zabs1(a[k - 1 + (k - 1) * lda]) != 0.0)
            {
                z[k - 1] /= a[k - 1 + (k - 1) * lda];
            }
            else
            {
                z[k - 1] = new Complex(1.0, 0.0);
            }

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