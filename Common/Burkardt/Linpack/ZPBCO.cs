using System;
using System.Numerics;
using Burkardt.BLAS;
using Burkardt.Types;

namespace Burkardt.Linpack;

public static class ZPBCO
{
    public static double zpbco(ref Complex[] abd, int lda, int n, int m, ref int info)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZPBCO factors a Complex hermitian positive definite band matrix.
        //
        //  Discussion:
        //
        //    The routine also estimates the condition number of the matrix.
        //
        //    If RCOND is not needed, ZPBFA is slightly faster.
        //
        //    To solve A*X = B, follow ZPBCO by ZPBSL.
        //
        //    To compute inverse(A)*C, follow ZPBCO by ZPBSL.
        //
        //    To compute determinant(A), follow ZPBCO by ZPBDI.
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
        //    Input/output, Complex ABD[LDA*N]; on input, the matrix to be factored.
        //    The columns of the upper triangle are stored in the columns of ABD,
        //    and the diagonals of the upper triangle are stored in the rows of ABD.
        //    On output, an upper triangular matrix R, stored in band form, so that
        //    A = hermitian(R) * R.  If INFO != 0, the factorization is not complete.
        //
        //    Input, int LDA, the leading dimension of ABD.
        //    LDA must be at least M+1.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int M, the number of diagonals above the main diagonal.
        //    0 <= M < N.
        //
        //    Output, double ZPBCO, an estimate of RCOND, the reciprocal condition of
        //    the matrix.  For the system A*X = B, relative perturbations in A and B
        //    of size EPSILON may cause relative perturbations in X of size
        //    (EPSILON/RCOND).  If RCOND is so small that the logical expression
        //      1.0 + RCOND == 1.0
        //    is true, then A may be singular to working precision.  In particular,
        //    RCOND is zero if exact singularity is detected or the estimate underflows.
        //
        //    Output, int *INFO.
        //    0, for normal return.
        //    K, signals an error condition.  The leading minor of order K is not
        //    positive definite.
        //
        //  Local Parameter:
        //
        //    Workspace, Complex Z[N], a work vector whose contents are usually
        //    unimportant.  If A is singular to working precision, then Z is
        //    an approximate null vector in the sense that
        //    norm ( A * Z ) = RCOND * norm ( A ) * norm ( Z ).
        //    If INFO != 0, Z is unchanged.
        //
    {
        int i;
        int j;
        int k;
        int la;
        int lb;
        int lm;
        double rcond;
        double s;
        Complex t;
        //
        //  Find the norm of A.
        //
        Complex[] z = new Complex [n];

        for (j = 1; j <= n; j++)
        {
            int l = Math.Min(j, m + 1);
            int mu = Math.Max(m + 2 - j, 1);
            z[j - 1] = new Complex(BLAS1Z.dzasum(l, abd, 1, index: +mu - 1 + (j - 1) * lda), 0.0);
            k = j - l;

            for (i = mu; i <= m; i++)
            {
                k += 1;
                z[k - 1] = new Complex(z[k - 1].Real
                                       + typeMethods.zabs1(abd[i - 1 + (j - 1) * lda]), 0.0);
            }
        }

        double anorm = 0.0;
        for (j = 0; j < n; j++)
        {
            anorm = Math.Max(anorm, z[j].Real);
        }

        //
        //  Factor.
        //
        info = ZPBFA.zpbfa(ref abd, lda, n, m);

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
        Complex ek = new Complex(1.0, 0.0);

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

            if (abd[m + (k - 1) * lda].Real < typeMethods.zabs1(ek - z[k - 1]))
            {
                s = abd[m + (k - 1) * lda].Real / typeMethods.zabs1(ek - z[k - 1]);
                BLAS1Z.zdscal(n, s, ref z, 1);
                ek = new Complex(s, 0.0) * ek;
            }

            Complex wk = ek - z[k - 1];
            Complex wkm = -ek - z[k - 1];
            s = typeMethods.zabs1(wk);
            double sm = typeMethods.zabs1(wkm);
            wk /= abd[m + (k - 1) * lda];
            wkm /= abd[m + (k - 1) * lda];
            int j2 = Math.Min(k + m, n);
            i = m + 1;

            if (k + 1 <= j2)
            {
                for (j = k + 1; j <= j2; j++)
                {
                    i -= 1;
                    sm += typeMethods.zabs1(z[j - 1] + wkm * Complex.Conjugate(abd[i - 1 + (j - 1) * lda]));
                    z[j - 1] += wk * Complex.Conjugate(abd[i - 1 + (j - 1) * lda]);
                    s += typeMethods.zabs1(z[j - 1]);
                }

                if (s < sm)
                {
                    t = wkm - wk;
                    wk = wkm;
                    i = m + 1;
                    for (j = k + 1; j <= j2; j++)
                    {
                        i -= 1;
                        z[j - 1] += t * Complex.Conjugate(abd[i - 1 + (j - 1) * lda]);
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
            if (abd[m + (k - 1) * lda].Real < typeMethods.zabs1(z[k - 1]))
            {
                s = abd[m + (k - 1) * lda].Real / typeMethods.zabs1(z[k - 1]);
                BLAS1Z.zdscal(n, s, ref z, 1);
            }

            z[k - 1] /= abd[m + (k - 1) * lda];
            lm = Math.Min(k - 1, m);
            la = m + 1 - lm;
            lb = k - lm;
            t = -z[k - 1];
            BLAS1Z.zaxpy(lm, t, abd, 1, ref z, 1, xIndex: +la - 1 + (k - 1) * lda, yIndex: +lb - 1);
        }

        s = 1.0 / BLAS1Z.dzasum(n, z, 1);
        BLAS1Z.zdscal(n, s, ref z, 1);
        double ynorm = 1.0;
        //
        //  Solve hermitian(R)*V = Y.
        //
        for (k = 1; k <= n; k++)
        {
            lm = Math.Min(k - 1, m);
            la = m + 1 - lm;
            lb = k - lm;
            z[k - 1] -= BLAS1Z.zdotc(lm, abd, 1, z, 1, xIndex: +la - 1 + (k - 1) * lda, yIndex: +lb - 1);

            if (abd[m + (k - 1) * lda].Real < typeMethods.zabs1(z[k - 1]))
            {
                s = abd[m + (k - 1) * lda].Real / typeMethods.zabs1(z[k - 1]);
                BLAS1Z.zdscal(n, s, ref z, 1);
                ynorm = s * ynorm;
            }

            z[k - 1] /= abd[m + (k - 1) * lda];
        }

        s = 1.0 / BLAS1Z.dzasum(n, z, 1);
        BLAS1Z.zdscal(n, s, ref z, 1);
        ynorm = s * ynorm;
        //
        //  Solve R * Z = W.
        //
        for (k = n; 1 <= k; k--)
        {
            if (abd[m + (k - 1) * lda].Real < typeMethods.zabs1(z[k - 1]))
            {
                s = abd[m + (k - 1) * lda].Real / typeMethods.zabs1(z[k - 1]);
                BLAS1Z.zdscal(n, s, ref z, 1);
                ynorm = s * ynorm;
            }

            z[k - 1] /= abd[m + (k - 1) * lda];
            lm = Math.Min(k - 1, m);
            la = m + 1 - lm;
            lb = k - lm;
            t = -z[k - 1];
            BLAS1Z.zaxpy(lm, t, abd, 1, ref z, 1, xIndex: +la - 1 + (k - 1) * lda, yIndex: +lb - 1);
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