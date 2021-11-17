using System;
using Burkardt.BLAS;
using Burkardt.Types;

namespace Burkardt.Linpack;

public static class DPBCO
{
    public static double dpbco(ref double[] abd, int lda, int n, int m, ref double[] z )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DPBCO factors a real symmetric positive definite banded matrix.
        //
        //  Discussion:
        //
        //    DPBCO also estimates the condition of the matrix.
        //
        //    If RCOND is not needed, DPBFA is slightly faster.
        //
        //    To solve A*X = B, follow DPBCO by DPBSL.
        //
        //    To compute inverse(A)*C, follow DPBCO by DPBSL.
        //
        //    To compute determinant(A), follow DPBCO by DPBDI.
        //
        //  Band storage:
        //
        //    If A is a symmetric positive definite band matrix, the following
        //    program segment will set up the input.
        //
        //      m = (band width above diagonal)
        //      do j = 1, n
        //        i1 = max (1, j-m)
        //        do i = i1, j
        //          k = i-j+m+1
        //          abd(k,j) = a(i,j)
        //        }
        //      }
        //
        //    This uses M + 1 rows of A, except for the M by M upper left triangle,
        //    which is ignored.
        //
        //    For example, if the original matrix is
        //
        //      11 12 13  0  0  0
        //      12 22 23 24  0  0
        //      13 23 33 34 35  0
        //       0 24 34 44 45 46
        //       0  0 35 45 55 56
        //       0  0  0 46 56 66
        //
        //    then N = 6, M = 2  and ABD should contain
        //
        //       *  * 13 24 35 46
        //       * 12 23 34 45 56
        //      11 22 33 44 55 66
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
        //    Input/output, double ABD[LDA*N].  On input, the matrix to be
        //    factored.  The columns of the upper triangle are stored in the columns
        //    of ABD and the diagonals of the upper triangle are stored in the rows
        //    of ABD.  On output, an upper triangular matrix R, stored in band form,
        //    so that A = R'*R.  If INFO /= 0, the factorization is not complete.
        //
        //    Input, int LDA, the leading dimension of the array ABD.
        //    M+1 <= LDA is required.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int M, the number of diagonals above the main diagonal.
        //
        //    Output, double Z[N], a work vector whose contents are usually
        //    unimportant.  If A is singular to working precision, then Z is an
        //    approximate null vector in the sense that
        //      norm(A*Z) = RCOND * norm(A) * norm(Z).
        //    If INFO /= 0, Z is unchanged.
        //
        //    Output, double DPBCO, an estimate of the reciprocal condition number
        //    RCOND.  For the system A*X = B, relative perturbations in A and B of size
        //    EPSILON may cause relative perturbations in X of size EPSILON/RCOND.
        //    If RCOND is so small that the logical expression
        //      1.0 + RCOND == 1.0D+00
        //    is true, then A may be singular to working precision.  In particular,
        //    RCOND is zero if exact singularity is detected or the estimate underflows.
        //
    {
        double anorm;
        double ek;
        int i;
        int info;
        int j;
        int j2;
        int k;
        int l;
        int la;
        int lb;
        int lm;
        int mu;
        double rcond;
        double s;
        double sm;
        double t;
        double wk;
        double wkm;
        double ynorm;
        //
        //  Find the norm of A.
        //
        for (j = 1; j <= n; j++)
        {
            l = Math.Min(j, m + 1);
            mu = Math.Max(m + 2 - j, 1);
            z[j - 1] = BLAS1D.dasum(l, abd, 1, index: + mu - 1 + (j - 1) * lda);
            k = j - l;
            for (i = mu; i <= m; i++)
            {
                k += 1;
                z[k - 1] += Math.Abs(abd[i - 1 + (j - 1) * lda]);
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
        info = DPBFA.dpbfa(ref abd, lda, n, m);

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
        //  growth in the elements of W where R'*W = E.
        //
        //  The vectors are frequently rescaled to avoid overflow.
        //
        //  Solve R' * W = E.
        //
        ek = 1.0;
        for (i = 1; i <= n; i++)
        {
            z[i - 1] = 0.0;
        }

        for (k = 1; k <= n; k++)
        {
            if (z[k - 1] != 0.0)
            {
                ek *= typeMethods.r8_sign(-z[k - 1]);
            }

            if (abd[m + (k - 1) * lda] < Math.Abs(ek - z[k - 1]))
            {
                s = abd[m + (k - 1) * lda] / Math.Abs(ek - z[k - 1]);
                for (i = 1; i <= n; i++)
                {
                    z[i - 1] = s * z[i - 1];
                }

                ek = s * ek;
            }

            wk = ek - z[k - 1];
            wkm = -ek - z[k - 1];
            s = Math.Abs(wk);
            sm = Math.Abs(wkm);
            wk /= abd[m + (k - 1) * lda];
            wkm /= abd[m + (k - 1) * lda];
            j2 = Math.Min(k + m, n);
            i = m + 1;

            if (k + 1 <= j2)
            {
                for (j = k + 1; j <= j2; j++)
                {
                    i -= 1;
                    sm += Math.Abs(z[j - 1] + wkm * abd[i - 1 + (j - 1) * lda]);
                    z[j - 1] += wk * abd[i - 1 + (j - 1) * lda];
                    s += Math.Abs(z[j - 1]);
                }

                if (s < sm)
                {
                    t = wkm - wk;
                    wk = wkm;
                    i = m + 1;

                    for (j = k + 1; j <= j2; j++)
                    {
                        i -= 1;
                        z[j - 1] += t * abd[i - 1 + (j - 1) * lda];
                    }
                }
            }

            z[k - 1] = wk;
        }

        s = BLAS1D.dasum(n, z, 1);
        for (i = 1; i <= n; i++)
        {
            z[i - 1] /= s;
        }

        //
        //  Solve R * Y = W.
        //
        for (k = n; 1 <= k; k--)
        {
            if (abd[m + (k - 1) * lda] < Math.Abs(z[k - 1]))
            {
                s = abd[m + (k - 1) * lda] / Math.Abs(z[k - 1]);
                for (i = 1; i <= n; i++)
                {
                    z[i - 1] = s * z[i - 1];
                }
            }

            z[k - 1] /= abd[m + (k - 1) * lda];
            lm = Math.Min(k - 1, m);
            la = m + 1 - lm;
            lb = k - lm;
            t = -z[k - 1];
            BLAS1D.daxpy(lm, t, abd, 1, ref z, 1, xIndex: + la - 1 + (k - 1) * lda, yIndex: + lb - 1);
        }

        s = BLAS1D.dasum(n, z, 1);
        for (i = 1; i <= n; i++)
        {
            z[i - 1] /= s;
        }

        ynorm = 1.0;
        //
        //  Solve R' * V = Y.
        //
        for (k = 1; k <= n; k++)
        {
            lm = Math.Min(k - 1, m);
            la = m + 1 - lm;
            lb = k - lm;

            z[k - 1] -= BLAS1D.ddot(lm, abd, 1, z, 1, xIndex:  + la - 1 + (k - 1) * lda, yIndex:  + lb - 1);

            if (abd[m + (k - 1) * lda] < Math.Abs(z[k - 1]))
            {
                s = abd[m + (k - 1) * lda] / Math.Abs(z[k - 1]);
                for (i = 1; i <= n; i++)
                {
                    z[i - 1] = s * z[i - 1];
                }

                ynorm = s * ynorm;
            }

            z[k - 1] /= abd[m + (k - 1) * lda];
        }

        s = 1.0 / BLAS1D.dasum(n, z, 1);
        for (i = 1; i <= n; i++)
        {
            z[i - 1] = s * z[i - 1];
        }

        ynorm = s * ynorm;
        //
        //  Solve R * Z = W.
        //
        for (k = n; 1 <= k; k--)
        {
            if (abd[m + (k - 1) * lda] < Math.Abs(z[k - 1]))
            {
                s = abd[m + (k - 1) * lda] / Math.Abs(z[k - 1]);
                for (i = 1; i <= n; i++)
                {
                    z[i - 1] = s * z[i - 1];
                }

                ynorm = s * ynorm;
            }

            z[k - 1] /= abd[m + (k - 1) * lda];
            lm = Math.Min(k - 1, m);
            la = m + 1 - lm;
            lb = k - lm;
            t = -z[k - 1];
            BLAS1D.daxpy(lm, t, abd, 1, ref z, 1, xIndex:  + la - 1 + (k - 1) * lda, yIndex:  + lb - 1);
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