﻿using System;
using Burkardt.BLAS;
using Burkardt.Types;

namespace Burkardt.Linpack;

public static class DPOCO
{
    public static double dpoco(ref double[] a, int lda, int n, ref double[] z )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DPOCO factors a real symmetric positive definite matrix and estimates its condition.
        //
        //  Discussion:
        //
        //    If RCOND is not needed, DPOFA is slightly faster.
        //
        //    To solve A*X = B, follow DPOCO by DPOSL.
        //
        //    To compute inverse(A)*C, follow DPOCO by DPOSL.
        //
        //    To compute determinant(A), follow DPOCO by DPODI.
        //
        //    To compute inverse(A), follow DPOCO by DPODI.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 June 2005
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
        //    On output, an upper triangular matrix R so that A = R'*R where R'
        //    is the transpose.  The strict lower triangle is unaltered.
        //    If INFO /= 0, the factorization is not complete.
        //
        //    Input, int LDA, the leading dimension of the array A.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Output, double Z[N], a work vector whose contents are usually
        //    unimportant.  If A is close to a singular matrix, then Z is an
        //    approximate null vector in the sense that
        //      norm(A*Z) = RCOND * norm(A) * norm(Z).
        //    If INFO /= 0, Z is unchanged.
        //
        //    Output, double DPOCO, an estimate of the reciprocal
        //    condition of A.  For the system A*X = B, relative perturbations in
        //    A and B of size EPSILON may cause relative perturbations in X of
        //    size EPSILON/RCOND.  If RCOND is so small that the logical expression
        //      1.0D+00 + RCOND == 1.0D+00
        //    is true, then A may be singular to working precision.  In particular,
        //    RCOND is zero if exact singularity is detected or the estimate underflows.
        //
    {
        int i;
        int j;
        int k;
        double rcond;
        double s;
        double t;
        //
        //  Find norm of A using only upper half.
        //
        for (j = 1; j <= n; j++)
        {
            z[j - 1] = BLAS1D.dasum(j, a, 1, index: + 0 + (j - 1) * lda);
            for (i = 1; i <= j - 1; i++)
            {
                z[i - 1] += Math.Abs(a[i - 1 + (j - 1) * lda]);
            }
        }

        double anorm = 0.0;
        for (i = 1; i <= n; i++)
        {
            anorm = Math.Max(anorm, z[i - 1]);
        }

        //
        //  Factor.
        //
        int info = DPOFA.dpofa(ref a, lda, n);

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
        double ek = 1.0;
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

            if (a[k - 1 + (k - 1) * lda] < Math.Abs(ek - z[k - 1]))
            {
                s = a[k - 1 + (k - 1) * lda] / Math.Abs(ek - z[k - 1]);
                for (i = 1; i <= n; i++)
                {
                    z[i - 1] = s * z[i - 1];
                }

                ek = s * ek;
            }

            double wk = ek - z[k - 1];
            double wkm = -ek - z[k - 1];
            s = Math.Abs(wk);
            double sm = Math.Abs(wkm);
            wk /= a[k - 1 + (k - 1) * lda];
            wkm /= a[k - 1 + (k - 1) * lda];

            if (k + 1 <= n)
            {
                for (j = k + 1; j <= n; j++)
                {
                    sm += Math.Abs(z[j - 1] + wkm * a[k - 1 + (j - 1) * lda]);
                    z[j - 1] += wk * a[k - 1 + (j - 1) * lda];
                    s += Math.Abs(z[j - 1]);
                }

                if (s < sm)
                {
                    t = wkm - wk;
                    wk = wkm;
                    for (j = k + 1; j <= n; j++)
                    {
                        z[j - 1] += t * a[k - 1 + (j - 1) * lda];
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
            if (a[k - 1 + (k - 1) * lda] < Math.Abs(z[k - 1]))
            {
                s = a[k - 1 + (k - 1) * lda] / Math.Abs(z[k - 1]);
                for (i = 1; i <= n; i++)
                {
                    z[i - 1] = s * z[i - 1];
                }
            }

            z[k - 1] /= a[k - 1 + (k - 1) * lda];
            t = -z[k - 1];
            BLAS1D.daxpy(k - 1, t, a, 1, ref z, 1, xIndex: + 0 + (k - 1) * lda);
        }

        s = BLAS1D.dasum(n, z, 1);
        for (i = 1; i <= n; i++)
        {
            z[i - 1] /= s;
        }

        double ynorm = 1.0;
        //
        //  Solve R' * V = Y.
        //
        for (k = 1; k <= n; k++)
        {
            z[k - 1] -= BLAS1D.ddot(k - 1, a, 1, z, 1, xIndex: + 0 + (k - 1) * lda);

            if (a[k - 1 + (k - 1) * lda] < Math.Abs(z[k - 1]))
            {
                s = a[k - 1 + (k - 1) * lda] / Math.Abs(z[k - 1]);
                for (i = 1; i <= n; i++)
                {
                    z[i - 1] = s * z[i - 1];
                }

                ynorm = s * ynorm;
            }

            z[k - 1] /= a[k - 1 + (k - 1) * lda];
        }

        s = 1.0 / BLAS1D.dasum(n, z, 1);
        for (i = 1; i <= n; i++)
        {
            z[i - 1] = s * z[i - 1];
        }

        ynorm = s * ynorm;
        //
        //  Solve R * Z = V.
        //
        for (k = n; 1 <= k; k--)
        {
            if (a[k - 1 + (k - 1) * lda] < Math.Abs(z[k - 1]))
            {
                s = a[k - 1 + (k - 1) * lda] / Math.Abs(z[k - 1]);
                for (i = 1; i <= n; i++)
                {
                    z[i - 1] = s * z[i - 1];
                }

                ynorm = s * ynorm;
            }

            z[k - 1] /= a[k - 1 + (k - 1) * lda];
            t = -z[k - 1];
            BLAS1D.daxpy(k - 1, t, a, 1, ref z, 1, xIndex: + 0 + (k - 1) * lda);
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