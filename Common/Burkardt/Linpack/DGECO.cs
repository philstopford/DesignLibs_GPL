using System;
using Burkardt.BLAS;
using Burkardt.Types;

namespace Burkardt.Linpack;

public static class DGECO
{
    public static double dgeco(ref double[] a, int lda, int n, ref int[] ipvt, ref double[] z )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DGECO factors a real matrix and estimates its condition number.
        //
        //  Discussion:
        //
        //    If RCOND is not needed, DGEFA is slightly faster.
        //
        //    To solve A * X = B, follow DGECO by DGESL.
        //
        //    To compute inverse ( A ) * C, follow DGECO by DGESL.
        //
        //    To compute determinant ( A ), follow DGECO by DGEDI.
        //
        //    To compute inverse ( A ), follow DGECO by DGEDI.
        //
        //    For the system A * X = B, relative perturbations in A and B
        //    of size EPSILON may cause relative perturbations in X of size
        //    EPSILON/RCOND.
        //
        //    If RCOND is so small that the logical expression
        //      1.0D+00 + RCOND == 1.0D+00
        //    is true, then A may be singular to working precision.  In particular,
        //    RCOND is zero if exact singularity is detected or the estimate
        //    underflows.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 May 2005
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
        //    Input/output, double A[LDA*N].  On input, a matrix to be
        //    factored.  On output, the LU factorization of the matrix.
        //
        //    Input, int LDA, the leading dimension of the array A.
        //
        //    Input, int N, the order of the matrix A.
        //
        //    Output, int IPVT[N], the pivot indices.
        //
        //    Output, double Z[N], a work vector whose contents are usually
        //    unimportant.  If A is close to a singular matrix, then Z is an
        //    approximate null vector in the sense that
        //      norm ( A * Z ) = RCOND * norm ( A ) * norm ( Z ).
        //
        //    Output, double DGECO, the value of RCOND, an estimate 
        //    of the reciprocal condition number of A.
        //
    {
        int i;
        int j;
        int k;
        int l;
        double rcond;
        double s;
        double t;
        //
        //  Compute the L1 norm of A.
        //
        double anorm = 0.0;
        for (j = 1; j <= n; j++)
        {
            anorm = Math.Max(anorm, BLAS1D.dasum(n, a, 1, index:  + 0 + (j - 1) * lda));
        }

        //
        //  Compute the LU factorization.
        //
        DGEFA.dgefa(ref a, lda, n, ref ipvt);
        //
        //  RCOND = 1 / ( norm(A) * (estimate of norm(inverse(A))) )
        //
        //  estimate of norm(inverse(A)) = norm(Z) / norm(Y)
        //
        //  where
        //    A * Z = Y
        //  and
        //    A' * Y = E
        //
        //  The components of E are chosen to cause maximum local growth in the
        //  elements of W, where U'*W = E.  The vectors are frequently rescaled
        //  to avoid overflow.
        //
        //  Solve U' * W = E.
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

            if (Math.Abs(a[k - 1 + (k - 1) * lda]) < Math.Abs(ek - z[k - 1]))
            {
                s = Math.Abs(a[k - 1 + (k - 1) * lda]) / Math.Abs(ek - z[k - 1]);
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

            if (a[k - 1 + (k - 1) * lda] != 0.0)
            {
                wk /= a[k - 1 + (k - 1) * lda];
                wkm /= a[k - 1 + (k - 1) * lda];
            }
            else
            {
                wk = 1.0;
                wkm = 1.0;
            }

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
                    for (i = k + 1; i <= n; i++)
                    {
                        z[i - 1] += t * a[k - 1 + (i - 1) * lda];
                    }
                }
            }

            z[k - 1] = wk;
        }

        t = BLAS1D.dasum(n, z, 1);
        for (i = 1; i <= n; i++)
        {
            z[i - 1] /= t;
        }

        //
        //  Solve L' * Y = W
        //
        for (k = n; 1 <= k; k--)
        {
            z[k - 1] += BLAS1D.ddot(n - k, a, 1, z, 1, xIndex:  + k + (k - 1) * lda, yIndex: + k);

            switch (Math.Abs(z[k - 1]))
            {
                case > 1.0:
                {
                    t = Math.Abs(z[k - 1]);
                    for (i = 1; i <= n; i++)
                    {
                        z[i - 1] /= t;
                    }

                    break;
                }
            }

            l = ipvt[k - 1];

            t = z[l - 1];
            z[l - 1] = z[k - 1];
            z[k - 1] = t;
        }

        t = BLAS1D.dasum(n, z, 1);
        for (i = 1; i <= n; i++)
        {
            z[i - 1] /= t;
        }

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

            for (i = k + 1; i <= n; i++)
            {
                z[i - 1] += t * a[i - 1 + (k - 1) * lda];
            }

            switch (Math.Abs(z[k - 1]))
            {
                case > 1.0:
                {
                    ynorm /= Math.Abs(z[k - 1]);
                    t = Math.Abs(z[k - 1]);
                    for (i = 1; i <= n; i++)
                    {
                        z[i - 1] /= t;
                    }

                    break;
                }
            }
        }

        s = BLAS1D.dasum(n, z, 1);
        for (i = 1; i <= n; i++)
        {
            z[i - 1] /= s;
        }

        ynorm /= s;
        //
        //  Solve U * Z = V.
        //
        for (k = n; 1 <= k; k--)
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
                z[k - 1] /= a[k - 1 + (k - 1) * lda];
            }
            else
            {
                z[k - 1] = 1.0;
            }

            for (i = 1; i <= k - 1; i++)
            {
                z[i - 1] -= z[k - 1] * a[i - 1 + (k - 1) * lda];
            }
        }

        //
        //  Normalize Z in the L1 norm.
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