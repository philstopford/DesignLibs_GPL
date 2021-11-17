using System;
using Burkardt.BLAS;
using Burkardt.Types;

namespace Burkardt.Linpack;

public static class DPPCO
{
    public static double dppco(ref double[] ap, int n, ref double[] z)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DPPCO factors a real symmetric positive definite matrix in packed form.
        //
        //  Discussion:
        //
        //    DPPCO also estimates the condition of the matrix.
        //
        //    If RCOND is not needed, DPPFA is slightly faster.
        //
        //    To solve A*X = B, follow DPPCO by DPPSL.
        //
        //    To compute inverse(A)*C, follow DPPCO by DPPSL.
        //
        //    To compute determinant(A), follow DPPCO by DPPDI.
        //
        //    To compute inverse(A), follow DPPCO by DPPDI.
        //
        //  Packed storage:
        //
        //    The following program segment will pack the upper triangle of
        //    a symmetric matrix.
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
        //    Input/output, double AP[N*(N+1)/2].  On input, the packed
        //    form of a symmetric matrix A.  The columns of the upper triangle are
        //    stored sequentially in a one-dimensional array.  On output, an upper
        //    triangular matrix R, stored in packed form, so that A = R'*R.
        //    If INFO /= 0, the factorization is not complete.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Output, double Z[N], a work vector whose contents are usually
        //    unimportant.  If A is singular to working precision, then Z is an
        //    approximate null vector in the sense that
        //      norm(A*Z) = RCOND * norm(A) * norm(Z).
        //    If INFO /= 0, Z is unchanged.
        //
        //    Output, double DPPCO, an estimate of the reciprocal condition number RCOND
        //    of A.  For the system A*X = B, relative perturbations in A and B of size
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
        int ij;
        int info;
        int j;
        int j1;
        int k;
        int kj;
        int kk;
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
        j1 = 1;
        for (j = 1; j <= n; j++)
        {
            z[j - 1] = BLAS1D.dasum(j, ap, 1, index: +j1 - 1);
            ij = j1;
            j1 += j;
            for (i = 1; i <= j - 1; i++)
            {
                z[i - 1] += Math.Abs(ap[ij - 1]);
                ij += 1;
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
        info = DPPFA.dppfa(ref ap, n);

        if (info != 0)
        {
            rcond = 0.0;
            return rcond;
        }

        //
        //  RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))).
        //
        //  Estimate = norm(Z)/norm(Y) where A * Z = Y and A * Y = E.
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

        kk = 0;

        for (k = 1; k <= n; k++)
        {
            kk += k;

            if (z[k - 1] != 0.0)
            {
                ek *= typeMethods.r8_sign(-z[k - 1]);
            }

            if (ap[kk - 1] < Math.Abs(ek - z[k - 1]))
            {
                s = ap[kk - 1] / Math.Abs(ek - z[k - 1]);
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
            wk /= ap[kk - 1];
            wkm /= ap[kk - 1];
            kj = kk + k;

            if (k + 1 <= n)
            {
                for (j = k + 1; j <= n; j++)
                {
                    sm += Math.Abs(z[j - 1] + wkm * ap[kj - 1]);
                    z[j - 1] += wk * ap[kj - 1];
                    s += Math.Abs(z[j - 1]);
                    kj += j;
                }

                if (s < sm)
                {
                    t = wkm - wk;
                    wk = wkm;
                    kj = kk + k;

                    for (j = k + 1; j <= n; j++)
                    {
                        z[j - 1] += t * ap[kj - 1];
                        kj += j;
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
            if (ap[kk - 1] < Math.Abs(z[k - 1]))
            {
                s = ap[kk - 1] / Math.Abs(z[k - 1]);
                for (i = 1; i <= n; i++)
                {
                    z[i - 1] = s * z[i - 1];
                }
            }

            z[k - 1] /= ap[kk - 1];
            kk -= k;
            t = -z[k - 1];
            BLAS1D.daxpy(k - 1, t, ap, 1, ref z, 1, xIndex: +kk);
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
            z[k - 1] -= BLAS1D.ddot(k - 1, ap, 1, z, 1, xIndex: +kk);
            kk += k;

            if (ap[kk - 1] < Math.Abs(z[k - 1]))
            {
                s = ap[kk - 1] / Math.Abs(z[k - 1]);
                for (i = 1; i <= n; i++)
                {
                    z[i - 1] = s * z[i - 1];
                }

                ynorm = s * ynorm;
            }

            z[k - 1] /= ap[kk - 1];
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
            if (ap[kk - 1] < Math.Abs(z[k - 1]))
            {
                s = ap[kk - 1] / Math.Abs(z[k - 1]);
                for (i = 1; i <= n; i++)
                {
                    z[i - 1] = s * z[i - 1];
                }

                ynorm = s * ynorm;
            }

            z[k - 1] /= ap[kk - 1];
            kk -= k;
            t = -z[k - 1];
            BLAS1D.daxpy(k - 1, t, ap, 1, ref z, 1, xIndex: +kk);
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