using System;
using System.Numerics;
using Burkardt.BLAS;
using Burkardt.Types;

namespace Burkardt.Linpack;

public static class ZPPCO
{
    public static double zppco(ref Complex[] ap, int n, ref int info )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZPPCO factors a complex <double> hermitian positive definite matrix.
        //
        //  Discussion:
        //
        //    The routine also estimates the condition of the matrix.
        //
        //    The matrix is stored in packed form.
        //
        //    If RCOND is not needed, ZPPFA is slightly faster.
        //
        //    To solve A*X = B, follow ZPPCO by ZPPSL.
        //
        //    To compute inverse(A)*C, follow ZPPCO by ZPPSL.
        //
        //    To compute determinant(A), follow ZPPCO by ZPPDI.
        //
        //    To compute inverse(A), follow ZPPCO by ZPPDI.
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
        //    Input/output, complex <double> AP[N*(N+1)/2]; on input, the packed form of a
        //    hermitian matrix.  The columns of the upper triangle are stored
        //    sequentially in a one-dimensional array.  On output, an upper
        //    triangular matrix R, stored in packed form, so that A = hermitian(R) * R.
        //    If INFO != 0 , the factorization is not complete.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Output, double ZPPCO, an estimate of RCOND, the reciprocal condition of
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
        //  Local Parameters:
        //
        //    Local, complex <double> Z[N], a work vector whose contents are usually
        //    unimportant.  If A is close to a singular matrix, then Z is an
        //    approximate null vector in the sense that
        //      norm(A*Z) = RCOND * norm(A) * norm(Z).
        //
    {
        double anorm;
        Complex ek;
        int i;
        int ij;
        int j;
        int j1;
        int k;
        int kj;
        int kk;
        double rcond;
        double s;
        double sm;
        Complex t;
        Complex wk;
        Complex wkm;
        double ynorm;
        Complex[] z;
        //
        //  Find norm of A.
        //
        z = new Complex [n];

        j1 = 1;

        for (j = 1; j <= n; j++)
        {
            z[j - 1] = new Complex(BLAS1Z.dzasum(j, ap, 1, index: + j1 - 1), 0.0);
            ij = j1;
            j1 += j;

            for (i = 1; i <= j - 1; i++)
            {
                z[i - 1] = new Complex(z[i - 1].Real + typeMethods.zabs1(ap[ij - 1]), 0.0);
                ij += 1;
            }
        }

        anorm = 0.0;
        for (j = 0; j < n; j++)
        {
            anorm = Math.Max(anorm, z[j].Real);
        }

        //
        //  Factor.
        //
        info = ZPPFA.zppfa(ref ap, n);

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

        kk = 0;

        for (k = 1; k <= n; k++)
        {
            kk += k;

            if (typeMethods.zabs1(z[k - 1]) != 0.0)
            {
                ek = typeMethods.zsign1(ek, -z[k - 1]);
            }

            if (ap[kk - 1].Real < typeMethods.zabs1(ek - z[k - 1]))
            {
                s = ap[kk - 1].Real / typeMethods.zabs1(ek - z[k - 1]);
                BLAS1Z.zdscal(n, s, ref z, 1);
                ek = new Complex(s, 0.0) * ek;
            }

            wk = ek - z[k - 1];
            wkm = -ek - z[k - 1];
            s = typeMethods.zabs1(wk);
            sm = typeMethods.zabs1(wkm);
            wk /= ap[kk - 1];
            wkm /= ap[kk - 1];
            kj = kk + k;

            if (k + 1 <= n)
            {
                for (j = k + 1; j <= n; j++)
                {
                    sm += typeMethods.zabs1(z[j - 1] + wkm * Complex.Conjugate(ap[kj - 1]));
                    z[j - 1] += wk * Complex.Conjugate(ap[kj - 1]);
                    s += typeMethods.zabs1(z[j - 1]);
                    kj += j;
                }

                if (s < sm)
                {
                    t = wkm - wk;
                    wk = wkm;
                    kj = kk + k;
                    for (j = k + 1; j <= n; j++)
                    {
                        z[j - 1] += t * Complex.Conjugate(ap[kj - 1]);
                        kj += j;
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
            if (ap[kk - 1].Real < typeMethods.zabs1(z[k - 1]))
            {
                s = ap[kk - 1].Real / typeMethods.zabs1(z[k - 1]);
                BLAS1Z.zdscal(n, s, ref z, 1);
            }

            z[k - 1] /= ap[kk - 1];
            kk -= k;
            t = -z[k - 1];
            BLAS1Z.zaxpy(k - 1, t, ap, 1, ref z, 1, xIndex: + kk);
        }

        s = 1.0 / BLAS1Z.dzasum(n, z, 1);
        BLAS1Z.zdscal(n, s, ref z, 1);
        ynorm = 1.0;
        //
        //  Solve hermitian(R) * V = Y.
        //
        for (k = 1; k <= n; k++)
        {
            z[k - 1] -= BLAS1Z.zdotc(k - 1, ap, 1, z, 1, xIndex: + kk);
            kk += k;

            if (ap[kk - 1].Real < typeMethods.zabs1(z[k - 1]))
            {
                s = ap[kk - 1].Real / typeMethods.zabs1(z[k - 1]);
                BLAS1Z.zdscal(n, s, ref z, 1);
                ynorm = s * ynorm;
            }

            z[k - 1] /= ap[kk - 1];
        }

        s = 1.0 / BLAS1Z.dzasum(n, z, 1);
        BLAS1Z.zdscal(n, s, ref z, 1);
        ynorm = s * ynorm;
        //
        //  Solve R * Z = V.
        //
        for (k = n; 1 <= k; k--)
        {
            if (ap[kk - 1].Real < typeMethods.zabs1(z[k - 1]))
            {
                s = ap[kk - 1].Real / typeMethods.zabs1(z[k - 1]);
                BLAS1Z.zdscal(n, s, ref z, 1);
                ynorm = s * ynorm;
            }

            z[k - 1] /= ap[kk - 1];
            kk -= k;
            t = -z[k - 1];
            BLAS1Z.zaxpy(k - 1, t, ap, 1, ref z, 1, xIndex: + kk);
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