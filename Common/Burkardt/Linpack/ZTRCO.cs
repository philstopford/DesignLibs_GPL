using System;
using System.Numerics;
using Burkardt.BLAS;
using Burkardt.Types;

namespace Burkardt.Linpack;

public static class ZTRCO
{
    public static double ztrco(Complex[] t, int ldt, int n, int job)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZTRCO estimates the condition of a complex triangular matrix.
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
        //    Input, Complex T[LDT*N], the triangular matrix.  The zero
        //    elements of the matrix are not referenced, and the corresponding
        //    elements of the array can be used to store other information.
        //
        //    Input, int LDT, the leading dimension of T.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int JOB, indicates if matrix is upper or lower triangular.
        //    0, lower triangular.
        //    nonzero, upper triangular.
        //
        //    Output, double ZTRCO, an estimate of RCOND, the reciprocal condition of T.
        //    For the system T*X = B, relative perturbations in T and B of size
        //    EPSILON may cause relative perturbations in X of size (EPSILON/RCOND).
        //    If RCOND is so small that the logical expression
        //      1.0 + RCOND == 1.0
        //    is true, then T may be singular to working precision.  In particular,
        //    RCOND is zero if exact singularity is detected or the estimate
        //    underflows.
        //
        //  Local Parameters:
        //
        //    Workspace, Complex Z[N], a work vector whose contents are usually
        //    unimportant.  If T is close to a singular matrix, then Z is
        //    an approximate null vector in the sense that
        //      norm(A*Z) = RCOND * norm(A) * norm(Z).
        //
    {
        int i;
        int i1;
        int j;
        int k;
        int kk;
        double rcond;
        double s;
        Complex w;

        bool lower = job == 0;
        //
        //  Compute 1-norm of T
        //
        double tnorm = 0.0;

        for (j = 1; j <= n; j++)
        {
            int l;
            switch (lower)
            {
                case true:
                    l = n + 1 - j;
                    i1 = j;
                    break;
                default:
                    l = j;
                    i1 = 1;
                    break;
            }

            tnorm = Math.Max(tnorm,
                BLAS1Z.dzasum(l, t, 1, index: +i1 - 1 + (j - 1) * ldt));
        }

        //
        //  RCOND = 1/(norm(T)*(estimate of norm(inverse(T)))).
        //
        //  Estimate = norm(Z)/norm(Y) where T*Z = Y and hermitian(T)*Y = E.
        //
        //  Hermitian(T) is the Complex.Conjugateugate transpose of T.
        //
        //  The components of E are chosen to cause maximum local
        //  growth in the elements of Y.
        //
        //  The vectors are frequently rescaled to avoid overflow.
        //
        //  Solve hermitian(T)*Y = E.
        //
        Complex ek = new Complex(1.0, 0.0);
        Complex[] z = new Complex[n];
        for (i = 0; i < n; i++)
        {
            z[i] = new Complex(0.0, 0.0);
        }

        for (kk = 1; kk <= n; kk++)
        {
            k = lower switch
            {
                true => n + 1 - kk,
                _ => kk
            };

            if (typeMethods.zabs1(z[k - 1]) != 0.0)
            {
                ek = typeMethods.zsign1(ek, -z[k - 1]);
            }

            if (typeMethods.zabs1(t[k - 1 + (k - 1) * ldt]) < typeMethods.zabs1(ek - z[k - 1]))
            {
                s = typeMethods.zabs1(t[k - 1 + (k - 1) * ldt]) / typeMethods.zabs1(ek - z[k - 1]);
                BLAS1Z.zdscal(n, s, ref z, 1);
                ek = new Complex(s, 0.0) * ek;
            }

            Complex wk = ek - z[k - 1];
            Complex wkm = -ek - z[k - 1];
            s = typeMethods.zabs1(wk);
            double sm = typeMethods.zabs1(wkm);

            if (typeMethods.zabs1(t[k - 1 + (k - 1) * ldt]) != 0.0)
            {
                wk /= Complex.Conjugate(t[k - 1 + (k - 1) * ldt]);
                wkm /= Complex.Conjugate(t[k - 1 + (k - 1) * ldt]);
            }
            else
            {
                wk = new Complex(1.0, 0.0);
                wkm = new Complex(1.0, 0.0);
            }

            if (kk != n)
            {
                int j2;
                int j1;
                switch (lower)
                {
                    case true:
                        j1 = 1;
                        j2 = k - 1;
                        break;
                    default:
                        j1 = k + 1;
                        j2 = n;
                        break;
                }

                for (j = j1; j <= j2; j++)
                {
                    sm += typeMethods.zabs1(z[j - 1] + wkm * Complex.Conjugate(t[k - 1 + (j - 1) * ldt]));
                    z[j - 1] += wk * Complex.Conjugate(t[k - 1 + (j - 1) * ldt]);
                    s += typeMethods.zabs1(z[j - 1]);
                }

                if (s < sm)
                {
                    w = wkm - wk;
                    wk = wkm;
                    for (j = j1; j <= j2; j++)
                    {
                        z[j - 1] += w * Complex.Conjugate(t[k - 1 + (j - 1) * ldt]);
                    }
                }
            }

            z[k - 1] = wk;
        }

        s = 1.0 / BLAS1Z.dzasum(n, z, 1);
        BLAS1Z.zdscal(n, s, ref z, 1);
        double ynorm = 1.0;
        //
        //  Solve T*Z = Y.
        //
        for (kk = 1; kk <= n; kk++)
        {
            k = lower switch
            {
                true => kk,
                _ => n + 1 - kk
            };

            if (typeMethods.zabs1(t[k - 1 + (k - 1) * ldt]) < typeMethods.zabs1(z[k - 1]))
            {
                s = typeMethods.zabs1(t[k - 1 + (k - 1) * ldt]) / typeMethods.zabs1(z[k - 1]);
                BLAS1Z.zdscal(n, s, ref z, 1);
                ynorm = s * ynorm;
            }

            if (typeMethods.zabs1(t[k - 1 + (k - 1) * ldt]) != 0.0)
            {
                z[k - 1] /= t[k - 1 + (k - 1) * ldt];
            }
            else
            {
                z[k - 1] = new Complex(1.0, 0.0);
            }

            i1 = lower switch
            {
                true => k + 1,
                _ => 1
            };

            if (kk < n)
            {
                w = -z[k - 1];
                BLAS1Z.zaxpy(n - kk, w, t, 1, ref z, 1, xIndex: +i1 - 1 + (k - 1) * ldt, yIndex: +i1 - 1);
            }
        }

        //
        //  Make ZNORM = 1.
        //
        s = 1.0 / BLAS1Z.dzasum(n, z, 1);
        BLAS1Z.zdscal(n, s, ref z, 1);
        ynorm = s * ynorm;

        if (tnorm != 0.0)
        {
            rcond = ynorm / tnorm;
        }
        else
        {
            rcond = 0.0;
        }

        return rcond;
    }

}