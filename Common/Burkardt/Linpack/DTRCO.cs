using System;
using Burkardt.BLAS;
using Burkardt.Types;

namespace Burkardt.Linpack;

public static class DTRCO
{
    public static double dtrco(double[] t, int ldt, int n, ref double[] z, int job)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DTRCO estimates the condition of a real triangular matrix.
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
        //    Input, double T[LDT*N], the triangular matrix.  The zero
        //    elements of the matrix are not referenced, and the corresponding
        //    elements of the array can be used to store other information.
        //
        //    Input, int LDT, the leading dimension of the array T.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Output, double Z[N] a work vector whose contents are usually
        //    unimportant.  If T is close to a singular matrix, then Z is an
        //    approximate null vector in the sense that
        //      norm(A*Z) = RCOND * norm(A) * norm(Z).
        //
        //    Input, int JOB, indicates the shape of T:
        //    0, T is lower triangular.
        //    nonzero, T is upper triangular.
        //
        //    Output, double DTRCO, an estimate of the reciprocal condition RCOND
        //    of T.  For the system T*X = B, relative perturbations in T and B of size
        //    EPSILON may cause relative perturbations in X of size EPSILON/RCOND.
        //    If RCOND is so small that the logical expression
        //      1.0D+00 + RCOND == 1.0D+00
        //    is true, then T may be singular to working precision.  In particular,
        //    RCOND is zero if exact singularity is detected or the estimate underflows.
        //
    {
        int i;
        int i1;
        int j;
        int k;
        int kk;
        double rcond;
        double s;
        double w;

        bool lower = job == 0;
        //
        //  Compute the 1-norm of T.
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

            tnorm = Math.Max(tnorm, BLAS1D.dasum(l, t, 1, index: +i1 - 1 + (j - 1) * ldt));
        }

        //
        //  RCOND = 1/(norm(T)*(estimate of norm(inverse(T)))).
        //
        //  Estimate = norm(Z)/norm(Y) where T * Z = Y and T' * Y = E.
        //
        //  T' is the transpose of T.
        //
        //  The components of E are chosen to cause maximum local
        //  growth in the elements of Y.
        //
        //  The vectors are frequently rescaled to avoid overflow.
        //
        //  Solve T' * Y = E.
        //
        double ek = 1.0;
        for (i = 1; i <= n; i++)
        {
            z[i - 1] = 0.0;
        }

        for (kk = 1; kk <= n; kk++)
        {
            k = lower switch
            {
                true => n + 1 - kk,
                _ => kk
            };

            if (z[k - 1] != 0.0)
            {
                ek = typeMethods.r8_sign(-z[k - 1]) * ek;
            }

            if (Math.Abs(t[k - 1 + (k - 1) * ldt]) < Math.Abs(ek - z[k - 1]))
            {
                s = Math.Abs(t[k - 1 + (k - 1) * ldt]) / Math.Abs(ek - z[k - 1]);
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

            if (t[k - 1 + (k - 1) * ldt] != 0.0)
            {
                wk /= t[k - 1 + (k - 1) * ldt];
                wkm /= t[k - 1 + (k - 1) * ldt];
            }
            else
            {
                wk = 1.0;
                wkm = 1.0;
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
                    sm += Math.Abs(z[j - 1] + wkm * t[k - 1 + (j - 1) * ldt]);
                    z[j - 1] += wk * t[k - 1 + (j - 1) * ldt];
                    s += Math.Abs(z[j - 1]);
                }

                if (s < sm)
                {
                    w = wkm - wk;
                    wk = wkm;
                    for (j = j1; j <= j2; j++)
                    {
                        z[j - 1] += w * t[k - 1 + (j - 1) * ldt];
                    }
                }
            }

            z[k - 1] = wk;
        }

        double temp = BLAS1D.dasum(n, z, 1);

        for (i = 1; i <= n; i++)
        {
            z[i - 1] /= temp;
        }

        double ynorm = 1.0;
        //
        //  Solve T * Z = Y.
        //
        for (kk = 1; kk <= n; kk++)
        {
            k = lower switch
            {
                true => kk,
                _ => n + 1 - kk
            };

            if (Math.Abs(t[k - 1 + (k - 1) * ldt]) < Math.Abs(z[k - 1]))
            {
                s = Math.Abs(t[k - 1 + (k - 1) * ldt]) / Math.Abs(z[k - 1]);
                for (i = 1; i <= n; i++)
                {
                    z[i - 1] = s * z[i - 1];
                }

                ynorm = s * ynorm;
            }

            if (t[k - 1 + (k - 1) * ldt] != 0.0)
            {
                z[k - 1] /= t[k - 1 + (k - 1) * ldt];
            }
            else
            {
                z[k - 1] = 1.0;
            }

            i1 = lower switch
            {
                true => k + 1,
                _ => 1
            };

            if (kk < n)
            {
                w = -z[k - 1];
                BLAS1D.daxpy(n - kk, w, t, 1, ref z, 1, xIndex: +i1 - 1 + (k - 1) * ldt, yIndex: +i1 - 1);
            }
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