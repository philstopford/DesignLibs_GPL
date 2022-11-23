using System;
using System.Numerics;
using Burkardt.BLAS;
using Burkardt.Types;

namespace Burkardt.Linpack;

public static class ZSVDC
{
    public static int zsvdc(ref Complex[] x, int ldx, int n, int p,
            ref Complex[] s, ref Complex[] e, ref Complex[] u, int ldu,
            ref Complex[] v, int ldv, int job)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZSVDC applies the singular value decompostion to an N by P matrix.
        //
        //  Discussion:
        //
        //    The routine reduces a Complex N by P matrix X, by unitary
        //    transformations U and V, to diagonal form.
        //
        //    The diagonal elements, S(I), are the singular values of Z.  The
        //    columns of U are the corresponding left singular vectors,
        //    and the columns of V the right singular vectors.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 May 2007
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
        //    Input/output, Complex X[LDX*P]; on input, the matrix whose singular
        //    value decomposition is to be computed.  X is destroyed on output.
        //
        //    Input, int LDX, the leading dimension of X.  N <= LDX.
        //
        //    Input, int N, the number of rows of the matrix.
        //
        //    Input, int P, the number of columns of the matrix X.
        //
        //    Output, Complex S[MM], where MM = min ( N + 1, P ), the first min ( N, P )
        //    entries of S contain the singular values of X arranged in descending
        //    order of magnitude.
        //
        //    Output, Complex E[MM], where MM = min ( N + 1, P ),
        //    ordinarily contains zeros on output.  However, see the discussion
        //    of INFO for exceptions.
        //
        //    Output, Complex U[LDU*K].  If JOBA == 1 then K == n; if JOBA >= 2,
        //    then K == min ( N, P ).  U contains the matrix of left singular vectors.
        //    U is not referenced if JOBA == 0.  If N <= P or if JOBA > 2,
        //    then U may be identified with X in the subroutine call.
        //
        //    Input, int LDU, the leading dimension of U.  N <= LDU.
        //
        //    Output, Complex V[LDV*P], if requested, the matrix of right singular
        //    vectors.  If P <= N, V may be identified with X in the subroutine call.
        //
        //    Input, int LDV, the leading dimension of V.  P <= LDV.
        //
        //    Input, int JOB, controls the computation of the singular vectors.
        //    It has the decimal expansion AB meaning:
        //    A =  0, do not compute the left singular vectors.
        //    A =  1, return the N left singular vectors in U.
        //    A >= 2, returns the first min ( N, P ) left singular vectors in U.
        //    B =  0, do not compute the right singular vectors.
        //    B =  1, return the right singular vectors in V.
        //
        //    Output, int ZSVDC, the value of INFO.  The singular values and their
        //    corresponding singular vectors are correct for entries,
        //    S(INFO+1), S(INFO+2), ..., S(M).  Here M = min ( N, P ).  Thus if
        //    INFO == 0, all the singular values and their vectors are correct.
        //    In any event, the matrix
        //      B = hermitian(U)*X*V
        //    is the bidiagonal matrix with the elements of S on its diagonal
        //    and the elements of E on its super-diagonal.  Hermitian(U)
        //    is the Complex.Conjugateugate-transpose of U.  Thus the singular values of X
        //    and B are the same.
        //
    {
        double cs = 0;
        int i;
        int j;
        int l;
        int ll;
        int lp1;
        int ls = 0;
        const int maxit = 30;
        double sn = 0;
        Complex t;

        Complex[] work = new Complex[n];
        //
        //  Determine what is to be computed.
        //
        bool wantu = false;
        bool wantv = false;
        int jobu = job % 100 / 10;

        int ncu = jobu switch
        {
            > 1 => Math.Min(n, p),
            _ => n
        };

        if (jobu != 0)
        {
            wantu = true;
        }

        if (job % 10 != 0)
        {
            wantv = true;
        }

        //
        //  Reduce X to bidiagonal form, storing the diagonal elements
        //  in S and the super-diagonal elements in E.
        //
        int info = 0;
        int nct = Math.Min(n - 1, p);
        int nrt = Math.Max(0, Math.Min(p - 2, n));
        int lu = Math.Max(nct, nrt);

        for (l = 1; l <= lu; l++)
        {
            lp1 = l + 1;
            //
            //  Compute the transformation for the L-th column and
            //  place the L-th diagonal in S(L).
            //
            if (l <= nct)
            {
                s[l - 1] = new Complex(BLAS1Z.dznrm2(n - l + 1, x, 1, index: +l - 1 + (l - 1) * ldx), 0.0);

                if (typeMethods.zabs1(s[l - 1]) != 0.0)
                {
                    if (typeMethods.zabs1(x[l - 1 + (l - 1) * ldx]) != 0.0)
                    {
                        s[l - 1] = typeMethods.zsign2(s[l - 1], x[l - 1 + (l - 1) * ldx]);
                    }

                    t = new Complex(1.0, 0.0) / s[l - 1];
                    BLAS1Z.zscal(n - l + 1, t, ref x, 1, index: +l - 1 + (l - 1) * ldx);
                    x[l - 1 + (l - 1) * ldx] = new Complex(1.0, 0.0) + x[l - 1 + (l - 1) * ldx];
                }

                s[l - 1] = -s[l - 1];
            }

            for (j = lp1; j <= p; j++)
            {
                if (l <= nct)
                {
                    if (typeMethods.zabs1(s[l - 1]) != 0.0)
                    {
                        t = -BLAS1Z.zdotc(n - l + 1, x, 1, x, 1, xIndex: +l - 1 + (l - 1) * ldx,
                                yIndex: +l - 1 + (j - 1) * ldx)
                            / x[l - 1 + (l - 1) * ldx];
                        BLAS1Z.zaxpy(n - l + 1, t, x, 1, ref x, 1, xIndex: +l - 1 + (l - 1) * ldx,
                            yIndex: +l - 1 + (j - 1) * ldx);
                    }
                }

                //
                //  Place the L-th row of X into E for the
                //  subsequent calculation of the row transformation.
                //
                e[j - 1] = Complex.Conjugate(x[l - 1 + (j - 1) * ldx]);
            }

            switch (wantu)
            {
                //
                //  Place the transformation in U for subsequent back multiplication.
                //
                case true when l <= nct:
                {
                    for (i = l; i <= n; i++)
                    {
                        u[i - 1 + (l - 1) * ldu] = x[i - 1 + (l - 1) * ldx];
                    }

                    break;
                }
            }

            if (l <= nrt)
            {
                //
                //  Compute the L-th row transformation and place the
                //  L-th super-diagonal in E(L).
                //
                e[l - 1] = new Complex(BLAS1Z.dznrm2(p - l, e, 1, index: +lp1 - 1), 0.0);

                if (typeMethods.zabs1(e[l - 1]) != 0.0)
                {
                    if (typeMethods.zabs1(e[lp1 - 1]) != 0.0)
                    {
                        e[l - 1] = typeMethods.zsign2(e[l - 1], e[lp1 - 1]);
                    }

                    t = new Complex(1.0, 0.0) / e[l - 1];
                    BLAS1Z.zscal(p - l, t, ref e, 1, index: +lp1 - 1);
                    e[lp1 - 1] = new Complex(1.0, 0.0) + e[lp1 - 1];
                }

                e[l - 1] = -Complex.Conjugate(e[l - 1]);
                //
                //  Apply the transformation.
                //
                if (lp1 <= n && typeMethods.zabs1(e[l - 1]) != 0.0)
                {
                    for (j = lp1; j <= n; j++)
                    {
                        work[j - 1] = new Complex(0.0, 0.0);
                    }

                    for (j = lp1; j <= p; j++)
                    {
                        BLAS1Z.zaxpy(n - l, e[j - 1], x, 1, ref work, 1, xIndex: +lp1 - 1 + (j - 1) * ldx,
                            yIndex: +lp1 - 1);
                    }

                    for (j = lp1; j <= p; j++)
                    {
                        BLAS1Z.zaxpy(n - l, Complex.Conjugate(-e[j - 1] / e[lp1 - 1]), work,
                            1, ref x, 1, xIndex: +lp1 - 1, yIndex: +lp1 - 1 + (j - 1) * ldx);
                    }
                }

                switch (wantv)
                {
                    //
                    //  Place the transformation in V for subsequent back multiplication.
                    //
                    case true:
                    {
                        for (i = lp1; i <= p; i++)
                        {
                            v[i - 1 + (l - 1) * ldv] = e[i - 1];
                        }

                        break;
                    }
                }
            }
        }

        //
        //  Set up the final bidiagonal matrix of order M.
        //
        int m = Math.Min(p, n + 1);
        int nctp1 = nct + 1;
        int nrtp1 = nrt + 1;

        if (nct < p)
        {
            s[nctp1 - 1] = x[nctp1 - 1 + (nctp1 - 1) * ldx];
        }

        if (n < m)
        {
            s[m - 1] = new Complex(0.0, 0.0);
        }

        if (nrtp1 < m)
        {
            e[nrtp1 - 1] = x[nrtp1 - 1 + (m - 1) * ldx];
        }

        e[m - 1] = new Complex(0.0, 0.0);
        switch (wantu)
        {
            //
            //  If required, generate U.
            //
            case true:
            {
                for (j = nctp1; j <= ncu; j++)
                {
                    for (i = 1; i <= n; i++)
                    {
                        u[i - 1 + (j - 1) * ldu] = new Complex(0.0, 0.0);
                    }

                    u[j - 1 + (j - 1) * ldu] = new Complex(1.0, 0.0);
                }

                for (ll = 1; ll <= nct; ll++)
                {
                    l = nct - ll + 1;

                    if (typeMethods.zabs1(s[l - 1]) != 0.0)
                    {
                        lp1 = l + 1;

                        for (j = l + 1; j <= ncu; j++)
                        {
                            t = -BLAS1Z.zdotc(n - l + 1, u, 1, u, 1, xIndex: +l - 1 + (l - 1) * ldu,
                                    yIndex: +l - 1 + (j - 1) * ldu)
                                / u[l - 1 + (l - 1) * ldu];
                            BLAS1Z.zaxpy(n - l + 1, t, u, 1, ref u, 1, xIndex: +l - 1 + (l - 1) * ldu,
                                yIndex: +l - 1 + (j - 1) * ldu);
                        }

                        BLAS1Z.zscal(n - l + 1, new Complex(-1.0, 0.0), ref u, 1, index: +l - 1 + (l - 1) * ldu);
                        u[l - 1 + (l - 1) * ldu] = new Complex(1.0, 0.0) + u[l - 1 + (l - 1) * ldu];
                        for (i = 1; i <= l - 1; i++)
                        {
                            u[i - 1 + (l - 1) * ldu] = new Complex(0.0, 0.0);
                        }
                    }
                    else
                    {
                        for (i = 1; i <= n; i++)
                        {
                            u[i - 1 + (l - 1) * ldu] = new Complex(0.0, 0.0);
                        }

                        u[l - 1 + (l - 1) * ldu] = new Complex(1.0, 0.0);
                    }
                }

                break;
            }
        }

        switch (wantv)
        {
            //
            //  If it is required, generate V.
            //
            case true:
            {
                for (ll = 1; ll <= p; ll++)
                {
                    l = p - ll + 1;
                    lp1 = l + 1;

                    if (l <= nrt)
                    {
                        if (typeMethods.zabs1(e[l - 1]) != 0.0)
                        {
                            for (j = lp1; j <= p; j++)
                            {
                                t = -BLAS1Z.zdotc(p - l, v, 1, v, 1, xIndex: +lp1 - 1 + (l - 1) * ldv,
                                        yIndex: +lp1 - 1 + (j - 1) * ldv)
                                    / v[lp1 - 1 + (l - 1) * ldv];
                                BLAS1Z.zaxpy(p - l, t, v, 1, ref v, 1, xIndex: +lp1 - 1 + (l - 1) * ldv,
                                    yIndex: +lp1 - 1 + (j - 1) * ldv);
                            }
                        }
                    }

                    for (i = 1; i <= p; i++)
                    {
                        v[i - 1 + (l - 1) * ldv] = new Complex(0.0, 0.0);
                    }

                    v[l - 1 + (l - 1) * ldv] = new Complex(1.0, 0.0);
                }

                break;
            }
        }

        //
        //  Transform S and E so that they are real.
        //
        for (i = 1; i <= m; i++)
        {
            Complex 
                r;
            if (typeMethods.zabs1(s[i - 1]) != 0.0)
            {
                t = new Complex(Complex.Abs(s[i - 1]), 0.0);
                r = s[i - 1] / t;
                s[i - 1] = t;

                if (i < m)
                {
                    e[i - 1] /= r;
                }

                switch (wantu)
                {
                    case true:
                        BLAS1Z.zscal(n, r, ref u, 1, index: +0 + (i - 1) * ldu);
                        break;
                }
            }

            if (i == m)
            {
                break;
            }

            if (typeMethods.zabs1(e[i - 1]) != 0.0)
            {
                t = new Complex(Complex.Abs(e[i - 1]), 0.0);
                r = t / e[i - 1];
                e[i - 1] = t;
                s[i] *= r;

                switch (wantv)
                {
                    case true:
                        BLAS1Z.zscal(p, r, ref v, 1, index: +0 + i * ldv);
                        break;
                }
            }
        }

        //
        //  Main iteration loop for the singular values.
        //
        int mm = m;
        int iter = 0;

        for (;;)
        {
            //
            //  Quit if all the singular values have been found.
            //
            if (m == 0)
            {
                break;
            }

            //
            //  If too many iterations have been performed, set flag and return.
            //
            if (maxit <= iter)
            {
                info = m;
                break;
            }

            //
            //  This section of the program inspects for negligible elements in S and E.
            //
            //  On completion, the variables KASE and L are set as follows.
            //
            //  KASE = 1     if S(M) and E(L-1) are negligible and L < M
            //  KASE = 2     if S(L) is negligible and L < M
            //  KASE = 3     if E(L-1) is negligible, L < M, and
            //               S(L), ..., S(M) are not negligible (QR step).
            //  KASE = 4     if E(M-1) is negligible (convergence).
            //
            double test;
            double ztest;
            for (ll = 1; ll <= m; ll++)
            {
                l = m - ll;

                if (l == 0)
                {
                    break;
                }

                test = Complex.Abs(s[l - 1]) + Complex.Abs(s[l]);
                ztest = test + Complex.Abs(e[l - 1]);

                if (!(Math.Abs(ztest - test) <= typeMethods.r8_epsilon()))
                {
                    continue;
                }

                e[l - 1] = new Complex(0.0, 0.0);
                break;
            }

            int kase;
            if (l == m - 1)
            {
                kase = 4;
            }
            else
            {
                lp1 = l + 1;
                int mp1 = m + 1;

                int lls;
                for (lls = lp1; lls <= mp1; lls++)
                {
                    ls = m - lls + lp1;

                    if (ls == l)
                    {
                        break;
                    }

                    test = 0.0;

                    if (ls != m)
                    {
                        test += Complex.Abs(e[ls - 1]);
                    }

                    if (ls != l + 1)
                    {
                        test += Complex.Abs(e[ls - 2]);
                    }

                    ztest = test + Complex.Abs(s[ls - 1]);

                    if (!(Math.Abs(ztest - test) <= typeMethods.r8_epsilon()))
                    {
                        continue;
                    }

                    s[ls - 1] = new Complex(0.0, 0.0);
                    break;
                }

                if (ls == l)
                {
                    kase = 3;
                }
                else if (ls == m)
                {
                    kase = 1;
                }
                else
                {
                    kase = 2;
                    l = ls;
                }
            }

            l += 1;
            int mm1;
            double t1;
            double f;
            int k;
            switch (kase)
            {
                //
                //  Deflate negligible S(M).
                //
                case 1:
                {
                    mm1 = m - 1;
                    f = e[m - 2].Real;
                    e[m - 2] = new Complex(0.0, 0.0);

                    int kk;
                    for (kk = 1; kk <= mm1; kk++)
                    {
                        k = mm1 - kk + l;
                        t1 = s[k - 1].Real;
                        BLAS1D.drotg(ref t1, ref f, ref cs, ref sn);
                        s[k - 1] = new Complex(t1, 0.0);

                        if (k != l)
                        {
                            f = -sn * e[k - 2].Real;
                            e[k - 2] = cs * e[k - 2];
                        }

                        switch (wantv)
                        {
                            case true:
                                BLAS1Z.zdrot(p, ref v, 1, ref v, 1, cs, sn, xIndex: +0 + (k - 1) * ldv,
                                    yIndex: +0 + (m - 1) * ldv);
                                break;
                        }
                    }

                    break;
                }
                //
                //  Split at negligible S(L).
                //
                case 2:
                {
                    f = e[l - 2].Real;
                    e[l - 2] = new Complex(0.0, 0.0);

                    for (k = l; k <= m; k++)
                    {
                        t1 = s[k - 1].Real;
                        BLAS1D.drotg(ref t1, ref f, ref cs, ref sn);
                        s[k - 1] = new Complex(t1, 0.0);
                        f = -sn * e[k - 1].Real;
                        e[k - 1] = cs * e[k - 1];

                        switch (wantu)
                        {
                            case true:
                                BLAS1Z.zdrot(n, ref u, 1, ref u, 1, cs, sn, xIndex: +0 + (k - 1) * ldu,
                                    yIndex: +0 + (l - 2) * ldu);
                                break;
                        }
                    }

                    break;
                }
                //
                //  Perform one QR step.
                //
                case 3:
                {
                    //
                    //  Calculate the shift.
                    //
                    double scale = Math.Max(Complex.Abs(s[m - 1]),
                        Math.Max(Complex.Abs(s[m - 2]),
                            Math.Max(Complex.Abs(e[m - 2]),
                                Math.Max(Complex.Abs(s[l - 1]), Complex.Abs(e[l - 1])))));

                    double sm = s[m - 1].Real / scale;
                    double 
                        smm1 = s[m - 2].Real / scale;
                    double 
                        emm1 = e[m - 2].Real / scale;
                    double sl = s[l - 1].Real / scale;
                    double el = e[l - 1].Real / scale;
                    double b = ((smm1 + sm) * (smm1 - sm) + emm1 * emm1) / 2.0;
                    double c = sm * emm1 * (sm * emm1);
                    double shift = 0.0;

                    if (b != 0.0 || c != 0.0)
                    {
                        shift = b switch
                        {
                            < 0.0 => -shift,
                            _ => Math.Sqrt(b * b + c)
                        };

                        shift = c / (b + shift);
                    }

                    f = (sl + sm) * (sl - sm) + shift;
                    double g = sl * el;
                    //
                    //  Chase zeros.
                    //
                    mm1 = m - 1;

                    for (k = l; k <= mm1; k++)
                    {
                        BLAS1D.drotg(ref f, ref g, ref cs, ref sn);

                        if (k != l)
                        {
                            e[k - 2] = new Complex(f, 0.0);
                        }

                        f = cs * s[k - 1].Real + sn * e[k - 1].Real;
                        e[k - 1] = cs * e[k - 1] - sn * s[k - 1];
                        g = sn * s[k].Real;
                        s[k] = cs * s[k];

                        switch (wantv)
                        {
                            case true:
                                BLAS1Z.zdrot(p, ref v, 1, ref v, 1, cs, sn, xIndex: +0 + (k - 1) * ldv,
                                    yIndex: +0 + k * ldv);
                                break;
                        }

                        BLAS1D.drotg(ref f, ref g, ref cs, ref sn);
                        s[k - 1] = new Complex(f, 0.0);
                        f = cs * e[k - 1].Real + sn * s[k].Real;
                        s[k] = -sn * e[k - 1] + cs * s[k];
                        g = sn * e[k].Real;
                        e[k] = cs * e[k];

                        switch (wantu)
                        {
                            case true when k < n:
                                BLAS1Z.zdrot(n, ref u, 1, ref u, 1, cs, sn, xIndex: +0 + (k - 1) * ldu,
                                    yIndex: +0 + k * ldu);
                                break;
                        }
                    }

                    e[m - 2] = new Complex(f, 0.0);
                    iter += 1;
                    break;
                }
                //
                //  Convergence.
                //
                case 4:
                {
                    switch (s[l - 1].Real)
                    {
                        //
                        //  Make the singular value positive.
                        //
                        case < 0.0:
                        {
                            s[l - 1] = -s[l - 1];
                            switch (wantv)
                            {
                                case true:
                                    BLAS1Z.zscal(p, new Complex(-1.0, 0.0), ref v, 1, index: +0 + (l - 1) * ldv);
                                    break;
                            }

                            break;
                        }
                    }

                    //
                    //  Order the singular values.
                    //
                    while (l != mm)
                    {
                        if (s[l].Real <= s[l - 1].Real)
                        {
                            break;
                        }

                        t = s[l - 1];
                        s[l - 1] = s[l];
                        s[l] = t;

                        switch (wantv)
                        {
                            case true when l < p:
                                BLAS1Z.zswap(p, ref v, 1, ref v, 1, xIndex: +0 + (l - 1) * ldv, yIndex: +0 + l * ldv);
                                break;
                        }

                        switch (wantu)
                        {
                            case true when l < n:
                                BLAS1Z.zswap(n, ref u, 1, ref u, 1, xIndex: +0 + (l - 1) * ldu, yIndex: +0 + l * ldu);
                                break;
                        }

                        l += 1;
                    }

                    iter = 0;
                    m -= 1;
                    break;
                }
            }
        }

        return info;
    }

}