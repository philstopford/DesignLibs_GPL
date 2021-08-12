using System;
using Burkardt.AppliedStatistics;
using Burkardt.BLAS;
using Burkardt.Types;

namespace Burkardt.MatrixNS
{
    public static class DSVDC
    {
        public static int dsvdc(ref double[] a, int lda, int m, int n, ref double[] s, ref double[] e,
        ref double[] u, int ldu, ref double[] v, int ldv, double[] work, int job )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DSVDC computes the singular value decomposition of a real rectangular matrix.
        //
        //  Discussion:
        //
        //    This routine reduces an M by N matrix A to diagonal form by orthogonal
        //    transformations U and V.  The diagonal elements S(I) are the singular
        //    values of A.  The columns of U are the corresponding left singular
        //    vectors, and the columns of V the right singular vectors.
        //
        //    The form of the singular value decomposition is then
        //
        //      A(MxN) = U(MxM) * S(MxN) * V(NxN)'
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
        //    Input/output, double A[LDA*N].  On input, the M by N matrix whose
        //    singular value decomposition is to be computed.  On output, the matrix
        //    has been destroyed.  Depending on the user's requests, the matrix may 
        //    contain other useful information.
        //
        //    Input, int LDA, the leading dimension of the array A.
        //    LDA must be at least M.
        //
        //    Input, int M, the number of rows of the matrix.
        //
        //    Input, int N, the number of columns of the matrix A.
        //
        //    Output, double S[MM], where MM = min(M+1,N).  The first
        //    min(M,N) entries of S contain the singular values of A arranged in
        //    descending order of magnitude.
        //
        //    Output, double E[MM], where MM = min(M+1,N), ordinarily contains zeros.
        //    However see the discussion of INFO for exceptions.
        //
        //    Output, double U[LDU*K].  If JOBA = 1 then K = M;
        //    if 2 <= JOBA, then K = min(M,N).  U contains the M by M matrix of left singular
        //    vectors.  U is not referenced if JOBA = 0.  If M <= N or if JOBA = 2, then
        //    U may be identified with A in the subroutine call.
        //
        //    Input, int LDU, the leading dimension of the array U.
        //    LDU must be at least M.
        //
        //    Output, double V[LDV*N], the N by N matrix of right singular vectors.
        //    V is not referenced if JOB is 0.  If N <= M, then V may be identified
        //    with A in the subroutine call.
        //
        //    Input, int LDV, the leading dimension of the array V.
        //    LDV must be at least N.
        //
        //    Workspace, double WORK[M].
        //
        //    Input, int JOB, controls the computation of the singular
        //    vectors.  It has the decimal expansion AB with the following meaning:
        //      A =  0, do not compute the left singular vectors.
        //      A =  1, return the M left singular vectors in U.
        //      A >= 2, return the first min(M,N) singular vectors in U.
        //      B =  0, do not compute the right singular vectors.
        //      B =  1, return the right singular vectors in V.
        //
        //    Output, int *DSVDC, status indicator INFO.
        //    The singular values (and their corresponding singular vectors)
        //    S(*INFO+1), S(*INFO+2),...,S(MN) are correct.  Here MN = min ( M, N ).
        //    Thus if *INFO is 0, all the singular values and their vectors are
        //    correct.  In any event, the matrix B = U' * A * V is the bidiagonal
        //    matrix with the elements of S on its diagonal and the elements of E on
        //    its superdiagonal.  Thus the singular values of A and B are the same.
        //
        {
            double b = 0;
            double c = 0;
            double cs = 0;
            double el = 0;
            double emm1 = 0;
            double f = 0;
            double g = 0;
            int i = 0;
            int info = 0;
            int iter = 0;
            int j = 0;
            int jobu = 0;
            int k = 0;
            int kase = 0;
            int kk = 0;
            int l = 0;
            int ll = 0;
            int lls = 0;
            int ls = 0;
            int lu = 0;
            int maxit = 30;
            int mm = 0;
            int mm1 = 0;
            int mn = 0;
            int nct = 0;
            int nctp1 = 0;
            int ncu = 0;
            int nrt = 0;
            int nrtp1 = 0;
            double scale = 0;
            double shift = 0;
            double sl = 0;
            double sm = 0;
            double smm1 = 0;
            double sn = 0;
            double t = 0;
            double t1 = 0;
            double test = 0;
            bool wantu = false;
            bool wantv = false;
            double ztest = 0;
            //
            //  Determine what is to be computed.
            //
            info = 0;
            wantu = false;
            wantv = false;
            jobu = (job % 100) / 10;

            if (1 < jobu)
            {
                ncu = Math.Min(m, n);
            }
            else
            {
                ncu = m;
            }

            if (jobu != 0)
            {
                wantu = true;
            }

            if ((job % 10) != 0)
            {
                wantv = true;
            }

            //
            //  Reduce A to bidiagonal form, storing the diagonal elements
            //  in S and the super-diagonal elements in E.
            //
            nct = Math.Min(m - 1, n);
            nrt = Math.Max(0, Math.Min(m, n - 2));
            lu = Math.Max(nct, nrt);

            for (l = 1; l <= lu; l++)
            {
                //
                //  Compute the transformation for the L-th column and
                //  place the L-th diagonal in S(L).
                //
                if (l <= nct)
                {
                    s[l - 1] = BLAS1D.dnrm2(m - l + 1, a, 1,  + l - 1 + (l - 1) * lda);

                    if (s[l - 1] != 0.0)
                    {
                        if (a[l - 1 + (l - 1) * lda] != 0.0)
                        {
                            s[l - 1] = typeMethods.r8_sign(a[l - 1 + (l - 1) * lda]) * Math.Abs(s[l - 1]);
                        }

                        BLAS1D.dscal(m - l + 1, 1.0 / s[l - 1], ref a, 1, + l - 1 + (l - 1) * lda);
                        a[l - 1 + (l - 1) * lda] = 1.0 + a[l - 1 + (l - 1) * lda];
                    }

                    s[l - 1] = -s[l - 1];
                }

                for (j = l + 1; j <= n; j++)
                {
                    //
                    //  Apply the transformation.
                    //
                    if (l <= nct && s[l - 1] != 0.0)
                    {
                        t = -BLAS1D.ddot(m - l + 1, a, 1, a, 1,  + l - 1 + (l - 1) * lda,  + l - 1 + (j - 1) * lda)
                            / a[l - 1 + (l - 1) * lda];
                        BLAS1D.daxpy(m - l + 1, t, a, 1, ref a, 1,  + l - 1 + (l - 1) * lda,  + l - 1 + (j - 1) * lda);
                    }

                    //
                    //  Place the L-th row of A into E for the
                    //  subsequent calculation of the row transformation.
                    //
                    e[j - 1] = a[l - 1 + (j - 1) * lda];
                }

                //
                //  Place the transformation in U for subsequent back multiplication.
                //
                if (wantu && l <= nct)
                {
                    for (i = l; i <= m; i++)
                    {
                        u[i - 1 + (l - 1) * ldu] = a[i - 1 + (l - 1) * lda];
                    }
                }

                if (l <= nrt)
                {
                    //
                    //  Compute the L-th row transformation and place the
                    //  L-th superdiagonal in E(L).
                    //
                    e[l - 1] = BLAS1D.dnrm2(n - l, e, 1,  + l);

                    if (e[l - 1] != 0.0)
                    {
                        if (e[l] != 0.0)
                        {
                            e[l - 1] = typeMethods.r8_sign(e[l]) * Math.Abs(e[l - 1]);
                        }

                        BLAS1D.dscal(n - l, 1.0 / e[l - 1], ref e, 1,  + l);
                        e[l] = 1.0 + e[l];
                    }

                    e[l - 1] = -e[l - 1];
                    //
                    //  Apply the transformation.
                    //
                    if (l + 1 <= m && e[l - 1] != 0.0)
                    {
                        for (j = l + 1; j <= m; j++)
                        {
                            work[j - 1] = 0.0;
                        }

                        for (j = l + 1; j <= n; j++)
                        {
                            BLAS1D.daxpy(m - l, e[j - 1], a, 1, ref work, 1,  + l + (j - 1) * lda,  + l);
                        }

                        for (j = l + 1; j <= n; j++)
                        {
                            BLAS1D.daxpy(m - l, -e[j - 1] / e[l], work, 1, ref a, 1,  + l,  + l + (j - 1) * lda);
                        }
                    }

                    //
                    //  Place the transformation in V for subsequent back multiplication.
                    //
                    if (wantv)
                    {
                        for (j = l + 1; j <= n; j++)
                        {
                            v[j - 1 + (l - 1) * ldv] = e[j - 1];
                        }
                    }
                }
            }

            //
            //  Set up the final bidiagonal matrix of order MN.
            //
            mn = Math.Min(m + 1, n);
            nctp1 = nct + 1;
            nrtp1 = nrt + 1;

            if (nct < n)
            {
                s[nctp1 - 1] = a[nctp1 - 1 + (nctp1 - 1) * lda];
            }

            if (m < mn)
            {
                s[mn - 1] = 0.0;
            }

            if (nrtp1 < mn)
            {
                e[nrtp1 - 1] = a[nrtp1 - 1 + (mn - 1) * lda];
            }

            e[mn - 1] = 0.0;
            //
            //  If required, generate U.
            //
            if (wantu)
            {
                for (i = 1; i <= m; i++)
                {
                    for (j = nctp1; j <= ncu; j++)
                    {
                        u[(i - 1) + (j - 1) * ldu] = 0.0;
                    }
                }

                for (j = nctp1; j <= ncu; j++)
                {
                    u[j - 1 + (j - 1) * ldu] = 1.0;
                }

                for (ll = 1; ll <= nct; ll++)
                {
                    l = nct - ll + 1;

                    if (s[l - 1] != 0.0)
                    {
                        for (j = l + 1; j <= ncu; j++)
                        {
                            t = -BLAS1D.ddot(m - l + 1, u, 1, u, 1,  + (l - 1) + (l - 1) * ldu,  + (l - 1) + (j - 1) * ldu)
                                / u[l - 1 + (l - 1) * ldu];
                            BLAS1D.daxpy(m - l + 1, t, u, 1, ref u, 1,  + (l - 1) + (l - 1) * ldu,  + (l - 1) + (j - 1) * ldu);
                        }

                        BLAS1D.dscal(m - l + 1, -1.0, ref u, 1,  + (l - 1) + (l - 1) * ldu);
                        u[l - 1 + (l - 1) * ldu] = 1.0 + u[l - 1 + (l - 1) * ldu];
                        for (i = 1; i <= l - 1; i++)
                        {
                            u[i - 1 + (l - 1) * ldu] = 0.0;
                        }
                    }
                    else
                    {
                        for (i = 1; i <= m; i++)
                        {
                            u[i - 1 + (l - 1) * ldu] = 0.0;
                        }

                        u[l - 1 + (l - 1) * ldu] = 1.0;
                    }
                }
            }

            //
            //  If it is required, generate V.
            //
            if (wantv)
            {
                for (ll = 1; ll <= n; ll++)
                {
                    l = n - ll + 1;

                    if (l <= nrt && e[l - 1] != 0.0)
                    {
                        for (j = l + 1; j <= n; j++)
                        {
                            t = -BLAS1D.ddot(n - l, v, 1, v, 1,  + l + (l - 1) * ldv,  + l + (j - 1) * ldv)
                                / v[l + (l - 1) * ldv];
                            BLAS1D.daxpy(n - l, t, v, 1, ref v, 1,  + l + (l - 1) * ldv,  + l + (j - 1) * ldv);
                        }

                    }

                    for (i = 1; i <= n; i++)
                    {
                        v[i - 1 + (l - 1) * ldv] = 0.0;
                    }

                    v[l - 1 + (l - 1) * ldv] = 1.0;
                }
            }

            //
            //  Main iteration loop for the singular values.
            //
            mm = mn;
            iter = 0;

            while (0 < mn)
            {
                //
                //  If too many iterations have been performed, set flag and return.
                //
                if (maxit <= iter)
                {
                    info = mn;
                    return info;
                }

                //
                //  This section of the program inspects for
                //  negligible elements in the S and E arrays.
                //
                //  On completion the variables KASE and L are set as follows:
                //
                //  KASE = 1     if S(MN) and E(L-1) are negligible and L < MN
                //  KASE = 2     if S(L) is negligible and L < MN
                //  KASE = 3     if E(L-1) is negligible, L < MN, and
                //               S(L), ..., S(MN) are not negligible (QR step).
                //  KASE = 4     if E(MN-1) is negligible (convergence).
                //
                for (ll = 1; ll <= mn; ll++)
                {
                    l = mn - ll;

                    if (l == 0)
                    {
                        break;
                    }

                    test = Math.Abs(s[l - 1]) + Math.Abs(s[l]);
                    ztest = test + Math.Abs(e[l - 1]);

                    if (ztest == test)
                    {
                        e[l - 1] = 0.0;
                        break;
                    }
                }

                if (l == mn - 1)
                {
                    kase = 4;
                }
                else
                {
                    for (lls = l + 1; lls <= mn + 1; lls++)
                    {
                        ls = mn - lls + l + 1;

                        if (ls == l)
                        {
                            break;
                        }

                        test = 0.0;
                        if (ls != mn)
                        {
                            test = test + Math.Abs(e[ls - 1]);
                        }

                        if (ls != l + 1)
                        {
                            test = test + Math.Abs(e[ls - 2]);
                        }

                        ztest = test + Math.Abs(s[ls - 1]);

                        if (ztest == test)
                        {
                            s[ls - 1] = 0.0;
                            break;
                        }

                    }

                    if (ls == l)
                    {
                        kase = 3;
                    }
                    else if (ls == mn)
                    {
                        kase = 1;
                    }
                    else
                    {
                        kase = 2;
                        l = ls;
                    }
                }

                l = l + 1;
                //
                //  Deflate negligible S(MN).
                //
                if (kase == 1)
                {
                    mm1 = mn - 1;
                    f = e[mn - 2];
                    e[mn - 2] = 0.0;

                    for (kk = 1; kk <= mm1; kk++)
                    {
                        k = mm1 - kk + l;
                        t1 = s[k - 1];
                        BLAS1D.drotg(ref t1, ref f, ref cs, ref sn);
                        s[k - 1] = t1;

                        if (k != l)
                        {
                            f = -sn * e[k - 2];
                            e[k - 2] = cs * e[k - 2];
                        }

                        if (wantv)
                        {
                            BLAS1D.drot(n, ref v, 1, ref v, 1, cs, sn,  + 0 + (k - 1) * ldv,  + 0 + (mn - 1) * ldv);
                        }
                    }
                }
                //
                //  Split at negligible S(L).
                //
                else if (kase == 2)
                {
                    f = e[l - 2];
                    e[l - 2] = 0.0;

                    for (k = l; k <= mn; k++)
                    {
                        t1 = s[k - 1];
                        BLAS1D.drotg(ref t1, ref f, ref cs, ref sn);
                        s[k - 1] = t1;
                        f = -sn * e[k - 1];
                        e[k - 1] = cs * e[k - 1];
                        if (wantu)
                        {
                            BLAS1D.drot(m, ref u, 1, ref u, 1, cs, sn,  + 0 + (k - 1) * ldu,  + 0 + (l - 2) * ldu);
                        }
                    }
                }
                //
                //  Perform one QR step.
                //
                else if (kase == 3)
                {
                    //
                    //  Calculate the shift.
                    //
                    scale = Math.Max(Math.Abs(s[mn - 1]),
                        Math.Max(Math.Abs(s[mn - 2]),
                            Math.Max(Math.Abs(e[mn - 2]),
                                Math.Max(Math.Abs(s[l - 1]), Math.Abs(e[l - 1])))));

                    sm = s[mn - 1] / scale;
                    smm1 = s[mn - 2] / scale;
                    emm1 = e[mn - 2] / scale;
                    sl = s[l - 1] / scale;
                    el = e[l - 1] / scale;
                    b = ((smm1 + sm) * (smm1 - sm) + emm1 * emm1) / 2.0;
                    c = (sm * emm1) * (sm * emm1);
                    shift = 0.0;

                    if (b != 0.0 || c != 0.0)
                    {
                        shift = Math.Sqrt(b * b + c);
                        if (b < 0.0)
                        {
                            shift = -shift;
                        }

                        shift = c / (b + shift);
                    }

                    f = (sl + sm) * (sl - sm) - shift;
                    g = sl * el;
                    //
                    //  Chase zeros.
                    //
                    mm1 = mn - 1;

                    for (k = l; k <= mm1; k++)
                    {
                        BLAS1D.drotg(ref f, ref g, ref cs, ref sn);

                        if (k != l)
                        {
                            e[k - 2] = f;
                        }

                        f = cs * s[k - 1] + sn * e[k - 1];
                        e[k - 1] = cs * e[k - 1] - sn * s[k - 1];
                        g = sn * s[k];
                        s[k] = cs * s[k];

                        if (wantv)
                        {
                            BLAS1D.drot(n, ref v, 1, ref v, 1, cs, sn,  + 0 + (k - 1) * ldv,  + 0 + k * ldv);
                        }

                        BLAS1D.drotg(ref f, ref g, ref cs, ref sn);
                        s[k - 1] = f;
                        f = cs * e[k - 1] + sn * s[k];
                        s[k] = -sn * e[k - 1] + cs * s[k];
                        g = sn * e[k];
                        e[k] = cs * e[k];

                        if (wantu && k < m)
                        {
                            BLAS1D.drot(m, ref u, 1, ref u, 1, cs, sn,  + 0 + (k - 1) * ldu,  + 0 + k * ldu);
                        }
                    }

                    e[mn - 2] = f;
                    iter = iter + 1;
                }
                //
                //  Convergence.
                //
                else if (kase == 4)
                {
                    //
                    //  Make the singular value nonnegative.
                    //
                    if (s[l - 1] < 0.0)
                    {
                        s[l - 1] = -s[l - 1];
                        if (wantv)
                        {
                            BLAS1D.dscal(n, -1.0, ref v, 1,  + 0 + (l - 1) * ldv);
                        }
                    }

                    //
                    //  Order the singular value.
                    //
                    for (;;)
                    {
                        if (l == mm)
                        {
                            break;
                        }

                        if (s[l] <= s[l - 1])
                        {
                            break;
                        }

                        t = s[l - 1];
                        s[l - 1] = s[l];
                        s[l] = t;

                        if (wantv && l < n)
                        {
                            BLAS1D.dswap(n, ref v, 1, ref v, 1,  + 0 + (l - 1) * ldv,  + 0 + l * ldv);
                        }

                        if (wantu && l < m)
                        {
                            BLAS1D.dswap(m, ref u, 1, ref u, 1,  + 0 + (l - 1) * ldu,  + 0 + l * ldu);
                        }

                        l = l + 1;
                    }

                    iter = 0;
                    mn = mn - 1;
                }
            }

            return info;
        }
    }
}