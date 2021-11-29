using System;
using System.Numerics;
using Burkardt.BLAS;
using Burkardt.Types;

namespace Burkardt.Linpack;

public static class ZQRSL
{
    public static int zqrsl(Complex[] x, int ldx, int n, int k,
            Complex[] qraux, Complex[] y, ref Complex[] qy,
            ref Complex[] qty, ref Complex[] b, ref Complex[] rsd,
            ref Complex[] xb, int job)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZQRSL solves, transforms or projects systems factored by ZQRDC.
        //
        //  Discussion:
        //
        //    The routine applies the output of ZQRDC to compute coordinate
        //    transformations, projections, and least squares solutions.
        //
        //    For K <= min ( N, P ), let XK be the matrix
        //
        //      XK = ( X(IPVT(1)), X(IPVT(2)), ... ,X(IPVT(k)) )
        //
        //    formed from columnns IPVT(1), ... ,IPVT(K) of the original
        //    N by P matrix X that was input to ZQRDC (if no pivoting was
        //    done, XK consists of the first K columns of X in their
        //    original order).  ZQRDC produces a factored unitary matrix Q
        //    and an upper triangular matrix R such that
        //
        //      XK = Q * ( R )
        //               ( 0 )
        //
        //    This information is contained in coded form in the arrays
        //    X and QRAUX.
        //
        //    The parameters QY, QTY, B, RSD, and XB are not referenced
        //    if their computation is not requested and in this case
        //    can be replaced by dummy variables in the calling program.
        //
        //    To save storage, the user may in some cases use the same
        //    array for different parameters in the calling sequence.  A
        //    frequently occuring example is when one wishes to compute
        //    any of B, RSD, or XB and does not need Y or QTY.  In this
        //    case one may identify Y, QTY, and one of B, RSD, or XB, while
        //    providing separate arrays for anything else that is to be
        //    computed.  Thus the calling sequence
        //
        //      zqrsl ( x, ldx, n, k, qraux, y, dum, y, b, y, dum, 110, info )
        //
        //    will result in the computation of B and RSD, with RSD
        //    overwriting Y.  More generally, each item in the following
        //    list contains groups of permissible identifications for
        //    a single callinng sequence.
        //
        //    1. ( Y, QTY, B )   ( RSD )      ( XB )  ( QY )
        //    2. ( Y, QTY, RSD ) ( B )        ( XB )  ( QY )
        //    3. ( Y, QTY, XB )  ( B )        ( RSD ) ( QY )
        //    4. ( Y, QY )       ( QTY, B )   ( RSD ) ( XB )
        //    5. ( Y, QY )       ( QTY, RSD ) ( B )   ( XB )
        //    6. ( Y, QY )       ( QTY, XB )  ( B )   ( RSD )
        //
        //    In any group the value returned in the array allocated to
        //    the group corresponds to the last member of the group.
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
        //    Input, Complex X[LDX*P], the output of ZQRDC.
        //
        //    Input, int LDX, the leading dimension of X.
        //
        //    Input, int N, the number of rows of the matrix XK, which
        //    must have the same value as N in ZQRDC.
        //
        //    Input, int K, the number of columns of the matrix XK.  K must not
        //    be greater than min ( N, P), where P is the same as in the calling
        //    sequence to ZQRDC.
        //
        //    Input, Complex QRAUX[P], the auxiliary output from ZQRDC.
        //
        //    Input, Complex Y[N], a vector that is to be manipulated by ZQRSL.
        //
        //    Output, Complex QY[N], contains Q*Y, if it has been requested.
        //
        //    Output, Complex QTY[N], contains hermitian(Q)*Y, if it has
        //    been requested.  Here hermitian(Q) is the conjugate transpose
        //    of the matrix Q.
        //
        //    Output, Complex B[K], the solution of the least squares problem
        //      minimize norm2 ( Y - XK * B ),
        //    if it has been requested.  If pivoting was requested in ZQRDC,
        //    the J-th component of B will be associated with column IPVT(J)
        //    of the original matrix X that was input into ZQRDC.
        //
        //    Output, Complex RSD[N], the least squares residual Y - XK*B,
        //    if it has been requested.  RSD is also the orthogonal projection
        //    of Y onto the orthogonal complement of the column space of XK.
        //
        //    Output, Complex XB[N], the least squares approximation XK*N,
        //    if its computation has been requested.  XB is also the orthogonal
        //    projection of Y onto the column space of X.
        //
        //    Input, int JOB, specifies what is to be computed.  JOB has
        //    the decimal expansion ABCDE, meaning:
        //    if A != 0, compute QY.
        //    if B, D, D, or E != 0, compute QTY.
        //    if C != 0, compute B.
        //    if D != 0, compute RSD.
        //    if E != 0, compute XB.
        //    A request to compute B, RSD, or XB automatically triggers the
        //    computation of QTY, for which an array must be provided in the
        //    calling sequence.
        //
        //    Output, int ZQRSL, the value of INFO, which is zero unless
        //    the computation of B has been requested and R is exactly singular.
        //    In this case, INFO is the index of the first zero diagonal element
        //    of R and B is left unaltered.
        //
    {
        int i;
        int j;
        int jj;
        Complex t;
        Complex temp;

        int info = 0;
        //
        //  Determine what is to be computed.
        //
        bool cqy = job / 10000 != 0;
        bool cqty = job % 10000 != 0;
        bool cb = job % 1000 / 100 != 0;
        bool cr = job % 100 / 10 != 0;
        bool cxb = job % 10 != 0;

        int ju = Math.Min(k, n - 1);
        switch (ju)
        {
            //
            //  Special action when N=1.
            //
            case 0:
            {
                qy[0] = cqy switch
                {
                    true => y[0],
                    _ => qy[0]
                };

                qty[0] = cqty switch
                {
                    true => y[0],
                    _ => qty[0]
                };

                xb[0] = cxb switch
                {
                    true => y[0],
                    _ => xb[0]
                };

                switch (cb)
                {
                    case true when typeMethods.zabs1(x[0 + 0 * ldx]) == 0.0:
                        info = 1;
                        break;
                    case true:
                        b[0] = y[0] / x[0 + 0 * ldx];
                        break;
                }

                rsd[0] = cr switch
                {
                    true => new Complex(0.0, 0.0),
                    _ => rsd[0]
                };

                return info;
            }
        }

        switch (cqy)
        {
            //
            //  Set up to compute QY or QTY.
            //
            case true:
            {
                for (i = 0; i < n; i++)
                {
                    qy[i] = y[i];
                }

                break;
            }
        }

        switch (cqty)
        {
            case true:
            {
                for (i = 0; i < n; i++)
                {
                    qty[i] = y[i];
                }

                break;
            }
        }

        switch (cqy)
        {
            //
            //  Compute QY.
            //
            case true:
            {
                for (jj = 1; jj <= ju; jj++)
                {
                    j = ju - jj + 1;

                    if (typeMethods.zabs1(qraux[j - 1]) == 0.0)
                    {
                        continue;
                    }

                    temp = x[j - 1 + (j - 1) * ldx];
                    x[j - 1 + (j - 1) * ldx] = qraux[j - 1];
                    t = -BLAS1Z.zdotc(n - j + 1, x, 1, qy, 1, xIndex: +j - 1 + (j - 1) * ldx, yIndex: +j - 1) /
                        x[j - 1 + (j - 1) * ldx];
                    BLAS1Z.zaxpy(n - j + 1, t, x, 1, ref qy, 1, xIndex: +j - 1 + (j - 1) * ldx, yIndex: +j - 1);
                    x[j - 1 + (j - 1) * ldx] = temp;
                }

                break;
            }
        }

        switch (cqty)
        {
            //
            //  Compute hermitian ( A ) * Y.
            //
            case true:
            {
                for (j = 1; j <= ju; j++)
                {
                    if (typeMethods.zabs1(qraux[j - 1]) == 0.0)
                    {
                        continue;
                    }

                    temp = x[j - 1 + (j - 1) * ldx];
                    x[j - 1 + (j - 1) * ldx] = qraux[j - 1];
                    t = -BLAS1Z.zdotc(n - j + 1, x, 1, qty, 1, xIndex: +j - 1 + (j - 1) * ldx, yIndex: +j - 1) /
                        x[j - 1 + (j - 1) * ldx];
                    BLAS1Z.zaxpy(n - j + 1, t, x, 1, ref qty, 1, xIndex: +j - 1 + (j - 1) * ldx, yIndex: +j - 1);
                    x[j - 1 + (j - 1) * ldx] = temp;
                }

                break;
            }
        }

        switch (cb)
        {
            //
            //  Set up to compute B, RSD, or XB.
            //
            case true:
            {
                for (i = 0; i < k; i++)
                {
                    b[i] = qty[i];
                }

                break;
            }
        }

        switch (cxb)
        {
            case true:
            {
                for (i = 0; i < k; i++)
                {
                    xb[i] = qty[i];
                }

                break;
            }
        }

        switch (cr)
        {
            case true when k < n:
            {
                for (i = k; i < n; i++)
                {
                    rsd[i] = qty[i];
                }

                break;
            }
        }

        switch (cxb)
        {
            case true:
            {
                for (i = k; i < n; i++)
                {
                    xb[i] = new Complex(0.0, 0.0);
                }

                break;
            }
        }

        switch (cr)
        {
            case true:
            {
                for (i = 0; i < k; i++)
                {
                    rsd[i] = new Complex(0.0, 0.0);
                }

                break;
            }
        }

        switch (cb)
        {
            //
            //  Compute B.
            //
            case true:
            {
                for (jj = 1; jj <= k; jj++)
                {
                    j = k - jj + 1;

                    if (typeMethods.zabs1(x[j - 1 + (j - 1) * ldx]) == 0.0)
                    {
                        info = j;
                        break;
                    }

                    b[j - 1] /= x[j - 1 + (j - 1) * ldx];

                    if (j == 1)
                    {
                        continue;
                    }

                    t = -b[j - 1];
                    BLAS1Z.zaxpy(j - 1, t, x, 1, ref b, 1, xIndex: +0 + (j - 1) * ldx);
                }

                break;
            }
        }

        if (!cr && !cxb)
        {
            return info;
        }

        //
        //  Compute RSD or XB as required.
        //
        for (jj = 1; jj <= ju; jj++)
        {
            j = ju - jj + 1;

            if (typeMethods.zabs1(qraux[j - 1]) == 0.0)
            {
                continue;
            }

            temp = x[j - 1 + (j - 1) * ldx];
            x[j - 1 + (j - 1) * ldx] = qraux[j - 1];

            switch (cr)
            {
                case true:
                    t = -BLAS1Z.zdotc(n - j + 1, x, 1, rsd, 1, xIndex: +j - 1 + (j - 1) * ldx, yIndex: +j - 1)
                        / x[j - 1 + (j - 1) * ldx];
                    BLAS1Z.zaxpy(n - j + 1, t, x, 1, ref rsd, 1, xIndex: +j - 1 + (j - 1) * ldx,
                        yIndex: +j - 1);
                    break;
            }

            switch (cxb)
            {
                case true:
                    t = -BLAS1Z.zdotc(n - j + 1, x, 1, xb, 1, xIndex: +j - 1 + (j - 1) * ldx, yIndex: +j - 1)
                        / x[j - 1 + (j - 1) * ldx];
                    BLAS1Z.zaxpy(n - j + 1, t, x, 1, ref xb, 1, xIndex: +j - 1 + (j - 1) * ldx, yIndex: +j - 1);
                    break;
            }

            x[j - 1 + (j - 1) * ldx] = temp;
        }

        return info;
    }

}