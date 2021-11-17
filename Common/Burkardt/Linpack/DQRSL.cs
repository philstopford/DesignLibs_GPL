using System;
using Burkardt.BLAS;

namespace Burkardt.Linpack;

public static class DQRSL
{
    public static int dqrsl(double[] a, int lda, int n, int k, double[] qraux, double[] y,
            ref double[] qy, ref double[] qty, ref double[] b, ref double[] rsd, ref double[] ab, int job)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DQRSL computes transformations, projections, and least squares solutions.
        //
        //  Discussion:
        //
        //    DQRSL requires the output of DQRDC.
        //
        //    For K <= min(N,P), let AK be the matrix
        //
        //      AK = ( A(JPVT[0]), A(JPVT(2)), ..., A(JPVT(K)) )
        //
        //    formed from columns JPVT[0], ..., JPVT(K) of the original
        //    N by P matrix A that was input to DQRDC.  If no pivoting was
        //    done, AK consists of the first K columns of A in their
        //    original order.  DQRDC produces a factored orthogonal matrix Q
        //    and an upper triangular matrix R such that
        //
        //      AK = Q * (R)
        //               (0)
        //
        //    This information is contained in coded form in the arrays
        //    A and QRAUX.
        //
        //    The parameters QY, QTY, B, RSD, and AB are not referenced
        //    if their computation is not requested and in this case
        //    can be replaced by dummy variables in the calling program.
        //    To save storage, the user may in some cases use the same
        //    array for different parameters in the calling sequence.  A
        //    frequently occuring example is when one wishes to compute
        //    any of B, RSD, or AB and does not need Y or QTY.  In this
        //    case one may identify Y, QTY, and one of B, RSD, or AB, while
        //    providing separate arrays for anything else that is to be
        //    computed.
        //
        //    Thus the calling sequence
        //
        //      dqrsl ( a, lda, n, k, qraux, y, dum, y, b, y, dum, 110, info )
        //
        //    will result in the computation of B and RSD, with RSD
        //    overwriting Y.  More generally, each item in the following
        //    list contains groups of permissible identifications for
        //    a single calling sequence.
        //
        //      1. (Y,QTY,B) (RSD) (AB) (QY)
        //
        //      2. (Y,QTY,RSD) (B) (AB) (QY)
        //
        //      3. (Y,QTY,AB) (B) (RSD) (QY)
        //
        //      4. (Y,QY) (QTY,B) (RSD) (AB)
        //
        //      5. (Y,QY) (QTY,RSD) (B) (AB)
        //
        //      6. (Y,QY) (QTY,AB) (B) (RSD)
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
        //    Input, double A[LDA*P], contains the output of DQRDC.
        //
        //    Input, int LDA, the leading dimension of the array A.
        //
        //    Input, int N, the number of rows of the matrix AK.  It must
        //    have the same value as N in DQRDC.
        //
        //    Input, int K, the number of columns of the matrix AK.  K
        //    must not be greater than min(N,P), where P is the same as in the
        //    calling sequence to DQRDC.
        //
        //    Input, double QRAUX[P], the auxiliary output from DQRDC.
        //
        //    Input, double Y[N], a vector to be manipulated by DQRSL.
        //
        //    Output, double QY[N], contains Q * Y, if requested.
        //
        //    Output, double QTY[N], contains Q' * Y, if requested.
        //
        //    Output, double B[K], the solution of the least squares problem
        //      minimize norm2 ( Y - AK * B),
        //    if its computation has been requested.  Note that if pivoting was
        //    requested in DQRDC, the J-th component of B will be associated with
        //    column JPVT(J) of the original matrix A that was input into DQRDC.
        //
        //    Output, double RSD[N], the least squares residual Y - AK * B,
        //    if its computation has been requested.  RSD is also the orthogonal
        //    projection of Y onto the orthogonal complement of the column space
        //    of AK.
        //
        //    Output, double AB[N], the least squares approximation Ak * B,
        //    if its computation has been requested.  AB is also the orthogonal
        //    projection of Y onto the column space of A.
        //
        //    Input, integer JOB, specifies what is to be computed.  JOB has
        //    the decimal expansion ABCDE, with the following meaning:
        //
        //      if A != 0, compute QY.
        //      if B != 0, compute QTY.
        //      if C != 0, compute QTY and B.
        //      if D != 0, compute QTY and RSD.
        //      if E != 0, compute QTY and AB.
        //
        //    Note that a request to compute B, RSD, or AB automatically triggers
        //    the computation of QTY, for which an array must be provided in the
        //    calling sequence.
        //
        //    Output, int DQRSL, is zero unless the computation of B has
        //    been requested and R is exactly singular.  In this case, INFO is the
        //    index of the first zero diagonal element of R, and B is left unaltered.
        //
    {
        bool cab;
        bool cb;
        bool cqty;
        bool cqy;
        bool cr;
        int i;
        int info;
        int j;
        int jj;
        int ju;
        double t;
        double temp;
        //
        //  set info flag.
        //
        info = 0;
        //
        //  Determine what is to be computed.
        //
        cqy = job / 10000 != 0;
        cqty = job % 10000 != 0;
        cb = job % 1000 / 100 != 0;
        cr = job % 100 / 10 != 0;
        cab = job % 10 != 0;

        ju = Math.Min(k, n - 1);
        switch (ju)
        {
            //
            //  Special action when N = 1.
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

                ab[0] = cab switch
                {
                    true => y[0],
                    _ => ab[0]
                };

                switch (cb)
                {
                    case true when a[0 + 0 * lda] == 0.0:
                        info = 1;
                        break;
                    case true:
                        b[0] = y[0] / a[0 + 0 * lda];
                        break;
                }

                rsd[0] = cr switch
                {
                    true => 0.0,
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
                for (i = 1; i <= n; i++)
                {
                    qy[i - 1] = y[i - 1];
                }

                break;
            }
        }

        switch (cqty)
        {
            case true:
            {
                for (i = 1; i <= n; i++)
                {
                    qty[i - 1] = y[i - 1];
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

                    if (qraux[j - 1] != 0.0)
                    {
                        temp = a[j - 1 + (j - 1) * lda];
                        a[j - 1 + (j - 1) * lda] = qraux[j - 1];
                        t = -BLAS1D.ddot(n - j + 1, a, 1, qy, 1, xIndex: +j - 1 + (j - 1) * lda, yIndex: +j - 1) /
                            a[j - 1 + (j - 1) * lda];
                        BLAS1D.daxpy(n - j + 1, t, a, 1, ref qy, 1, xIndex: +j - 1 + (j - 1) * lda, yIndex: +j - 1);
                        a[j - 1 + (j - 1) * lda] = temp;
                    }
                }

                break;
            }
        }

        switch (cqty)
        {
            //
            //  Compute Q'*Y.
            //
            case true:
            {
                for (j = 1; j <= ju; j++)
                {
                    if (qraux[j - 1] != 0.0)
                    {
                        temp = a[j - 1 + (j - 1) * lda];
                        a[j - 1 + (j - 1) * lda] = qraux[j - 1];
                        t = -BLAS1D.ddot(n - j + 1, a, 1, qty, 1, xIndex: +j - 1 + (j - 1) * lda, yIndex: +j - 1) /
                            a[j - 1 + (j - 1) * lda];
                        BLAS1D.daxpy(n - j + 1, t, a, 1, ref qty, 1, xIndex: +j - 1 + (j - 1) * lda, yIndex: +j - 1);
                        a[j - 1 + (j - 1) * lda] = temp;
                    }
                }

                break;
            }
        }

        switch (cb)
        {
            //
            //  Set up to compute B, RSD, or AB.
            //
            case true:
            {
                for (i = 1; i <= k; i++)
                {
                    b[i - 1] = qty[i - 1];
                }

                break;
            }
        }

        switch (cab)
        {
            case true:
            {
                for (i = 1; i <= k; i++)
                {
                    ab[i - 1] = qty[i - 1];
                }

                break;
            }
        }

        switch (cr)
        {
            case true when k < n:
            {
                for (i = k + 1; i <= n; i++)
                {
                    rsd[i - 1] = qty[i - 1];
                }

                break;
            }
        }

        switch (cab)
        {
            case true when k + 1 <= n:
            {
                for (i = k + 1; i <= n; i++)
                {
                    ab[i - 1] = 0.0;
                }

                break;
            }
        }

        switch (cr)
        {
            case true:
            {
                for (i = 1; i <= k; i++)
                {
                    rsd[i - 1] = 0.0;
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

                    if (a[j - 1 + (j - 1) * lda] == 0.0)
                    {
                        info = j;
                        break;
                    }

                    b[j - 1] /= a[j - 1 + (j - 1) * lda];

                    if (j != 1)
                    {
                        t = -b[j - 1];
                        BLAS1D.daxpy(j - 1, t, a, 1, ref b, 1, xIndex: +0 + (j - 1) * lda);
                    }
                }

                break;
            }
        }

        //
        //  Compute RSD or AB as required.
        //
        if (cr || cab)
        {
            for (jj = 1; jj <= ju; jj++)
            {
                j = ju - jj + 1;

                if (qraux[j - 1] != 0.0)
                {
                    temp = a[j - 1 + (j - 1) * lda];
                    a[j - 1 + (j - 1) * lda] = qraux[j - 1];

                    switch (cr)
                    {
                        case true:
                            t = -BLAS1D.ddot(n - j + 1, a, 1, rsd, 1, xIndex: +j - 1 + (j - 1) * lda, yIndex: +j - 1)
                                / a[j - 1 + (j - 1) * lda];
                            BLAS1D.daxpy(n - j + 1, t, a, 1, ref rsd, 1, xIndex: +j - 1 + (j - 1) * lda,
                                yIndex: +j - 1);
                            break;
                    }

                    switch (cab)
                    {
                        case true:
                            t = -BLAS1D.ddot(n - j + 1, a, 1, ab, 1, xIndex: +j - 1 + (j - 1) * lda, yIndex: +j - 1)
                                / a[j - 1 + (j - 1) * lda];
                            BLAS1D.daxpy(n - j + 1, t, a, 1, ref ab, 1, xIndex: +j - 1 + (j - 1) * lda, yIndex: +j - 1);
                            break;
                    }

                    a[j - 1 + (j - 1) * lda] = temp;
                }
            }
        }

        return info;
    }
}