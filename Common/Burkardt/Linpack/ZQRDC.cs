using System;
using System.Numerics;
using Burkardt.BLAS;
using Burkardt.Types;

namespace Burkardt.Linpack;

public static class ZQRDC
{
    public static void zqrdc(ref Complex[] x, int ldx, int n, int p,
            ref Complex[] qraux, ref int[] ipvt, int job )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZQRDC computes the QR factorization of an N by P complex <double> matrix.
        //
        //  Discussion:
        //
        //    ZQRDC uses Householder transformations to compute the QR factorization
        //    of an N by P matrix X.  Column pivoting based on the 2-norms of the
        //    reduced columns may be performed at the user's option.
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
        //    Input/output, complex <double> X[LDX*P]; on input, the matrix whose decomposition
        //    is to be computed.  On output, the upper triangle contains the upper
        //    triangular matrix R of the QR factorization.  Below its diagonal, X
        //    contains information from which the unitary part of the decomposition
        //    can be recovered.  If pivoting has been requested, the decomposition is
        //    not that of the original matrix X, but that of X with its columns
        //    permuted as described by IPVT.
        //
        //    Input, int LDX, the leading dimension of X.  N <= LDX.
        //
        //    Input, int N, the number of rows of the matrix.
        //
        //    Input, int P, the number of columns in the matrix X.
        //
        //    Output, complex <double> QRAUX[P], further information required to recover
        //    the unitary part of the decomposition.
        //
        //    Input/output, int IPVT[P]; on input, ints that control the
        //    selection of the pivot columns.  The K-th column X(K) of X is placed
        //    in one of three classes according to the value of IPVT(K):
        //      IPVT(K) > 0, then X(K) is an initial column.
        //      IPVT(K) == 0, then X(K) is a free column.
        //      IPVT(K) < 0, then X(K) is a final column.
        //    Before the decomposition is computed, initial columns are moved to the
        //    beginning of the array X and final columns to the end.  Both initial
        //    and final columns are frozen in place during the computation and only
        //    free columns are moved.  At the K-th stage of the reduction, if X(K)
        //    is occupied by a free column it is interchanged with the free column
        //    of largest reduced norm.
        //    On output, IPVT(K) contains the index of the column of the
        //    original matrix that has been interchanged into
        //    the K-th column, if pivoting was requested.
        //    IPVT is not referenced if JOB == 0.
        //
        //    Input, int JOB, initiates column pivoting.
        //    0, no pivoting is done.
        //    nonzero, pivoting is done.
        //
    {
        int itemp;
        int j;
        int l;

        int pl = 1;
        int pu = 0;
        Complex[] work = new Complex [p];

        if (job != 0)
        {
            //
            //  Pivoting has been requested.  Rearrange the columns according to IPVT.
            //
            for (j = 1; j <= p; j++)
            {
                bool swapj = 0 < ipvt[j - 1];
                bool negj = ipvt[j - 1] < 0;

                ipvt[j - 1] = negj switch
                {
                    true => -j,
                    _ => j
                };

                switch (swapj)
                {
                    case true:
                    {
                        if (j != pl)
                        {
                            BLAS1Z.zswap(n, ref x, 1, ref x, 1, xIndex: + 0 + (pl - 1) * ldx, yIndex: + 0 + (j - 1) * ldx);
                        }

                        ipvt[j - 1] = ipvt[pl - 1];
                        ipvt[pl - 1] = j;
                        pl += 1;
                        break;
                    }
                }
            }

            pu = p;

            int jj;
            for (jj = 1; jj <= p; jj++)
            {
                j = p - jj + 1;

                switch (ipvt[j - 1])
                {
                    case < 0:
                    {
                        ipvt[j - 1] = -ipvt[j - 1];

                        if (j != pu)
                        {
                            BLAS1Z.zswap(n, ref x, 1, ref x, 1, xIndex: + 0 + (pu - 1) * ldx, yIndex: + 0 + (j - 1) * ldx);

                            itemp = ipvt[pu - 1];
                            ipvt[pu - 1] = ipvt[j - 1];
                            ipvt[j - 1] = itemp;
                        }

                        pu -= 1;
                        break;
                    }
                }
            }
        }

        //
        //  Compute the norms of the free columns.
        //
        for (j = pl; j <= pu; j++)
        {
            qraux[j - 1] = new Complex(BLAS1Z.dznrm2(n, x, 1, index: + 0 + (j - 1) * ldx), 0.0);
            work[j - 1] = qraux[j - 1];
        }

        //
        //  Perform the Householder reduction of X.
        //
        int lup = Math.Min(n, p);

        for (l = 1; l <= lup; l++)
        {
            //
            //  Locate the column of largest norm and bring it
            //  into the pivot position.
            //
            if (pl <= l && l < pu)
            {
                double maxnrm = 0.0;
                int maxj = l;

                for (j = l; j <= pu; j++)
                {
                    if (!(maxnrm < qraux[j - 1].Real))
                    {
                        continue;
                    }

                    maxnrm = qraux[j - 1].Real;
                    maxj = j;
                }

                if (maxj != l)
                {
                    BLAS1Z.zswap(n, ref x, 1, ref x, 1, xIndex: + 0 + (l - 1) * ldx, yIndex: + 0 + (maxj - 1) * ldx);
                    qraux[maxj - 1] = qraux[l - 1];
                    work[maxj - 1] = work[l - 1];

                    itemp = ipvt[maxj - 1];
                    ipvt[maxj - 1] = ipvt[l - 1];
                    ipvt[l - 1] = itemp;
                }
            }

            qraux[l - 1] = new Complex(0.0, 0.0);

            if (l == n)
            {
                continue;
            }

            //
            //  Compute the Householder transformation for column L.
            //
            Complex nrmxl = new(BLAS1Z.dznrm2(n - l + 1, x, 1, index: + l - 1 + (l - 1) * ldx), 0.0);

            if (typeMethods.zabs1(nrmxl) == 0.0)
            {
                continue;
            }

            if (typeMethods.zabs1(x[l - 1 + (l - 1) * ldx]) != 0.0)
            {
                nrmxl = typeMethods.zsign2(nrmxl, x[l - 1 + (l - 1) * ldx]);
            }

            Complex t = new Complex(1.0, 0.0) / nrmxl;
            BLAS1Z.zscal(n - l + 1, t, ref x, 1, index: + l - 1 + (l - 1) * ldx);
            x[l - 1 + (l - 1) * ldx] = new Complex(1.0, 0.0) + x[l - 1 + (l - 1) * ldx];
            //
            //  Apply the transformation to the remaining columns,
            //  updating the norms.
            //
            for (j = l + 1; j <= p; j++)
            {
                t = -BLAS1Z.zdotc(n - l + 1, x, 1, x, 1, xIndex: + l - 1 + (l - 1) * ldx, yIndex: + l - 1 + (j - 1) * ldx)
                    / x[l - 1 + (l - 1) * ldx];
                BLAS1Z.zaxpy(n - l + 1, t, x, 1, ref x, 1, xIndex: + l - 1 + (l - 1) * ldx, yIndex: + l - 1 + (j - 1) * ldx);

                if (j < pl || pu < j)
                {
                    continue;
                }

                if (typeMethods.zabs1(qraux[j - 1]) == 0.0)
                {
                    continue;
                }

                double tt = 1.0 - Math.Pow(Complex.Abs(x[l - 1 + (j - 1) * ldx]) / qraux[j - 1].Real, 2);
                tt = Math.Max(tt, 0.0);
                t = new Complex(tt, 0.0);
                tt = 1.0 + 0.05 * tt
                                * Math.Pow(qraux[j - 1].Real / work[j - 1].Real, 2);

                if (Math.Abs(tt - 1.0) > typeMethods.r8_epsilon())
                {
                    qraux[j - 1] *= Complex.Sqrt(t);
                }
                else
                {
                    qraux[j - 1] =
                        new Complex(BLAS1Z.dznrm2(n - l, x, 1, index: + l + (j - 1) * ldx), 0.0);
                    work[j - 1] = qraux[j - 1];
                }
            }

            //
            //  Save the transformation.
            //
            qraux[l - 1] = x[l - 1 + (l - 1) * ldx];
            x[l - 1 + (l - 1) * ldx] = -nrmxl;
        }
    }

}