using System;
using Burkardt.BLAS;
using Burkardt.Types;

namespace Burkardt.Linpack;

public static class DQRDC
{
    public static void dqrdc(ref double[] a, int lda, int n, int p, ref double[] qraux, ref int[] jpvt,
            double[] work, int job)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DQRDC computes the QR factorization of a real rectangular matrix.
        //
        //  Discussion:
        //
        //    DQRDC uses Householder transformations.
        //
        //    Column pivoting based on the 2-norms of the reduced columns may be
        //    performed at the user's option.
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
        //    Input/output, double A(LDA,P).  On input, the N by P matrix
        //    whose decomposition is to be computed.  On output, A contains in
        //    its upper triangle the upper triangular matrix R of the QR
        //    factorization.  Below its diagonal A contains information from
        //    which the orthogonal part of the decomposition can be recovered.
        //    Note that if pivoting has been requested, the decomposition is not that
        //    of the original matrix A but that of A with its columns permuted
        //    as described by JPVT.
        //
        //    Input, int LDA, the leading dimension of the array A.  LDA must
        //    be at least N.
        //
        //    Input, int N, the number of rows of the matrix A.
        //
        //    Input, int P, the number of columns of the matrix A.
        //
        //    Output, double QRAUX[P], contains further information required
        //    to recover the orthogonal part of the decomposition.
        //
        //    Input/output, integer JPVT[P].  On input, JPVT contains integers that
        //    control the selection of the pivot columns.  The K-th column A(*,K) of A
        //    is placed in one of three classes according to the value of JPVT(K).
        //      > 0, then A(K) is an initial column.
        //      = 0, then A(K) is a free column.
        //      < 0, then A(K) is a final column.
        //    Before the decomposition is computed, initial columns are moved to
        //    the beginning of the array A and final columns to the end.  Both
        //    initial and final columns are frozen in place during the computation
        //    and only free columns are moved.  At the K-th stage of the
        //    reduction, if A(*,K) is occupied by a free column it is interchanged
        //    with the free column of largest reduced norm.  JPVT is not referenced
        //    if JOB == 0.  On output, JPVT(K) contains the index of the column of the
        //    original matrix that has been interchanged into the K-th column, if
        //    pivoting was requested.
        //
        //    Workspace, double WORK[P].  WORK is not referenced if JOB == 0.
        //
        //    Input, int JOB, initiates column pivoting.
        //    0, no pivoting is done.
        //    nonzero, pivoting is done.
        //
    {
        int j;
        int jp;
        int l;
        int lup;
        int maxj;
        double maxnrm;
        double nrmxl;
        int pl;
        int pu;
        bool swapj;
        double t;
        double tt;

        pl = 1;
        pu = 0;
        //
        //  If pivoting is requested, rearrange the columns.
        //
        if (job != 0)
        {
            for (j = 1; j <= p; j++)
            {
                swapj = 0 < jpvt[j - 1];

                jpvt[j - 1] = jpvt[j - 1] switch
                {
                    < 0 => -j,
                    _ => j
                };

                switch (swapj)
                {
                    case true:
                    {
                        if (j != pl)
                        {
                            BLAS1D.dswap(n, ref a, 1, ref a, 1, xIndex: +0 + (pl - 1) * lda, yIndex: +0 + (j - 1));
                        }

                        jpvt[j - 1] = jpvt[pl - 1];
                        jpvt[pl - 1] = j;
                        pl += 1;
                        break;
                    }
                }
            }

            pu = p;

            for (j = p; 1 <= j; j--)
            {
                switch (jpvt[j - 1])
                {
                    case < 0:
                    {
                        jpvt[j - 1] = -jpvt[j - 1];

                        if (j != pu)
                        {
                            BLAS1D.dswap(n, ref a, 1, ref a, 1, xIndex: +0 + (pu - 1) * lda,
                                yIndex: +0 + (j - 1) * lda);
                            jp = jpvt[pu - 1];
                            jpvt[pu - 1] = jpvt[j - 1];
                            jpvt[j - 1] = jp;
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
            qraux[j - 1] = BLAS1D.dnrm2(n, a, 1, index: +0 + (j - 1) * lda);
        }

        for (j = pl; j <= pu; j++)
        {
            work[j - 1] = qraux[j - 1];
        }

        //
        //  Perform the Householder reduction of A.
        //
        lup = Math.Min(n, p);

        for (l = 1; l <= lup; l++)
        {
            //
            //  Bring the column of largest norm into the pivot position.
            //
            if (pl <= l && l < pu)
            {
                maxnrm = 0.0;
                maxj = l;
                for (j = l; j <= pu; j++)
                {
                    if (maxnrm < qraux[j - 1])
                    {
                        maxnrm = qraux[j - 1];
                        maxj = j;
                    }
                }

                if (maxj != l)
                {
                    BLAS1D.dswap(n, ref a, 1, ref a, 1, xIndex: +0 + (l - 1) * lda, yIndex: +0 + (maxj - 1) * lda);
                    qraux[maxj - 1] = qraux[l - 1];
                    work[maxj - 1] = work[l - 1];
                    jp = jpvt[maxj - 1];
                    jpvt[maxj - 1] = jpvt[l - 1];
                    jpvt[l - 1] = jp;
                }
            }

            //
            //  Compute the Householder transformation for column L.
            //
            qraux[l - 1] = 0.0;

            if (l != n)
            {
                nrmxl = BLAS1D.dnrm2(n - l + 1, a, 1, index: +l - 1 + (l - 1) * lda);

                if (nrmxl != 0.0)
                {
                    if (a[l - 1 + (l - 1) * lda] != 0.0)
                    {
                        nrmxl *= typeMethods.r8_sign(a[l - 1 + (l - 1) * lda]);
                    }

                    BLAS1D.dscal(n - l + 1, 1.0 / nrmxl, ref a, 1, index: +l - 1 + (l - 1) * lda);
                    a[l - 1 + (l - 1) * lda] = 1.0 + a[l - 1 + (l - 1) * lda];
                    //
                    //  Apply the transformation to the remaining columns, updating the norms.
                    //
                    for (j = l + 1; j <= p; j++)
                    {
                        t = -BLAS1D.ddot(n - l + 1, a, 1, a, 1, xIndex: +l - 1 + (l - 1) * lda,
                                yIndex: +l - 1 + (j - 1) * lda)
                            / a[l - 1 + (l - 1) * lda];
                        BLAS1D.daxpy(n - l + 1, t, a, 1, ref a, 1, xIndex: +l - 1 + (l - 1) * lda,
                            yIndex: +l - 1 + (j - 1) * lda);

                        if (pl <= j && j <= pu)
                        {
                            if (qraux[j - 1] != 0.0)
                            {
                                tt = 1.0 - Math.Pow(Math.Abs(a[l - 1 + (j - 1) * lda]) / qraux[j - 1], 2);
                                tt = Math.Max(tt, 0.0);
                                t = tt;
                                tt = 1.0 + 0.05 * tt * Math.Pow(qraux[j - 1] / work[j - 1], 2);

                                if (Math.Abs(tt - 1.0) > double.Epsilon)
                                {
                                    qraux[j - 1] *= Math.Sqrt(t);
                                }
                                else
                                {
                                    qraux[j - 1] = BLAS1D.dnrm2(n - l, a, 1, index: +l + (j - 1) * lda);
                                    work[j - 1] = qraux[j - 1];
                                }
                            }
                        }
                    }

                    //
                    //  Save the transformation.
                    //
                    qraux[l - 1] = a[l - 1 + (l - 1) * lda];
                    a[l - 1 + (l - 1) * lda] = -nrmxl;
                }
            }
        }
    }

}