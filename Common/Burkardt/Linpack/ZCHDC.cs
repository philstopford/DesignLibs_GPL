using System;
using System.Numerics;
using Burkardt.BLAS;

namespace Burkardt.Linpack
{
    public static class ZCHDC
    {
        public static int zchdc(ref Complex[] a, int lda, int p, ref int[] ipvt, int job)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ZCHDC: Cholesky decomposition of a Hermitian positive definite matrix.
            //
            //  Discussion:
            //
            //    A pivoting option allows the user to estimate the condition of a
            //    Hermitian positive definite matrix or determine the rank of a
            //    Hermitian positive semidefinite matrix.
            //
            //    For Hermitian positive definite matrices, INFO = P is the normal return.
            //
            //    For pivoting with Hermitian positive semidefinite matrices, INFO will
            //    in general be less than P.  However, INFO may be greater than
            //    the rank of A, since rounding error can cause an otherwise zero
            //    element to be positive.  Indefinite systems will always cause
            //    INFO to be less than P.
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
            //    Input/output, complex <double> A[LDA*P].  On input, A contains the matrix
            //    whose decomposition is to be computed.  Only the upper half of A
            //    need be stored.  The lower part of the array A is not referenced.
            //    On output, A contains in its upper half the Cholesky factor
            //    of the matrix A as it has been permuted by pivoting.
            //
            //    Input, int LDA, the leading dimension of A.
            //
            //    Input, int P, the order of the matrix.
            //
            //    Input/output, int IPVT[P].  IPVT is not referenced if JOB == 0.
            //    On input, IPVT contains integers that control the selection of the
            //    pivot elements, if pivoting has been requested.  Each diagonal element
            //    A(K,K) is placed in one of three classes according to the input
            //    value of IPVT(K):
            //      IPVT(K) >  0, X(K) is an initial element.
            //      IPVT(K) == 0, X(K) is a free element.
            //      IPVT(K) <  0, X(K) is a final element.
            //    Before the decomposition is computed, initial elements are moved by
            //    symmetric row and column interchanges to the beginning of the array A
            //    and final elements to the end.  Both initial and final elements
            //    are frozen in place during the computation and only free elements
            //    are moved.  At the K-th stage of the reduction, if A(K,K) is occupied
            //    by a free element, it is interchanged with the largest free element
            //    A(L,L) with K <= L.
            //    On output, IPVT(K) contains the index of the diagonal element
            //    of A that was moved into the J-th position, if pivoting was requested.
            //
            //    Input, int JOB, specifies whether column pivoting is to be done.
            //    0, no pivoting is done.
            //    nonzero, pivoting is done.
            //
            //    Output, int ZCHDC, contains the index of the last positive
            //    diagonal element of the Cholesky factor.
            //
        {
            int i_temp;
            int info;
            int j;
            int k;
            int kb;
            int l;
            double maxdia;
            int maxl;
            bool negk;
            int pl;
            int plp1;
            int pu;
            bool swapk;
            Complex temp;
            Complex[] work;

            pl = 1;
            pu = 0;
            info = p;

            work = new Complex[p];

            if (job != 0)
            {
                //
                //  Pivoting has been requested.  Rearrange the elements according to IPVT.
                //
                for (k = 1; k <= p; k++)
                {
                    swapk = (0 < ipvt[k - 1]);
                    negk = (ipvt[k - 1] < 0);

                    if (negk)
                    {
                        ipvt[k - 1] = -k;
                    }
                    else
                    {
                        ipvt[k - 1] = k;
                    }

                    if (swapk)
                    {
                        if (k != pl)
                        {
                            BLAS1Z.zswap(pl - 1, ref a, 1, ref a, 1, xIndex: +0 + (k - 1) * lda,
                                yIndex: +0 + (pl - 1) * lda);

                            temp = a[k - 1 + (k - 1) * lda];
                            a[k - 1 + (k - 1) * lda] = a[pl - 1 + (pl - 1) * lda];
                            a[pl - 1 + (pl - 1) * lda] = temp;

                            a[pl - 1 + (k - 1) * lda] = Complex.Conjugate(a[pl - 1 + (k - 1) * lda]);
                            plp1 = pl + 1;

                            for (j = plp1; j <= p; j++)
                            {
                                if (j < k)
                                {
                                    temp = Complex.Conjugate(a[pl - 1 + (j - 1) * lda]);
                                    a[pl - 1 + (j - 1) * lda] = Complex.Conjugate(a[j - 1 + (k - 1) * lda]);
                                    a[j - 1 + (k - 1) * lda] = temp;
                                }
                                else if (j != k)
                                {
                                    temp = a[pl - 1 + (j - 1) * lda];
                                    a[pl - 1 + (j - 1) * lda] = a[k - 1 + (j - 1) * lda];
                                    a[k - 1 + (j - 1) * lda] = temp;
                                }
                            }

                            ipvt[k - 1] = ipvt[pl - 1];
                            ipvt[pl - 1] = k;
                        }

                        pl = pl + 1;
                    }
                }

                pu = p;

                for (kb = pl; kb <= p; kb++)
                {
                    k = p - kb + pl;

                    if (ipvt[k - 1] < 0)
                    {
                        ipvt[k - 1] = -ipvt[k - 1];

                        if (pu != k)
                        {
                            BLAS1Z.zswap(k - 1, ref a, 1, ref a, 1, xIndex: +0 + (k - 1) * lda,
                                yIndex: +0 + (pu - 1) * lda);

                            temp = a[k - 1 + (k - 1) * lda];
                            a[k - 1 + (k - 1) * lda] = a[pu - 1 + (pu - 1) * lda];
                            a[pu - 1 + (pu - 1) * lda] = temp;

                            a[k - 1 + (pu - 1) * lda] = Complex.Conjugate(a[k - 1 + (pu - 1) * lda]);

                            for (j = k + 1; j <= p; j++)
                            {
                                if (j < pu)
                                {
                                    temp = Complex.Conjugate(a[k - 1 + (j - 1) * lda]);
                                    a[k - 1 + (j - 1) * lda] = Complex.Conjugate(a[j - 1 + (pu - 1) * lda]);
                                    a[j - 1 + (pu - 1) * lda] = temp;
                                }
                                else if (j != pu)
                                {
                                    temp = a[k - 1 + (j - 1) * lda];
                                    a[k - 1 + (j - 1) * lda] = a[pu - 1 + (j - 1) * lda];
                                    a[pu - 1 + (j - 1) * lda] = temp;
                                }
                            }

                            i_temp = ipvt[k - 1];
                            ipvt[k - 1] = ipvt[pu - 1];
                            ipvt[pu - 1] = i_temp;
                        }

                        pu = pu - 1;
                    }
                }
            }

            for (k = 1; k <= p; k++)
            {
                //
                //  Reduction loop.
                //
                maxdia = a[k - 1 + (k - 1) * lda].Real;
                maxl = k;
                //
                //  Determine the pivot element.
                //
                if (pl <= k && k < pu)
                {
                    for (l = k + 1; l <= pu; l++)
                    {
                        if (maxdia < a[l - 1 + (l - 1) * lda].Real)
                        {
                            maxdia = a[l - 1 + (l - 1) * lda].Real;
                            maxl = l;
                        }
                    }
                }

                //
                //  Quit if the pivot element is not positive.
                //
                if (maxdia <= 0.0)
                {
                    info = k - 1;
                    return info;
                }

                //
                //  Start the pivoting and update IPVT.
                //
                if (k != maxl)
                {
                    BLAS1Z.zswap(k - 1, ref a, 1, ref a, 1, xIndex: +0 + (k - 1) * lda, yIndex: +0 + (maxl - 1) * lda);
                    a[maxl - 1 + (maxl - 1) * lda] = a[k - 1 + (k - 1) * lda];
                    a[k - 1 + (k - 1) * lda] = new Complex(maxdia, 0.0);

                    i_temp = ipvt[maxl - 1];
                    ipvt[maxl - 1] = ipvt[k - 1];
                    ipvt[k - 1] = i_temp;

                    a[k - 1 + (maxl - 1) * lda] = Complex.Conjugate(a[k - 1 + (maxl - 1) * lda]);
                }

                //
                //  Reduction step.  Pivoting is contained across the rows.
                //
                work[k - 1] = new Complex(Math.Sqrt(a[k - 1 + (k - 1) * lda].Real), 0.0);
                a[k - 1 + (k - 1) * lda] = work[k - 1];

                for (j = k + 1; j <= p; j++)
                {
                    if (k != maxl)
                    {
                        if (j < maxl)
                        {
                            temp = Complex.Conjugate(a[k - 1 + (j - 1) * lda]);
                            a[k - 1 + (j - 1) * lda] = Complex.Conjugate(a[j - 1 + (maxl - 1) * lda]);
                            a[j - 1 + (maxl - 1) * lda] = temp;
                        }
                        else if (j != maxl)
                        {
                            temp = a[k - 1 + (j - 1) * lda];
                            a[k - 1 + (j - 1) * lda] = a[maxl - 1 + (j - 1) * lda];
                            a[maxl - 1 + (j - 1) * lda] = temp;
                        }
                    }

                    a[k - 1 + (j - 1) * lda] = a[k - 1 + (j - 1) * lda] / work[k - 1];
                    work[j - 1] = Complex.Conjugate(a[k - 1 + (j - 1) * lda]);
                    temp = -a[k - 1 + (j - 1) * lda];
                    BLAS1Z.zaxpy(j - k, temp, work, 1, ref a, 1, xIndex: +k, yIndex: +k + (j - 1) * lda);
                }
            }

            return info;
        }

    }
}