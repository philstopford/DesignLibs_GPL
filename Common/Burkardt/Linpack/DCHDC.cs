using System;
using Burkardt.BLAS;

namespace Burkardt.Linpack;

public static class DCHDC
{
    public static int dchdc(ref double[] a, int lda, int p, double[] work, ref int[] ipvt, int job )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DCHDC computes the Cholesky decomposition of a positive definite matrix.
        //
        //  Discussion:
        //
        //    A pivoting option allows the user to estimate the condition of a
        //    positive definite matrix or determine the rank of a positive
        //    semidefinite matrix.
        //
        //    For positive definite matrices, INFO = P is the normal return.
        //
        //    For pivoting with positive semidefinite matrices, INFO will
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
        //    23 June 2009
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
        //    Input/output, double A[LDA*P].
        //    On input, A contains the matrix whose decomposition is to
        //    be computed.  Only the upper half of A need be stored.
        //    The lower part of the array a is not referenced.
        //    On output, A contains in its upper half the Cholesky factor
        //    of the input matrix, as it has been permuted by pivoting.
        //
        //    Input, int LDA, the leading dimension of the array A.
        //
        //    Input, int P, the order of the matrix.
        //
        //    Input, double WORK[P] is a work array.
        //
        //    Input/output, int IPVT[P].
        //    On input, IPVT contains integers that control the selection
        //    of the pivot elements, if pivoting has been requested.
        //    Each diagonal element A(K,K) is placed in one of three classes
        //    according to the value of IPVT(K).
        //
        //      > 0, then X(K) is an initial element.
        //      = 0, then X(K) is a free element.
        //      < 0, then X(K) is a final element.
        //
        //    Before the decomposition is computed, initial elements are moved by
        //    symmetric row and column interchanges to the beginning of the array A
        //    and final elements to the end.  Both initial and final elements are
        //    frozen in place during the computation and only free elements are moved.
        //    At the K-th stage of the reduction, if A(K,K) is occupied by a free
        //    element, it is interchanged with the largest free element A(L,L) with
        //    K <= L.  IPVT is not referenced if JOB is 0.
        //
        //    On output, IPVT(J) contains the index of the diagonal element
        //    of A that was moved into the J-th position, if pivoting was requested.
        //
        //    Input, int JOB, initiates column pivoting.
        //    0, no pivoting is done.
        //    nonzero, pivoting is done.
        //
        //    Output, int DCHDC, contains the index of the last positive diagonal
        //    element of the Cholesky factor.
        //
    {
        int j;
        int k;
        double temp;

        int pl = 1;
        int pu = 0;
        int info = p;
        //
        //  Pivoting has been requested.
        //  Rearrange the the elements according to IPVT.
        //
        if (job != 0)
        {
            for (k = 1; k <= p; k++)
            {
                bool swapk = 0 < ipvt[k - 1];
                bool negk = ipvt[k - 1] < 0;

                ipvt[k - 1] = negk switch
                {
                    true => -k,
                    _ => k
                };

                switch (swapk)
                {
                    case true:
                    {
                        if (k != pl)
                        {
                            BLAS1D.dswap(pl - 1, ref a, 1, ref a, 1, xIndex:  + 0 + (k - 1) * lda, yIndex:  + 0 + (pl - 1) * lda);

                            temp = a[k - 1 + (k - 1) * lda];
                            a[k - 1 + (k - 1) * lda] = a[pl - 1 + (pl - 1) * lda];
                            a[pl - 1 + (pl - 1) * lda] = temp;

                            for (j = pl + 1; j <= p; j++)
                            {
                                if (j < k)
                                {
                                    temp = a[pl - 1 + (j - 1) * lda];
                                    a[pl - 1 + (j - 1) * lda] = a[j - 1 + (k - 1) * lda];
                                    a[j - 1 + (k - 1) * lda] = temp;
                                }
                                else if (k < j)
                                {
                                    temp = a[k - 1 + (j - 1) * lda];
                                    a[k - 1 + (j - 1) * lda] = a[pl - 1 + (j - 1) * lda];
                                    a[pl - 1 + (j - 1) * lda] = temp;
                                }
                            }

                            ipvt[k - 1] = ipvt[pl - 1];
                            ipvt[pl - 1] = k;
                        }

                        pl += 1;
                        break;
                    }
                }
            }

            pu = p;

            for (k = p; pl <= k; k--)
            {
                switch (ipvt[k - 1])
                {
                    case < 0:
                    {
                        ipvt[k - 1] = -ipvt[k - 1];

                        if (pu != k)
                        {
                            BLAS1D.dswap(k - 1, ref a, 1, ref a, 1, xIndex:  + 0 + (k - 1) * lda, yIndex: + 0 + (pu - 1) * lda);

                            temp = a[k - 1 + (k - 1) * lda];
                            a[k - 1 + (k - 1) * lda] = a[pu - 1 + (pu - 1) * lda];
                            a[pu - 1 + (pu - 1) * lda] = temp;

                            for (j = k + 1; j <= p; j++)
                            {
                                if (j < pu)
                                {
                                    temp = a[k - 1 + (j - 1) * lda];
                                    a[k - 1 + (j - 1) * lda] = a[j - 1 + (pu - 1) * lda];
                                    a[j - 1 + (pu - 1) * lda] = temp;
                                }
                                else if (pu < j)
                                {
                                    temp = a[k - 1 + (j - 1) * lda];
                                    a[k - 1 + (j - 1) * lda] = a[pu - 1 + (j - 1) * lda];
                                    a[pu - 1 + (j - 1) * lda] = temp;
                                }
                            }

                            (ipvt[k - 1], ipvt[pu - 1]) = (ipvt[pu - 1], ipvt[k - 1]);
                        }

                        pu -= 1;
                        break;
                    }
                }
            }
        }

        for (k = 1; k <= p; k++)
        {
            //
            //  Reduction loop.
            //
            double maxdia = a[k - 1 + (k - 1) * lda];
            int maxl = k;
            //
            //  Determine the pivot element.
            //
            if (pl <= k && k < pu)
            {
                int l;
                for (l = k + 1; l <= pu; l++)
                {
                    if (!(maxdia < a[l - 1 + (l - 1) * lda]))
                    {
                        continue;
                    }

                    maxdia = a[l - 1 + (l - 1) * lda];
                    maxl = l;
                }
            }

            switch (maxdia)
            {
                //
                //  Quit if the pivot element is not positive.
                //
                case <= 0.0:
                    info = k - 1;
                    return info;
            }

            //
            //  Start the pivoting and update IPVT.
            //
            if (k != maxl)
            {
                BLAS1D.dswap(k - 1, ref a, 1, ref a, 1, xIndex:  + 0 + (k - 1) * lda, yIndex:  + 0 + (maxl - 1) * lda);
                a[maxl - 1 + (maxl - 1) * lda] = a[k - 1 + (k - 1) * lda];
                a[k - 1 + (k - 1) * lda] = maxdia;
                (ipvt[maxl - 1], ipvt[k - 1]) = (ipvt[k - 1], ipvt[maxl - 1]);
            }

            //
            //  Reduction step.
            //  Pivoting is contained across the rows.
            //
            work[k - 1] = Math.Sqrt(a[k - 1 + (k - 1) * lda]);
            a[k - 1 + (k - 1) * lda] = work[k - 1];

            for (j = k + 1; j <= p; j++)
            {
                if (k != maxl)
                {
                    if (j < maxl)
                    {
                        temp = a[k - 1 + (j - 1) * lda];
                        a[k - 1 + (j - 1) * lda] = a[j - 1 + (maxl - 1) * lda];
                        a[j - 1 + (maxl - 1) * lda] = temp;
                    }
                    else if (maxl < j)
                    {
                        temp = a[k - 1 + (j - 1) * lda];
                        a[k - 1 + (j - 1) * lda] = a[maxl - 1 + (j - 1) * lda];
                        a[maxl - 1 + (j - 1) * lda] = temp;
                    }
                }

                a[k - 1 + (j - 1) * lda] /= work[k - 1];
                work[j - 1] = a[k - 1 + (j - 1) * lda];
                temp = -a[k - 1 + (j - 1) * lda];
                BLAS1D.daxpy(j - k, temp, work, 1, ref a, 1, xIndex: + k, yIndex:  + k + (j - 1) * lda);
            }
        }

        return info;
    }
}