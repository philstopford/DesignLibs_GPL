using System;
using Burkardt.BLAS;

namespace Burkardt.Linpack;

public static class DSIFA
{
    public static int dsifa(ref double[] a, int lda, int n, ref int[] kpvt)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DSIFA factors a real symmetric matrix.
        //
        //  Discussion:
        //
        //    To solve A*X = B, follow DSIFA by DSISL.
        //
        //    To compute inverse(A)*C, follow DSIFA by DSISL.
        //
        //    To compute determinant(A), follow DSIFA by DSIDI.
        //
        //    To compute inertia(A), follow DSIFA by DSIDI.
        //
        //    To compute inverse(A), follow DSIFA by DSIDI.
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
        //    Input/output, double A[LDA*N].  On input, the symmetric matrix
        //    to be factored.  Only the diagonal and upper triangle are used.
        //    On output, a block diagonal matrix and the multipliers which
        //    were used to obtain it.  The factorization can be written A = U*D*U'
        //    where U is a product of permutation and unit upper triangular
        //    matrices, U' is the transpose of U, and D is block diagonal
        //    with 1 by 1 and 2 by 2 blocks.
        //
        //    Input, int LDA, the leading dimension of the array A.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Output, int KPVT[N], the pivot indices.
        //
        //    Output, integer DSIFA, error flag.
        //    0, normal value.
        //    K, if the K-th pivot block is singular.  This is not an error
        //    condition for this subroutine, but it does indicate that DSISL
        //    or DSIDI may divide by zero if called.
        //
    {
        //
        //  ALPHA is used in choosing pivot block size.
        //
        double alpha = (1.0 + Math.Sqrt(17.0)) / 8.0;

        int info = 0;
        //
        //  Main loop on K, which goes from N to 1.
        //
        int k = n;

        while (0 < k)
        {
            switch (k)
            {
                case 1:
                {
                    kpvt[0] = 1;
                    info = a[0 + 0 * lda] switch
                    {
                        0.0 => 1,
                        _ => info
                    };

                    return info;
                }
            }

            //
            //  This section of code determines the kind of
            //  elimination to be performed.  When it is completed,
            //  KSTEP will be set to the size of the pivot block, and
            //  SWAP will be set to .true. if an interchange is required.
            //
            double absakk = Math.Abs(a[k - 1 + (k - 1) * lda]);
            //
            //  Determine the largest off-diagonal element in column K.
            //
            int imax = BLAS1D.idamax(k - 1, a, 1, index: +0 + (k - 1) * lda);
            double colmax = Math.Abs(a[imax - 1 + (k - 1) * lda]);

            int j;
            int kstep;
            bool swap;
            if (alpha * colmax <= absakk)
            {
                kstep = 1;
                swap = false;
            }
            //
            //  Determine the largest off-diagonal element in row IMAX.
            //
            else
            {
                double rowmax = 0.0;
                int imaxp1 = imax + 1;
                for (j = imaxp1; j <= k; j++)
                {
                    rowmax = Math.Max(rowmax, Math.Abs(a[imax - 1 + (j - 1) * lda]));
                }

                if (imax != 1)
                {
                    int jmax = BLAS1D.idamax(imax - 1, a, 1, index: +0 + (imax - 1) * lda);
                    rowmax = Math.Max(rowmax, Math.Abs(a[jmax - 1 + (imax - 1) * lda]));
                }

                if (alpha * rowmax <= Math.Abs(a[imax - 1 + (imax - 1) * lda]))
                {
                    kstep = 1;
                    swap = true;
                }
                else if (alpha * colmax * (colmax / rowmax) <= absakk)
                {
                    kstep = 1;
                    swap = false;
                }
                else
                {
                    kstep = 2;
                    swap = imax != k - 1;
                }
            }

            switch (Math.Max(absakk, colmax))
            {
                //
                //  Column K is zero.
                //  Set INFO and iterate the loop.
                //
                case 0.0:
                    kpvt[k - 1] = k;
                    info = k;
                    break;
                //
                default:
                {
                    int jj;
                    double mulk;
                    double t;
                    if (kstep != 2)
                    {
                        switch (swap)
                        {
                            case true:
                            {
                                BLAS1D.dswap(imax, ref a, 1, ref a, 1, xIndex: +0 + (imax - 1) * lda,
                                    yIndex: +0 + (k - 1) * lda);

                                for (jj = imax; jj <= k; jj++)
                                {
                                    j = k + imax - jj;
                                    t = a[j - 1 + (k - 1) * lda];
                                    a[j - 1 + (k - 1) * lda] = a[imax - 1 + (j - 1) * lda];
                                    a[imax - 1 + (j - 1) * lda] = t;
                                }

                                break;
                            }
                        }

                        //
                        //  Perform the elimination.
                        //
                        for (jj = 1; jj <= k - 1; jj++)
                        {
                            j = k - jj;
                            mulk = -a[j - 1 + (k - 1) * lda] / a[k - 1 + (k - 1) * lda];
                            t = mulk;
                            BLAS1D.daxpy(j, t, a, 1, ref a, 1, xIndex: +0 + (k - 1) * lda, yIndex: +0 + (j - 1) * lda);
                            a[j - 1 + (k - 1) * lda] = mulk;
                        }

                        kpvt[k - 1] = swap switch
                        {
                            //
                            //  Set the pivot array.
                            //
                            true => imax,
                            _ => k
                        };
                    }
                    //
                    //  2 x 2 pivot block.
                    //
                    //  Perform an interchange.
                    //
                    else
                    {
                        switch (swap)
                        {
                            case true:
                            {
                                BLAS1D.dswap(imax, ref a, 1, ref a, 1, xIndex: +0 + (imax - 1) * lda,
                                    yIndex: +0 + (k - 2) * lda);

                                for (jj = imax; jj <= k - 1; jj++)
                                {
                                    j = k - 1 + imax - jj;
                                    t = a[j - 1 + (k - 1) * lda];
                                    a[j - 1 + (k - 1) * lda] = a[imax - 1 + (j - 1) * lda];
                                    a[imax - 1 + (j - 1) * lda] = t;
                                }

                                t = a[k - 2 + (k - 1) * lda];
                                a[k - 2 + (k - 1) * lda] = a[imax - 1 + (k - 1) * lda];
                                a[imax - 1 + (k - 1) * lda] = t;
                                break;
                            }
                        }

                        //
                        //  Perform the elimination.
                        //
                        if (k - 2 != 0)
                        {
                            double ak = a[k - 1 + (k - 1) * lda] / a[k - 2 + (k - 1) * lda];
                            double akm1 = a[k - 2 + (k - 2) * lda] / a[k - 2 + (k - 1) * lda];
                            double denom = 1.0 - ak * akm1;

                            for (jj = 1; jj <= k - 2; jj++)
                            {
                                j = k - 1 - jj;
                                double bk = a[j - 1 + (k - 1) * lda] / a[k - 2 + (k - 1) * lda];
                                double bkm1 = a[j - 1 + (k - 2) * lda] / a[k - 2 + (k - 1) * lda];
                                mulk = (akm1 * bk - bkm1) / denom;
                                double mulkm1 = (ak * bkm1 - bk) / denom;
                                t = mulk;
                                BLAS1D.daxpy(j, t, a, 1, ref a, 1, xIndex: +0 + (k - 1) * lda, yIndex: +0 + (j - 1) * lda);
                                t = mulkm1;
                                BLAS1D.daxpy(j, t, a, 1, ref a, 1, xIndex: +0 + (k - 2) * lda, yIndex: +0 + (j - 1) * lda);
                                a[j - 1 + (k - 1) * lda] = mulk;
                                a[j - 1 + (k - 2) * lda] = mulkm1;
                            }
                        }

                        kpvt[k - 1] = swap switch
                        {
                            //
                            //  Set the pivot array.
                            //
                            true => -imax,
                            _ => 1 - k
                        };

                        kpvt[k - 2] = kpvt[k - 1];
                    }

                    break;
                }
            }

            k -= kstep;
        }

        return info;
    }

}