using System;
using Burkardt.BLAS;

namespace Burkardt.Linpack
{
    public static class DSPFA
    {
        public static int dspfa(ref double[] ap, int n, ref int[] kpvt)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DSPFA factors a real symmetric matrix stored in packed form.
            //
            //  Discussion:
            //
            //    To solve A*X = B, follow DSPFA by DSPSL.
            //
            //    To compute inverse(A)*C, follow DSPFA by DSPSL.
            //
            //    To compute determinant(A), follow DSPFA by DSPDI.
            //
            //    To compute inertia(A), follow DSPFA by DSPDI.
            //
            //    To compute inverse(A), follow DSPFA by DSPDI.
            //
            //  Packed storage:
            //
            //    The following program segment will pack the upper triangle of a
            //    symmetric matrix.
            //
            //      k = 0
            //      do j = 1, n
            //        do i = 1, j
            //          k = k + 1
            //          ap(k) = a(i,j)
            //        end do
            //      end do
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
            //    Input/output, double AP[(N*(N+1))/2].  On input, the packed form of a
            //    symmetric matrix A.  The columns of the upper triangle are stored
            //    sequentially in a one-dimensional array.  On output, a block diagonal
            //    matrix and the multipliers which were used to obtain it stored in
            //    packed form.  The factorization can be written A = U*D*U' where U
            //    is a product of permutation and unit upper triangular matrices, U'
            //    is the transpose of U, and D is block diagonal with 1 by 1 and 2
            //    by 2 blocks.
            //
            //    Input, int N, the order of the matrix.
            //
            //    Output, int KPVT[N], the pivot indices.
            //
            //    Output, int DSPFA, error flag.
            //    0, normal value.
            //    K, if the K-th pivot block is singular.  This is not an error
            //    condition for this subroutine, but it does indicate that DSPSL or
            //    DSPDI may divide by zero if called.
            //
        {
            double absakk;
            double ak;
            double akm1;
            double alpha;
            double bk;
            double bkm1;
            double colmax;
            double denom;
            int ij;
            int ik;
            int ikm1;
            int im = 0;
            int imax;
            int imaxp1;
            int imim;
            int imj;
            int imk;
            int info;
            int j;
            int jj;
            int jk;
            int jkm1;
            int jmax;
            int jmim;
            int k;
            int kk;
            int km1;
            int km1k;
            int km1km1;
            int kstep;
            double mulk;
            double mulkm1;
            double rowmax;
            bool swap;
            double t;
            //
            //  ALPHA is used in choosing pivot block size.
            //
            alpha = (1.0 + Math.Sqrt(17.0)) / 8.0;

            info = 0;
            //
            //  Main loop on K, which goes from N to 1.
            //
            k = n;
            ik = (n * (n - 1)) / 2;

            for (;;)
            {
                //
                //  Leave the loop if K = 0 or K = 1.
                //
                if (k == 0)
                {
                    break;
                }

                if (k == 1)
                {
                    kpvt[0] = 1;
                    if (ap[0] == 0.0)
                    {
                        info = 1;
                    }

                    break;
                }

                //
                //  This section of code determines the kind of elimination to be performed.
                //  When it is completed, KSTEP will be set to the size of the pivot block,
                //  and SWAP will be set to .true. if an interchange is required.
                //
                km1 = k - 1;
                kk = ik + k;
                absakk = Math.Abs(ap[kk - 1]);
                //
                //  Determine the largest off-diagonal element in column K.
                //
                imax = BLAS1D.idamax(k - 1, ap, 1, index: +ik);
                imk = ik + imax;
                colmax = Math.Abs(ap[imk - 1]);

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
                    rowmax = 0.0;
                    imaxp1 = imax + 1;
                    im = (imax * (imax - 1)) / 2;
                    imj = im + 2 * imax;

                    for (j = imaxp1; j <= k; j++)
                    {
                        rowmax = Math.Max(rowmax, Math.Abs(ap[imj - 1]));
                        imj = imj + j;
                    }

                    if (imax != 1)
                    {
                        jmax = BLAS1D.idamax(imax - 1, ap, 1, index: +im);
                        jmim = jmax + im;
                        rowmax = Math.Max(rowmax, Math.Abs(ap[jmim - 1]));
                    }

                    imim = imax + im;

                    if (alpha * rowmax <= Math.Abs(ap[imim - 1]))
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
                        swap = imax != km1;
                    }
                }

                //
                //  Column K is zero.  Set INFO and iterate the loop.
                //
                if (Math.Max(absakk, colmax) == 0.0)
                {
                    kpvt[k - 1] = k;
                    info = k;
                }
                else
                {
                    if (kstep != 2)
                    {
                        //
                        //  1 x 1 pivot block.
                        //
                        if (swap)
                        {
                            //
                            //  Perform an interchange.
                            //
                            BLAS1D.dswap(imax, ref ap, 1, ref ap, 1, xIndex: +im, yIndex: +ik);
                            imj = ik + imax;

                            for (jj = imax; jj <= k; jj++)
                            {
                                j = k + imax - jj;
                                jk = ik + j;
                                t = ap[jk - 1];
                                ap[jk - 1] = ap[imj - 1];
                                ap[imj - 1] = t;
                                imj = imj - (j - 1);
                            }
                        }

                        //
                        //  Perform the elimination.
                        //
                        ij = ik - (k - 1);

                        for (jj = 1; jj <= km1; jj++)
                        {
                            j = k - jj;
                            jk = ik + j;
                            mulk = -ap[jk - 1] / ap[kk - 1];
                            t = mulk;
                            BLAS1D.daxpy(j, t, ap, 1, ref ap, 1, xIndex: +ik, yIndex: +ij);
                            ap[jk - 1] = mulk;
                            ij = ij - (j - 1);
                        }

                        //
                        //  Set the pivot array.
                        //
                        if (swap)
                        {
                            kpvt[k - 1] = imax;
                        }
                        else
                        {
                            kpvt[k - 1] = k;
                        }
                    }
                    else
                    {
                        //
                        //  2 x 2 pivot block.
                        //
                        km1k = ik + k - 1;
                        ikm1 = ik - (k - 1);
                        //
                        //  Perform an interchange.
                        //
                        if (swap)
                        {
                            BLAS1D.dswap(imax, ref ap, 1, ref ap, 1, xIndex: +im, yIndex: +ikm1);
                            imj = ikm1 + imax;

                            for (jj = imax; jj <= km1; jj++)
                            {
                                j = km1 + imax - jj;
                                jkm1 = ikm1 + j;
                                t = ap[jkm1 - 1];
                                ap[jkm1 - 1] = ap[imj - 1];
                                ap[imj - 1] = t;
                                imj = imj - (j - 1);
                            }

                            t = ap[km1k - 1];
                            ap[km1k - 1] = ap[imk - 1];
                            ap[imk - 1] = t;
                        }

                        //
                        //  Perform the elimination.
                        //
                        if (k - 2 != 0)
                        {
                            ak = ap[kk - 1] / ap[km1k - 1];
                            km1km1 = ikm1 + k - 1;
                            akm1 = ap[km1km1 - 1] / ap[km1k - 1];
                            denom = 1.0 - ak * akm1;
                            ij = ik - (k - 1) - (k - 2);

                            for (jj = 1; jj <= k - 2; jj++)
                            {
                                j = km1 - jj;
                                jk = ik + j;
                                bk = ap[jk - 1] / ap[km1k - 1];
                                jkm1 = ikm1 + j;
                                bkm1 = ap[jkm1 - 1] / ap[km1k - 1];
                                mulk = (akm1 * bk - bkm1) / denom;
                                mulkm1 = (ak * bkm1 - bk) / denom;
                                t = mulk;
                                BLAS1D.daxpy(j, t, ap, 1, ref ap, 1, xIndex: +ik, yIndex: +ij);
                                t = mulkm1;
                                BLAS1D.daxpy(j, t, ap, 1, ref ap, 1, xIndex: +ikm1, yIndex: +ij);
                                ap[jk - 1] = mulk;
                                ap[jkm1 - 1] = mulkm1;
                                ij = ij - (j - 1);
                            }
                        }

                        //
                        //  Set the pivot array.
                        //
                        if (swap)
                        {
                            kpvt[k - 1] = -imax;
                        }
                        else
                        {
                            kpvt[k - 1] = 1 - k;
                        }

                        kpvt[k - 2] = kpvt[k - 1];
                    }
                }

                ik = ik - (k - 1);
                if (kstep == 2)
                {
                    ik = ik - (k - 2);
                }

                k = k - kstep;

            }

            return info;
        }

    }
}