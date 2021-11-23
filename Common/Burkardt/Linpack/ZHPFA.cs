﻿using System;
using System.Numerics;
using Burkardt.BLAS;
using Burkardt.Types;

namespace Burkardt.Linpack;

public static class ZHPFA
{
    public static int zhpfa(ref Complex[] ap, int n, ref int[] ipvt)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZHPFA factors a complex hermitian packed matrix.
        //
        //  Discussion:
        //
        //    To solve A*X = B, follow ZHPFA by ZHPSL.
        //
        //    To compute inverse(A)*C, follow ZHPFA by ZHPSL.
        //
        //    To compute determinant(A), follow ZHPFA by ZHPDI.
        //
        //    To compute inertia(A), follow ZHPFA by ZHPDI.
        //
        //    To compute inverse(A), follow ZHPFA by ZHPDI.
        //
        //  Packed storage:
        //
        //    The following program segment will pack the upper
        //    triangle of a hermitian matrix.
        //
        //      k = 0
        //      do j = 1, n
        //        do i = 1, j
        //          k = k + 1
        //          ap(k) = a(i,j)
        //        }
        //      }
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
        //    Input/output, Complex AP[N*(N+1)/2]; on input, the packed form
        //    of a hermitian matrix.  The columns of the upper triangle are
        //    stored sequentially in a one-dimensional array.  On output, a
        //    block diagonal matrix and the multipliers which were used to
        //    obtain it stored in packed form.  The factorization can be
        //    written A = U*D*hermitian(U) where U is a product of permutation
        //    and unit upper triangular matrices , hermitian(U) is the
        //    Complex.Conjugateugate transpose of U, and D is block diagonal with 1 by 1
        //    and 2 by 2 blocks.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Output, int IPVT[N], the pivot indices.
        //
        //    Output, int ZHPFA.
        //    0, normal value.
        //    K, if the K-th pivot block is singular.  This is not an error condition
        //    for this subroutine, but it does indicate that ZHPSL or ZHPDI may divide
        //    by zero if called.
        //
    {
        int im = 0;
        //
        //  Initialize.
        //
        //  ALPHA is used in choosing pivot block size.
        //
        double alpha = (1.0 + Math.Sqrt(17.0)) / 8.0;

        int info = 0;
        //
        //  Main loop on K, which goes from N to 1.
        //
        int k = n;
        int ik = n * (n - 1) / 2;

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
                ipvt[0] = 1;
                if (typeMethods.zabs1(ap[0]) == 0.0)
                {
                    info = 1;
                }

                break;
            }

            //
            //  This section of code determines the kind of
            //  elimination to be performed.  When it is completed,
            //  KSTEP will be set to the size of the pivot block, and
            //  SWAP will be set to .true. if an interchange is
            //  required.
            //
            int km1 = k - 1;
            int kk = ik + k;
            double absakk = typeMethods.zabs1(ap[kk - 1]);
            //
            //  Determine the largest off-diagonal element in column K.
            //
            int imax = BLAS1Z.izamax(k - 1, ap, 1, index: +ik);
            int imk = ik + imax;
            double colmax = typeMethods.zabs1(ap[imk - 1]);

            int kstep;
            int j;
            bool swap;
            int imj;
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
                im = imax * (imax - 1) / 2;
                imj = im + 2 * imax;

                for (j = imax + 1; j <= k; j++)
                {
                    rowmax = Math.Max(rowmax, typeMethods.zabs1(ap[imj - 1]));
                    imj += j;
                }

                if (imax != 1)
                {
                    int jmax = BLAS1Z.izamax(imax - 1, ap, 1, index: +im);
                    int jmim = jmax + im;
                    rowmax = Math.Max(rowmax, typeMethods.zabs1(ap[jmim - 1]));
                }

                int imim = imax + im;

                if (alpha * rowmax <= typeMethods.zabs1(ap[imim - 1]))
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

            switch (Math.Max(absakk, colmax))
            {
                //
                //  Column K is zero.  Set INFO and iterate the loop.
                //
                case 0.0:
                {
                    ipvt[k - 1] = k;
                    info = k;
                    ik -= k - 1;
                    switch (kstep)
                    {
                        case 2:
                            ik -= k - 2;
                            break;
                    }

                    k -= kstep;
                    continue;
                }
            }

            int jj;
            int jk;
            Complex mulk;
            Complex t;
            int ij;
            int ijj;
            if (kstep != 2)
            {
                switch (swap)
                {
                    //
                    //  1 x 1 pivot block.
                    //
                    case true:
                    {
                        BLAS1Z.zswap(imax, ref ap, 1, ref ap, 1, xIndex: +im, yIndex: +ik);
                        imj = ik + imax;

                        for (jj = imax; jj <= k; jj++)
                        {
                            j = k + imax - jj;
                            jk = ik + j;

                            t = Complex.Conjugate(ap[jk - 1]);
                            ap[jk - 1] = Complex.Conjugate(ap[imj - 1]);
                            ap[imj - 1] = t;

                            imj -= j - 1;
                        }

                        break;
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
                    t = Complex.Conjugate(mulk);
                    BLAS1Z.zaxpy(j, t, ap, 1, ref ap, 1, xIndex: +ik, yIndex: +ij);
                    ijj = ij + j;
                    ap[ijj - 1] = new Complex(ap[ijj - 1].Real, 0.0);
                    ap[jk - 1] = mulk;
                    ij -= j - 1;
                }

                ipvt[k - 1] = swap switch
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
            else
            {
                int km1k = ik + k - 1;
                int ikm1 = ik - (k - 1);

                int jkm1;
                switch (swap)
                {
                    case true:
                    {
                        BLAS1Z.zswap(imax, ref ap, 1, ref ap, 1, xIndex: +im, yIndex: +ikm1);
                        imj = ikm1 + imax;

                        for (jj = imax; jj <= km1; jj++)
                        {
                            j = km1 + imax - jj;
                            jkm1 = ikm1 + j;

                            t = Complex.Conjugate(ap[jkm1 - 1]);
                            ap[jkm1 - 1] = Complex.Conjugate(ap[imj - 1]);
                            ap[imj - 1] = t;

                            imj -= j - 1;
                        }

                        t = ap[km1k - 1];
                        ap[km1k - 1] = ap[imk - 1];
                        ap[imk - 1] = t;
                        break;
                    }
                }

                //
                //  Perform the elimination.
                //
                int km2 = k - 2;

                if (km2 != 0)
                {
                    Complex ak = ap[kk - 1] / ap[km1k - 1];
                    int km1km1 = ikm1 + k - 1;
                    Complex akm1 = ap[km1km1 - 1] / Complex.Conjugate(ap[km1k - 1]);
                    Complex denom = new Complex(1.0, 0.0) - ak * akm1;
                    ij = ik - (k - 1) - (k - 2);

                    for (jj = 1; jj <= km2; jj++)
                    {
                        j = km1 - jj;
                        jk = ik + j;
                        Complex bk = ap[jk - 1] / ap[km1k - 1];
                        jkm1 = ikm1 + j;
                        Complex bkm1 = ap[jkm1 - 1] / Complex.Conjugate(ap[km1k - 1]);
                        mulk = (akm1 * bk - bkm1) / denom;
                        Complex mulkm1 = (ak * bkm1 - bk) / denom;
                        t = Complex.Conjugate(mulk);
                        BLAS1Z.zaxpy(j, t, ap, 1, ref ap, 1, xIndex: +ik, yIndex: +ij);
                        t = Complex.Conjugate(mulkm1);
                        BLAS1Z.zaxpy(j, t, ap, 1, ref ap, 1, xIndex: +ikm1, yIndex: +ij);
                        ap[jk - 1] = mulk;
                        ap[jkm1 - 1] = mulkm1;
                        ijj = ij + j;
                        ap[ijj - 1] = new Complex(ap[ijj - 1].Real, 0.0);
                        ij -= j - 1;
                    }
                }

                ipvt[k - 1] = swap switch
                {
                    //
                    //  Set the pivot array.
                    //
                    true => -imax,
                    _ => 1 - k
                };

                ipvt[k - 2] = ipvt[k - 1];
            }

            ik -= k - 1;
            switch (kstep)
            {
                case 2:
                    ik -= k - 2;
                    break;
            }

            k -= kstep;
        }

        return info;
    }

}