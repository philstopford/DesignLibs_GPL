﻿using System;
using Burkardt.BLAS;

namespace Burkardt.Linpack;

public static class DSPSL
{
    public static void dspsl(double[] ap, int n, int[] kpvt, ref double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DSPSL solves the real symmetric system factored by DSPFA.
        //
        //  Discussion:
        //
        //    To compute inverse(A) * C where C is a matrix with P columns:
        //
        //      call dspfa ( ap, n, kpvt, info )
        //
        //      if ( info /= 0 ) go to ...
        //
        //      do j = 1, p
        //        call dspsl ( ap, n, kpvt, c(1,j) )
        //      end do
        //
        //    A division by zero may occur if DSPCO has set RCOND == 0.0D+00
        //    or DSPFA has set INFO /= 0.
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
        //    Input, double AP[(N*(N+1))/2], the output from DSPFA.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int KPVT[N], the pivot vector from DSPFA.
        //
        //    Input/output, double B[N].  On input, the right hand side.
        //    On output, the solution.
        //
    {
        int kp;
        double temp;
        //
        //  Loop backward applying the transformations and D inverse to B.
        //
        int k = n;
        int ik = n * (n - 1) / 2;

        while (0 < k)
        {
            int kk = ik + k;

            switch (kpvt[k - 1])
            {
                case >= 0:
                {
                    //
                    //  1 x 1 pivot block.
                    //
                    if (k != 1)
                    {
                        kp = kpvt[k - 1];
                        //
                        //  Interchange.
                        //
                        if (kp != k)
                        {
                            temp = b[k - 1];
                            b[k - 1] = b[kp - 1];
                            b[kp - 1] = temp;
                        }

                        //
                        //  Apply the transformation.
                        //
                        BLAS1D.daxpy(k - 1, b[k - 1], ap, 1, ref b, 1, xIndex: +ik);
                    }

                    //
                    //  Apply D inverse.
                    //
                    b[k - 1] /= ap[kk - 1];
                    k -= 1;
                    ik -= k;
                    break;
                }
                default:
                {
                    //
                    //  2 x 2 pivot block.
                    //
                    int ikm1 = ik - (k - 1);

                    if (k != 2)
                    {
                        kp = Math.Abs(kpvt[k - 1]);
                        //
                        //  Interchange.
                        //
                        if (kp != k - 1)
                        {
                            temp = b[k - 2];
                            b[k - 2] = b[kp - 1];
                            b[kp - 1] = temp;
                        }

                        //
                        //  Apply the transformation.
                        //
                        BLAS1D.daxpy(k - 2, b[k - 1], ap, 1, ref b, 1, xIndex: +ik);
                        BLAS1D.daxpy(k - 2, b[k - 2], ap, 1, ref b, 1, xIndex: +ikm1);
                    }

                    //
                    //  Apply D inverse.
                    //
                    int km1k = ik + k - 1;
                    kk = ik + k;
                    double ak = ap[kk - 1] / ap[km1k - 1];
                    int km1km1 = ikm1 + k - 1;
                    double akm1 = ap[km1km1 - 1] / ap[km1k - 1];
                    double bk = b[k - 1] / ap[km1k - 1];
                    double bkm1 = b[k - 2] / ap[km1k - 1];
                    double denom = ak * akm1 - 1.0;
                    b[k - 1] = (akm1 * bk - bkm1) / denom;
                    b[k - 2] = (ak * bkm1 - bk) / denom;
                    k -= 2;
                    ik = ik - (k + 1) - k;
                    break;
                }
            }
        }

        //
        //  Loop forward applying the transformations.
        //
        k = 1;
        ik = 0;

        while (k <= n)
        {
            switch (kpvt[k - 1])
            {
                case >= 0:
                {
                    //
                    //  1 x 1 pivot block.
                    //
                    if (k != 1)
                    {
                        //
                        //  Apply the transformation.
                        //
                        b[k - 1] += BLAS1D.ddot(k - 1, ap, 1, b, 1, xIndex: +ik);
                        kp = kpvt[k - 1];
                        //
                        //  Interchange.
                        //
                        if (kp != k)
                        {
                            temp = b[k - 1];
                            b[k - 1] = b[kp - 1];
                            b[kp - 1] = temp;
                        }
                    }

                    ik += k;
                    k += 1;
                    break;
                }
                default:
                {
                    //
                    //  2 x 2 pivot block.
                    //
                    if (k != 1)
                    {
                        //
                        //  Apply the transformation.
                        //
                        b[k - 1] += BLAS1D.ddot(k - 1, ap, 1, b, 1, xIndex: +ik);
                        int ikp1 = ik + k;
                        b[k] += BLAS1D.ddot(k - 1, ap, 1, b, 1, xIndex: +ikp1);
                        kp = Math.Abs(kpvt[k - 1]);
                        //
                        //  Interchange.
                        //
                        if (kp != k)
                        {
                            temp = b[k - 1];
                            b[k - 1] = b[kp - 1];
                            b[kp - 1] = temp;
                        }
                    }

                    ik = ik + k + k + 1;
                    k += 2;
                    break;
                }
            }
        }
    }

}