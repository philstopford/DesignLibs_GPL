using System;
using Burkardt.BLAS;

namespace Burkardt.Linpack
{
    public static class DSISL
    {
        public static void dsisl(ref double[] a, int lda, int n, int[] kpvt, ref double[] b)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DSISL solves a real symmetric system factored by DSIFA.
            //
            //  Discussion:
            //
            //    To compute inverse(A) * C where C is a matrix with P columns
            //
            //      call dsifa ( a, lda, n, kpvt, info )
            //
            //      if ( info == 0 ) then
            //        do j = 1, p
            //          call dsisl ( a, lda, n, kpvt, c(1,j) )
            //        end do
            //      end if
            //
            //    A division by zero may occur if the inverse is requested
            //    and DSICO has set RCOND == 0.0D+00 or DSIFA has set INFO /= 0.
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
            //    Input, double A[LDA*N], the output from DSIFA.
            //
            //    Input, int LDA, the leading dimension of the array A.
            //
            //    Input, int N, the order of the matrix.
            //
            //    Input, int KPVT[N], the pivot vector from DSIFA.
            //
            //    Input/output, double B[N].  On input, the right hand side.
            //    On output, the solution.
            //
        {
            double ak;
            double akm1;
            double bk;
            double bkm1;
            double denom;
            int k;
            int kp;
            double temp;
            //
            //  Loop backward applying the transformations and D inverse to B.
            //
            k = n;

            while (0 < k)
            {
                if (0 <= kpvt[k - 1])
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
                        BLAS1D.daxpy(k - 1, b[k - 1], a, 1, ref b, 1, xIndex: +0 + (k - 1) * lda);
                    }

                    //
                    //  Apply D inverse.
                    //
                    b[k - 1] = b[k - 1] / a[k - 1 + (k - 1) * lda];
                    k = k - 1;
                }
                else
                {
                    //
                    //  2 x 2 pivot block.
                    //
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
                        BLAS1D.daxpy(k - 2, b[k - 1], a, 1, ref b, 1, xIndex: +0 + (k - 1) * lda);
                        BLAS1D.daxpy(k - 2, b[k - 2], a, 1, ref b, 1, xIndex: +0 + (k - 2) * lda);
                    }

                    //
                    //  Apply D inverse.
                    //
                    ak = a[k - 1 + (k - 1) * lda] / a[k - 2 + (k - 1) * lda];
                    akm1 = a[k - 2 + (k - 2) * lda] / a[k - 2 + (k - 1) * lda];
                    bk = b[k - 1] / a[k - 2 + (k - 1) * lda];
                    bkm1 = b[k - 2] / a[k - 2 + (k - 1) * lda];
                    denom = ak * akm1 - 1.0;
                    b[k - 1] = (akm1 * bk - bkm1) / denom;
                    b[k - 2] = (ak * bkm1 - bk) / denom;
                    k = k - 2;
                }
            }

            //
            //  Loop forward applying the transformations.
            //
            k = 1;

            while (k <= n)
            {
                if (0 <= kpvt[k - 1])
                {
                    //
                    //  1 x 1 pivot block.
                    //
                    if (k != 1)
                    {
                        //
                        //  Apply the transformation.
                        //
                        b[k - 1] = b[k - 1] + BLAS1D.ddot(k - 1, a, 1, b, 1, xIndex: +0 + (k - 1) * lda);
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

                    k = k + 1;
                }
                else
                {
                    //
                    //  2 x 2 pivot block.
                    //
                    if (k != 1)
                    {
                        //
                        //  Apply the transformation.
                        //
                        b[k - 1] = b[k - 1] + BLAS1D.ddot(k - 1, a, 1, b, 1, xIndex: +0 + (k - 1) * lda);
                        b[k] = b[k] + BLAS1D.ddot(k - 1, a, 1, b, 1, xIndex: +0 + k * lda);
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

                    k = k + 2;
                }
            }
        }

    }
}