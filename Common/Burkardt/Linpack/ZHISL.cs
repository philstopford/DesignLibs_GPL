﻿using System;
using System.Numerics;
using Burkardt.BLAS;

namespace Burkardt.Linpack
{
    public static class ZHISL
    {
        public static void zhisl(Complex[] a, int lda, int n, int[] ipvt,
                ref Complex[] b)

            //*****************************************************************************
            //
            //  Purpose:
            //
            //    ZHISL solves a complex hermitian system factored by ZHIFA.
            //
            //  Discussion:
            //
            //    A division by zero may occur if ZHICO has set RCOND == 0.0
            //    or ZHIFA has set INFO != 0.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 May 2006
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
            //    Input, Complex A[LDA*N], the output from ZHIFA.
            //
            //    Input, int LDA, the leading dimension of A.
            //
            //    Input, int N, the order of the matrix.
            //
            //    Input, int IPVT[N], the pivot vector from ZHIFA.
            //
            //    Input/output, Complex B[N].  On input, the right hand side.
            //    On output, the solution.
            //
        {
            Complex ak;
            Complex akm1;
            Complex bk;
            Complex bkm1;
            Complex denom;
            int k;
            int kp;
            Complex t;
            //
            //  Loop backward applying the transformations and D inverse to B.
            //
            k = n;

            while (0 < k)
            {
                //
                //  1 x 1 pivot block.
                //
                if (0 <= ipvt[k - 1])
                {
                    if (k != 1)
                    {
                        kp = ipvt[k - 1];

                        if (kp != k)
                        {
                            t = b[k - 1];
                            b[k - 1] = b[kp - 1];
                            b[kp - 1] = t;
                        }

                        BLAS1Z.zaxpy(k - 1, b[k - 1], a, 1, ref b, 1, xIndex: +0 + (k - 1) * lda);
                    }

                    //
                    //  Apply D inverse.
                    //
                    b[k - 1] = b[k - 1] / a[k - 1 + (k - 1) * lda];
                    k = k - 1;
                }
                //
                //  2 x 2 pivot block.
                //
                else
                {
                    if (k != 2)
                    {
                        kp = Math.Abs(ipvt[k - 1]);

                        if (kp != k - 1)
                        {
                            t = b[k - 2];
                            b[k - 2] = b[kp - 1];
                            b[kp - 1] = t;
                        }

                        BLAS1Z.zaxpy(k - 2, b[k - 1], a, 1, ref b, 1, xIndex: +0 + (k - 1) * lda);
                        BLAS1Z.zaxpy(k - 2, b[k - 2], a, 1, ref b, 1, xIndex: +0 + (k - 2) * lda);
                    }

                    //
                    //  Apply D inverse.
                    //
                    ak = a[k - 1 + (k - 1) * lda] / Complex.Conjugate(a[k - 2 + (k - 1) * lda]);
                    akm1 = a[k - 2 + (k - 2) * lda] / a[k - 2 + (k - 1) * lda];
                    bk = b[k - 1] / Complex.Conjugate(a[k - 2 + (k - 1) * lda]);
                    bkm1 = b[k - 2] / a[k - 2 + (k - 1) * lda];
                    denom = ak * akm1 - new Complex(1.0, 0.0);
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
                //
                //  1 x 1 pivot block.
                //
                if (0 <= ipvt[k - 1])
                {
                    if (k != 1)
                    {
                        b[k - 1] = b[k - 1] + BLAS1Z.zdotc(k - 1, a, 1, b, 1, xIndex: +0 + (k - 1) * lda);
                        kp = ipvt[k - 1];

                        if (kp != k)
                        {
                            t = b[k - 1];
                            b[k - 1] = b[kp - 1];
                            b[kp - 1] = t;
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
                        b[k - 1] = b[k - 1] + BLAS1Z.zdotc(k - 1, a, 1, b, 1, xIndex: +0 + (k - 1) * lda);
                        b[k] = b[k] + BLAS1Z.zdotc(k - 1, a, 1, b, 1, xIndex: +0 + k * lda);
                        kp = Math.Abs(ipvt[k - 1]);

                        if (kp != k)
                        {
                            t = b[k - 1];
                            b[k - 1] = b[kp - 1];
                            b[kp - 1] = t;
                        }
                    }

                    k = k + 2;
                }
            }
        }

    }
}