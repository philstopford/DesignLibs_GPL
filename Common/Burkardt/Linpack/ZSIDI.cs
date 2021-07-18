using System;
using System.Numerics;
using Burkardt.BLAS;
using Burkardt.Types;

namespace Burkardt.Linpack
{
    public static class ZSIDI
    {
        public static void zsidi(ref Complex[] a, int lda, int n, int[] ipvt,
                ref Complex[] det, int job)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ZSIDI computes the determinant and inverse of a matrix factored by ZSIFA.
            //
            //  Discussion:
            //
            //    It is assumed the complex symmetric matrix has already been factored
            //    by ZSIFA.
            //
            //    A division by zero may occur if the inverse is requested
            //    and ZSICO set RCOND == 0.0 or ZSIFA set INFO nonzero.
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
            //    Input/output, Complex A[LDA*N]; on input, the output from ZSIFA.
            //    If the inverse was requested, then on output, A contains the upper triangle
            //    of the inverse of the original matrix.  The strict lower triangle
            //    is never referenced.
            //
            //    Input, int LDA, the leading dimension of A.
            //
            //    Input, int N, the order of the matrix.
            //
            //    Input, int IPVT[N], the pivot vector from ZSIFA.
            //
            //    Output, Complex DET[2], if requested, the determinant of the matrix.
            //    Determinant = DET(1) * 10.0**DET(2) with 1.0 <= abs ( DET(1) ) < 10.0
            //    or DET(1) = 0.0.  Also, DET(2) is strictly real.
            //
            //    Input, int JOB, has the decimal expansion AB where
            //    if B != 0, the inverse is computed,
            //    if A != 0, the determinant is computed,
            //    For example, JOB = 11 gives both.
            //
        {
            Complex ak;
            Complex akkp1;
            Complex akp1;
            Complex d;
            int i;
            int j;
            int jb;
            int k;
            int km1;
            int ks;
            int kstep;
            bool nodet;
            bool noinv;
            Complex t;
            Complex[] work;

            noinv = (job % 10) == 0;
            nodet = (job % 100) / 10 == 0;

            if (!nodet)
            {
                det[0] = new Complex(1.0, 0.0);
                det[1] = new Complex(0.0, 0.0);
                t = new Complex(0.0, 0.0);

                for (k = 1; k <= n; k++)
                {
                    d = a[k - 1 + (k - 1) * lda];
                    //
                    //   2 by 2 block.
                    //   Use det ( D  T ) = ( D / T * C - T ) * T
                    //           ( T  C )
                    //   to avoid underflow/overflow troubles.
                    //   Take two passes through scaling.  Use T for flag.
                    //
                    if (ipvt[k - 1] <= 0)
                    {
                        if (typeMethods.zabs1(t) == 0.0)
                        {
                            t = a[k - 1 + k * lda];
                            d = (d / t) * a[k + k * lda] - t;
                        }
                        else
                        {
                            d = t;
                            t = new Complex(0.0, 0.0);
                        }
                    }

                    det[0] = det[0] * d;

                    if (typeMethods.zabs1(det[0]) != 0.0)
                    {
                        while (typeMethods.zabs1(det[0]) < 1.0)
                        {
                            det[0] = det[0] * new Complex(10.0, 0.0);
                            det[1] = det[1] - new Complex(1.0, 0.0);
                        }

                        while (10.0 <= typeMethods.zabs1(det[0]))
                        {
                            det[0] = det[0] / new Complex(10.0, 0.0);
                            det[1] = det[1] + new Complex(1.0, 0.0);
                        }
                    }
                }
            }

            //
            //  Compute inverse ( A ).
            //
            if (!noinv)
            {
                work = new Complex [n];

                k = 1;

                while (k <= n)
                {
                    km1 = k - 1;
                    //
                    //  1 by 1
                    //
                    if (0 <= ipvt[k - 1])
                    {
                        a[k - 1 + (k - 1) * lda] = new Complex(1.0, 0.0) / a[k - 1 + (k - 1) * lda];

                        if (1 <= km1)
                        {
                            for (i = 1; i <= km1; i++)
                            {
                                work[i - 1] = a[i - 1 + (k - 1) * lda];
                            }

                            for (j = 1; j <= km1; j++)
                            {
                                a[j - 1 + (k - 1) * lda] = BLAS1Z.zdotu(j, a, 1, work, 1, xIndex: +0 + (j - 1) * lda);
                                BLAS1Z.zaxpy(j - 1, work[j - 1], a, 1, ref a, 1, xIndex: +0 + (j - 1) * lda,
                                    yIndex: +0 + (k - 1) * lda);
                            }

                            a[k - 1 + (k - 1) * lda] = a[k - 1 + (k - 1) * lda]
                                                       + BLAS1Z.zdotu(km1, work, 1, a, 1, yIndex: +0 + (k - 1) * lda);
                        }

                        kstep = 1;
                    }
                    //
                    //  2 by 2
                    //
                    else
                    {
                        t = a[k - 1 + k * lda];
                        ak = a[k - 1 + (k - 1) * lda] / t;
                        akp1 = a[k + k * lda] / t;
                        akkp1 = a[k - 1 + k * lda] / t;
                        d = t * (ak * akp1 - new Complex(1.0, 0.0));
                        a[k - 1 + (k - 1) * lda] = akp1 / d;
                        a[k + k * lda] = ak / d;
                        a[k - 1 + k * lda] = -akkp1 / d;

                        if (1 <= km1)
                        {
                            for (i = 1; i <= km1; i++)
                            {
                                work[i - 1] = a[i - 1 + k * lda];
                            }

                            for (j = 1; j <= km1; j++)
                            {
                                a[j - 1 + k * lda] = BLAS1Z.zdotu(j, a, 1, work, 1, xIndex: +0 + (j - 1) * lda);
                                BLAS1Z.zaxpy(j - 1, work[j - 1], a, 1, ref a, 1, xIndex: +0 + (j - 1) * lda,
                                    yIndex: +0 + k * lda);
                            }

                            a[k + k * lda] = a[k + k * lda] + BLAS1Z.zdotu(km1, work, 1, a, 1, yIndex: +0 + k * lda);
                            a[k - 1 + k * lda] = a[k - 1 + k * lda]
                                                 + BLAS1Z.zdotu(km1, a, 1, a, 1, xIndex: +0 + (k - 1) * lda,
                                                     yIndex: +0 + k * lda);

                            for (i = 1; i <= km1; i++)
                            {
                                work[i - 1] = a[i - 1 + (k - 1) * lda];
                            }

                            for (j = 1; j <= km1; j++)
                            {
                                a[j - 1 + (k - 1) * lda] = BLAS1Z.zdotu(j, a, 1, work, 1, xIndex: +0 + (j - 1) * lda);
                                BLAS1Z.zaxpy(j - 1, work[j - 1], a, 1, ref a, 1, xIndex: +0 + (j - 1) * lda,
                                    yIndex: +0 + (k - 1) * lda);
                            }

                            a[k - 1 + (k - 1) * lda] = a[k - 1 + (k - 1) * lda]
                                                       + BLAS1Z.zdotu(km1, work, 1, a, 1, yIndex: +0 + (k - 1) * lda);
                        }

                        kstep = 2;
                    }

                    //
                    //  Swap.
                    //
                    ks = Math.Abs(ipvt[k - 1]);

                    if (ks != k)
                    {
                        BLAS1Z.zswap(ks, ref a, 1, ref a, 1, xIndex: +0 + (ks - 1) * lda, yIndex: +0 + (k - 1) * lda);

                        for (jb = ks; jb <= k; jb++)
                        {
                            j = k + ks - jb;

                            t = a[j - 1 + (k - 1) * lda];
                            a[j - 1 + (k - 1) * lda] = a[ks - 1 + (j - 1) * lda];
                            a[ks - 1 + (j - 1) * lda] = t;
                        }

                        if (kstep != 1)
                        {
                            t = a[ks - 1 + k * lda];
                            a[ks - 1 + k * lda] = a[k - 1 + k * lda];
                            a[k - 1 + k * lda] = t;
                        }
                    }

                    k = k + kstep;
                }
            }
        }

    }
}