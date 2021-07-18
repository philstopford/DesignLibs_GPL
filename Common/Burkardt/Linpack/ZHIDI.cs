using System;
using System.Numerics;
using Burkardt.BLAS;

namespace Burkardt.Linpack
{
    public static class ZHIDI
    {
        public static void zhidi(ref Complex[] a, int lda, int n, int[] ipvt, ref double[] det,
                ref int[] inert, int job)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ZHIDI computes the determinant and inverse of a matrix factored by ZHIFA.
            //
            //  Discussion:
            //
            //    ZHIDI computes the determinant, inertia (number of positive, zero,
            //    and negative eigenvalues) and inverse of a complex hermitian matrix
            //    using the factors from ZHIFA.
            //
            //    A division by zero may occur if the inverse is requested
            //    and ZHICO has set RCOND == 0.0 or ZHIFA has set INFO /= 0.
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
            //    Input/output, Complex A[LDA*N]; on input, the factored matrix
            //    from ZHIFA.  On output, if the inverse was requested, A contains
            //    the inverse matrix.  The strict lower triangle of A is never
            //    referenced.
            //
            //    Input, int LDA, the leading dimension of A.
            //
            //    Input, int N, the order of the matrix.
            //
            //    Input, int IPVT[N], the pivot vector from ZHIFA.
            //
            //    Output, double DET[2], the determinant of the original matrix.
            //    Determinant = det[0] * 10.0**det[1] with 1.0 <= Math.Abs ( det[0] ) < 10.0
            //    or det[0] = 0.0.
            //
            //    Output, int INERT[3], the inertia of the original matrix.
            //    INERT(1) = number of positive eigenvalues.
            //    INERT(2) = number of negative eigenvalues.
            //    INERT(3) = number of zero eigenvalues.
            //
            //    Input, int JOB, has the decimal expansion ABC where:
            //    if C /= 0, the inverse is computed,
            //    if B /= 0, the determinant is computed,
            //    if A /= 0, the inertia is computed.
            //    For example, JOB = 111 gives all three.
            //
        {
            double ak;
            Complex akkp1;
            double akp1;
            double d;
            int i;
            int j;
            int k;
            int km1;
            int ks;
            int kstep;
            bool nodet;
            bool noert;
            bool noinv;
            double t;
            Complex t2;
            Complex[] work;

            noinv = (job % 10) == 0;
            nodet = (job % 100) / 10 == 0;
            noert = (job % 1000) / 100 == 0;

            if (!nodet || !noert)
            {
                if (!noert)
                {
                    for (i = 0; i < 3; i++)
                    {
                        inert[i] = 0;
                    }
                }

                if (!nodet)
                {
                    det[0] = 1.0;
                    det[1] = 0.0;
                }

                t = 0.0;

                for (k = 0; k < n; k++)
                {
                    d = (a[k + k * lda].Real);
                    //
                    //  Check if 1 by 1.
                    //
                    if (ipvt[k] <= 0)
                    {
                        //
                        //  2 by 2 block
                        //  Use DET = ( D / T * C - T ) * T, T = Math.Abs ( S )
                        //  to avoid underflow/overflow troubles.
                        //  Take two passes through scaling.  Use T for flag.
                        //
                        if (t == 0.0)
                        {
                            t = Complex.Abs(a[k + (k + 1) * lda]);
                            d = (d / t) * (a[k + 1 + (k + 1) * lda].Real) - t;
                        }
                        else
                        {
                            d = t;
                            t = 0.0;
                        }
                    }

                    if (!noert)
                    {
                        if (0.0 < d)
                        {
                            inert[0] = inert[0] + 1;
                        }
                        else if (d < 0.0)
                        {
                            inert[1] = inert[1] + 1;
                        }
                        else if (d == 0.0)
                        {
                            inert[2] = inert[2] + 1;
                        }
                    }

                    if (!nodet)
                    {
                        det[0] = det[0] * d;

                        if (det[0] != 0.0)
                        {
                            while (Math.Abs(det[0]) < 1.0)
                            {
                                det[0] = det[0] * 10.0;
                                det[1] = det[1] - 1.0;
                            }

                            while (10.0 <= Math.Abs(det[0]))
                            {
                                det[0] = det[0] / 10.0;
                                det[1] = det[1] + 1.0;
                            }
                        }
                    }
                }
            }

            //
            //  Compute inverse(A).
            //
            if (!noinv)
            {
                work = new Complex [n];

                k = 1;

                while (k <= n)
                {
                    km1 = k - 1;

                    if (0 <= ipvt[k - 1])
                    {
                        //
                        //  1 by 1
                        //
                        a[k - 1 + (k - 1) * lda] =
                            new Complex(1.0 / (a[k - 1 + (k - 1) * lda].Real), 0.0);

                        if (1 <= km1)
                        {
                            for (i = 1; i <= km1; i++)
                            {
                                work[i - 1] = a[i - 1 + (k - 1) * lda];
                            }

                            for (j = 1; j <= km1; j++)
                            {
                                a[j - 1 + (k - 1) * lda] = BLAS1Z.zdotc(j, a, 1, work, 1, xIndex: +0 + (j - 1) * lda);
                                BLAS1Z.zaxpy(j - 1, work[j - 1], a, 1, ref a, 1, xIndex: +0 + (j - 1) * lda,
                                    yIndex: +0 + (k - 1) * lda);
                            }

                            a[k - 1 + (k - 1) * lda] = a[k - 1 + (k - 1) * lda] + new Complex(
                                (BLAS1Z.zdotc(km1, work, 1, a, 1, yIndex: +0 + (k - 1) * lda).Real), 0.0);
                        }

                        kstep = 1;
                    }
                    else
                    {
                        //
                        //  2 by 2
                        //
                        t = Complex.Abs(a[k - 1 + k * lda]);
                        ak = (a[k - 1 + (k - 1) * lda].Real) / t;
                        akp1 = (a[k + k * lda].Real) / t;
                        akkp1 = a[k - 1 + k * lda] / t;
                        d = t * (ak * akp1 - 1.0);
                        a[k - 1 + (k - 1) * lda] = new Complex(akp1 / d, 0.0);
                        a[k + k * lda] = new Complex(ak / d, 0.0);
                        a[k - 1 + k * lda] = -akkp1 / d;

                        if (1 <= km1)
                        {
                            for (i = 1; i <= km1; i++)
                            {
                                work[i - 1] = a[i - 1 + k * lda];
                            }

                            for (j = 1; j <= km1; j++)
                            {
                                a[j - 1 + k * lda] = BLAS1Z.zdotc(j, a, 1, work, 1, xIndex: +0 + (j - 1) * lda);
                                BLAS1Z.zaxpy(j - 1, work[j - 1], a, 1, ref a, 1, xIndex: +0 + (j - 1) * lda,
                                    yIndex: +0 + k * lda);
                            }

                            a[k + k * lda] = a[k + k * lda] + new Complex(
                                (BLAS1Z.zdotc(km1, work, 1, a, 1, yIndex: +0 + k * lda).Real), 0.0);

                            a[k - 1 + k * lda] = a[k - 1 + k * lda]
                                                 + BLAS1Z.zdotc(km1, a, 1, a, 1, xIndex: +0 + (k - 1) * lda,
                                                     yIndex: +0 + k * lda);

                            for (i = 1; i <= km1; i++)
                            {
                                work[i - 1] = a[i - 1 + (k - 1) * lda];
                            }

                            for (j = 1; j <= km1; j++)
                            {
                                a[j - 1 + (k - 1) * lda] = BLAS1Z.zdotc(j, a, 1, work, 1, xIndex: +0 + (j - 1) * lda);
                                BLAS1Z.zaxpy(j - 1, work[j - 1], a, 1, ref a, 1, xIndex: +0 + (j - 1) * lda,
                                    yIndex: +0 + (k - 1) * lda);
                            }

                            a[k - 1 + (k - 1) * lda] = a[k - 1 + (k - 1) * lda] + new Complex(
                                (BLAS1Z.zdotc(km1, work, 1, a, 1, yIndex: +0 + (k - 1) * lda).Real), 0.0);
                        }

                        kstep = 2;
                    }

                    //
                    //  Swap
                    //
                    ks = Math.Abs(ipvt[k - 1]);

                    if (ks != k)
                    {
                        BLAS1Z.zswap(ks, ref a, 1, ref a, 1, xIndex: +0 + (ks - 1) * lda, yIndex: +0 + (k - 1) * lda);

                        for (j = k; ks <= j; j--)
                        {
                            t2 = Complex.Conjugate(a[j - 1 + (k - 1) * lda]);
                            a[j - 1 + (k - 1) * lda] = Complex.Conjugate(a[ks - 1 + (j - 1) * lda]);
                            a[ks - 1 + (j - 1) * lda] = t2;
                        }

                        if (kstep != 1)
                        {
                            t2 = a[ks - 1 + k * lda];
                            a[ks - 1 + k * lda] = a[k - 1 + k * lda];
                            a[k - 1 + k * lda] = t2;
                        }
                    }

                    k = k + kstep;
                }
            }
        }

    }
}