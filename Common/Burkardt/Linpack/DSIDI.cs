using System;
using Burkardt.BLAS;

namespace Burkardt.Linpack
{
    public static class DSIDI
    {
        public static void dsidi(ref double[] a, int lda, int n, int[] kpvt, ref double[] det,
                ref int[] inert, double[] work, int job)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DSIDI computes the determinant, inertia and inverse of a real symmetric matrix.
            //
            //  Discussion:
            //
            //    DSIDI uses the factors from DSIFA.
            //
            //    A division by zero may occur if the inverse is requested
            //    and DSICO has set RCOND == 0.0D+00 or DSIFA has set INFO /= 0.
            //
            //    Variables not requested by JOB are not used.
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
            //    Input/output, double A(LDA,N).  On input, the output from DSIFA.
            //    On output, the upper triangle of the inverse of the original matrix,
            //    if requested.  The strict lower triangle is never referenced.
            //
            //    Input, int LDA, the leading dimension of the array A.
            //
            //    Input, int N, the order of the matrix.
            //
            //    Input, int KPVT[N], the pivot vector from DSIFA.
            //
            //    Output, double DET[2], the determinant of the original matrix,
            //    if requested.
            //      determinant = DET[0] * 10.0**DET[1]
            //    with 1.0D+00 <= abs ( DET[0] ) < 10.0D+00 or DET[0] = 0.0.
            //
            //    Output, int INERT(3), the inertia of the original matrix,
            //    if requested.
            //    INERT(1) = number of positive eigenvalues.
            //    INERT(2) = number of negative eigenvalues.
            //    INERT(3) = number of zero eigenvalues.
            //
            //    Workspace, double WORK[N].
            //
            //    Input, int JOB, specifies the tasks.
            //    JOB has the decimal expansion ABC where
            //    If C /= 0, the inverse is computed,
            //    If B /= 0, the determinant is computed,
            //    If A /= 0, the inertia is computed.
            //    For example, JOB = 111 gives all three.
            //
        {
            double ak;
            double akkp1;
            double akp1;
            double d;
            bool dodet;
            bool doert;
            bool doinv;
            int j;
            int jb;
            int k;
            int ks;
            int kstep;
            double t;
            double temp;

            doinv = (job % 10) != 0;
            dodet = (job % 100) / 10 != 0;
            doert = (job % 1000) / 100 != 0;

            if (dodet || doert)
            {
                if (doert)
                {
                    inert[0] = 0;
                    inert[1] = 0;
                    inert[2] = 0;
                }

                if (dodet)
                {
                    det[0] = 1.0;
                    det[1] = 0.0;
                }

                t = 0.0;

                for (k = 1; k <= n; k++)
                {
                    d = a[k - 1 + (k - 1) * lda];
                    //
                    //  2 by 2 block.
                    //
                    //  use det (d  s)  =  (d/t * c - t) * t,  t = abs ( s )
                    //          (s  c)
                    //  to avoid underflow/overflow troubles.
                    //
                    //  Take two passes through scaling.  Use T for flag.
                    //
                    if (kpvt[k - 1] <= 0)
                    {
                        if (t == 0.0)
                        {
                            t = Math.Abs(a[k - 1 + k * lda]);
                            d = (d / t) * a[k + k * lda] - t;
                        }
                        else
                        {
                            d = t;
                            t = 0.0;
                        }
                    }

                    if (doert)
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

                    if (dodet)
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
            if (doinv)
            {
                k = 1;

                while (k <= n)
                {
                    if (0 <= kpvt[k - 1])
                    {
                        //
                        //  1 by 1.
                        //
                        a[k - 1 + (k - 1) * lda] = 1.0 / a[k - 1 + (k - 1) * lda];

                        if (2 <= k)
                        {
                            BLAS1D.dcopy(k - 1, a, 1, ref work, 1, xIndex: +0 + (k - 1) * lda);

                            for (j = 1; j <= k - 1; j++)
                            {
                                a[j - 1 + (k - 1) * lda] = BLAS1D.ddot(j, a, 1, work, 1, xIndex: +0 + (j - 1) * lda);
                                BLAS1D.daxpy(j - 1, work[j - 1], a, 1, ref a, 1, xIndex: +0 + (j - 1) * lda,
                                    yIndex: +0 + (k - 1) * lda);
                            }

                            a[k - 1 + (k - 1) * lda] = a[k - 1 + (k - 1) * lda]
                                                       + BLAS1D.ddot(k - 1, work, 1, a, 1, yIndex: +0 + (k - 1) * lda);
                        }

                        kstep = 1;
                    }
                    //
                    //  2 by 2.
                    //
                    else
                    {
                        t = Math.Abs(a[k - 1 + k * lda]);
                        ak = a[k - 1 + (k - 1) * lda] / t;
                        akp1 = a[k + k * lda] / t;
                        akkp1 = a[k - 1 + k * lda] / t;
                        d = t * (ak * akp1 - 1.0);
                        a[k - 1 + (k - 1) * lda] = akp1 / d;
                        a[k + k * lda] = ak / d;
                        a[k - 1 + k * lda] = -akkp1 / d;

                        if (2 <= k)
                        {
                            BLAS1D.dcopy(k - 1, a, 1, ref work, 1, xIndex: +0 + k * lda);

                            for (j = 1; j <= k - 1; j++)
                            {
                                a[j - 1 + k * lda] = BLAS1D.ddot(j, a, 1, work, 1, xIndex: +0 + (j - 1) * lda);
                                BLAS1D.daxpy(j - 1, work[j - 1], a, 1, ref a, 1, xIndex: +0 + (j - 1) * lda,
                                    yIndex: +0 + k * lda);
                            }

                            a[k + k * lda] = a[k + k * lda] + BLAS1D.ddot(k - 1, work, 1, a, 1, yIndex: +0 + k * lda);
                            a[k - 1 + k * lda] = a[k - 1 + k * lda]
                                                 + BLAS1D.ddot(k - 1, a, 1, a, 1, xIndex: +0 + (k - 1) * lda,
                                                     yIndex: +0 + k * lda);
                            BLAS1D.dcopy(k - 1, a, 1, ref work, 1, xIndex: +0 + (k - 1) * lda);

                            for (j = 1; j <= k - 1; j++)
                            {
                                a[j - 1 + (k - 1) * lda] = BLAS1D.ddot(j, a, 1, work, 1, xIndex: +0 + (j - 1) * lda);
                                BLAS1D.daxpy(j - 1, work[j - 1], a, 1, ref a, 1, xIndex: +0 + (j - 1) * lda,
                                    yIndex: +0 + (k - 1) * lda);
                            }

                            a[k - 1 + (k - 1) * lda] = a[k - 1 + (k - 1) * lda]
                                                       + BLAS1D.ddot(k - 1, work, 1, a, 1, yIndex: +0 + (k - 1) * lda);
                        }

                        kstep = 2;
                    }

                    //
                    //  Swap.
                    //
                    ks = Math.Abs(kpvt[k - 1]);

                    if (ks != k)
                    {
                        BLAS1D.dswap(ks, ref a, 1, ref a, 1, xIndex: +0 + (ks - 1) * lda, yIndex: +0 + (k - 1) * lda);

                        for (jb = ks; jb <= k; jb++)
                        {
                            j = k + ks - jb;
                            temp = a[j - 1 + (k - 1) * lda];
                            a[j - 1 + (k - 1) * lda] = a[ks - 1 + (j - 1) * lda];
                            a[ks - 1 + (j - 1) * lda] = temp;
                        }

                        if (kstep != 1)
                        {
                            temp = a[ks - 1 + k * lda];
                            a[ks - 1 + k * lda] = a[k - 1 + k * lda];
                            a[k - 1 + k * lda] = temp;
                        }
                    }

                    k = k + kstep;
                }
            }
        }

    }
}