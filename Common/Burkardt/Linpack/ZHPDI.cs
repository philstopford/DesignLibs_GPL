using System;
using System.Numerics;
using Burkardt.BLAS;

namespace Burkardt.Linpack
{
    public static class ZHPDI
    {
        public static void zhpdi(ref Complex[] ap, int n, int[] ipvt, ref double[] det,
                ref int[] inert, int job)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ZHPDI: determinant, inertia and inverse of a complex hermitian matrix.
            //
            //  Discussion:
            //
            //    The routine uses the factors from ZHPFA.
            //
            //    The matrix is stored in packed form.
            //
            //    A division by zero will occur if the inverse is requested and ZHPCO has
            //    set RCOND == 0.0 or ZHPFA has set INFO != 0.
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
            //    Input/output, Complex AP[N*(N+1)/2]; on input, the factored matrix
            //    from ZHPFA.  If the inverse was requested, then on output, AP contains
            //    the upper triangle of the inverse of the original matrix, stored in packed
            //    form.  The columns of the upper triangle are stored sequentially in a
            //    one-dimensional array.
            //
            //    Input, int N, the order of the matrix.
            //
            //    Input, int IPVT[N], the pivot vector from ZHPFA.
            //
            //    Output, double DET[2], if requested, the determinant of the original
            //    matrix.  Determinant = DET(1) * 10.0**DET(2) with
            //    1.0 <= abs ( DET(1) ) < 10.0 or DET(1) = 0.0.
            //
            //    Output, int INERT[3], if requested, the inertia of the original matrix.
            //    INERT(1) = number of positive eigenvalues.
            //    INERT(2) = number of negative eigenvalues.
            //    INERT(3) = number of zero eigenvalues.
            //
            //    Input, int JOB, has the decimal expansion ABC where:
            //    if C != 0, the inverse is computed,
            //    if B != 0, the determinant is computed,
            //    if A != 0, the inertia is computed.
            //    For example, JOB = 111 gives all three.
            //
        {
            double ak;
            Complex akkp1;
            double akp1;
            double d;
            int ij;
            int ik;
            int ikp1;
            int iks;
            int j;
            int jb;
            int jk;
            int jkp1;
            int k;
            int kk;
            int kkp1;
            int km1;
            int ks;
            int ksj;
            int kskp1;
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
                    inert[0] = 0;
                    inert[1] = 0;
                    inert[2] = 0;
                }

                if (!nodet)
                {
                    det[0] = 1.0;
                    det[1] = 0.0;
                }

                t = 0.0;
                ik = 0;

                for (k = 1; k <= n; k++)
                {
                    kk = ik + k;
                    d = (ap[kk - 1].Real);
                    //
                    //  Check if 1 by 1
                    //
                    if (ipvt[k - 1] <= 0)
                    {
                        //
                        //  2 by 2 block
                        //  Use DET (D  S; S  C)  =  ( D / T * C - T ) * T, T = abs ( S )
                        //  to avoid underflow/overflow troubles.
                        //  Take two passes through scaling.  Use T for flag.
                        //
                        if (t == 0.0)
                        {
                            ikp1 = ik + k;
                            kkp1 = ikp1 + k;
                            t = Complex.Abs(ap[kkp1 - 1]);
                            d = (d / t) * (ap[kkp1].Real) - t;
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

                    ik = ik + k;
                }
            }

            //
            //  Compute inverse(A).
            //
            if (!noinv)
            {
                work = new Complex [n];

                k = 1;
                ik = 0;

                while (k <= n)
                {
                    km1 = k - 1;
                    kk = ik + k;
                    ikp1 = ik + k;
                    kkp1 = ikp1 + k;
                    //
                    //  1 by 1
                    //
                    if (0 <= ipvt[k - 1])
                    {
                        ap[kk - 1] = new Complex(1.0 / (ap[kk - 1].Real), 0.0);

                        if (1 <= km1)
                        {
                            for (j = 1; j <= km1; j++)
                            {
                                work[j - 1] = ap[ik + j - 1];
                            }

                            ij = 0;
                            for (j = 1; j <= km1; j++)
                            {
                                jk = ik + j;
                                ap[jk - 1] = BLAS1Z.zdotc(j, ap, 1, work, 1, xIndex: +ij);
                                BLAS1Z.zaxpy(j - 1, work[j - 1], ap, 1, ref ap, 1, xIndex: +ij, yIndex: +ik);
                                ij = ij + j;
                            }

                            ap[kk - 1] = ap[kk - 1] + new Complex
                                ((BLAS1Z.zdotc(km1, work, 1, ap, 1, yIndex: +ik).Real), 0.0);
                        }

                        kstep = 1;
                    }
                    //
                    //  2 by 2
                    //
                    else
                    {
                        t = Complex.Abs(ap[kkp1 - 1]);
                        ak = (ap[kk - 1].Real) / t;
                        akp1 = (ap[kkp1].Real) / t;
                        akkp1 = ap[kkp1 - 1] / t;
                        d = t * (ak * akp1 - 1.0);
                        ap[kk - 1] = new Complex(akp1 / d, 0.0);
                        ap[kkp1] = new Complex(ak / d, 0.0);
                        ap[kkp1 - 1] = -akkp1 / d;

                        if (1 <= km1)
                        {
                            for (j = 1; j <= km1; j++)
                            {
                                work[j - 1] = ap[ikp1 + j - 1];
                            }

                            ij = 0;
                            for (j = 1; j <= km1; j++)
                            {
                                jkp1 = ikp1 + j;
                                ap[jkp1 - 1] = BLAS1Z.zdotc(j, ap, 1, work, 1, xIndex: +ij);
                                BLAS1Z.zaxpy(j - 1, work[j - 1], ap, 1, ref ap, 1, xIndex: +ij, yIndex: +ikp1);
                                ij = ij + j;
                            }

                            ap[kkp1] = ap[kkp1] + new Complex
                                ((BLAS1Z.zdotc(km1, work, 1, ap, 1, xIndex: +ikp1).Real), 0.0);

                            ap[kkp1 - 1] = ap[kkp1 - 1] + BLAS1Z.zdotc(km1, ap, 1, ap, 1, xIndex: +ik, yIndex: +ikp1);
                            for (j = 1; j <= km1; j++)
                            {
                                work[j - 1] = ap[ik + j - 1];
                            }

                            ij = 0;

                            for (j = 1; j <= km1; j++)
                            {
                                jk = ik + j;
                                ap[jk - 1] = BLAS1Z.zdotc(j, ap, 1, work, 1, xIndex: +ij);
                                BLAS1Z.zaxpy(j - 1, work[j - 1], ap, 1, ref ap, 1, xIndex: +ij, yIndex: +ik);
                                ij = ij + j;
                            }

                            ap[kk - 1] = ap[kk - 1] + new Complex
                                ((BLAS1Z.zdotc(km1, work, 1, ap, 1, yIndex: +ik).Real), 0.0);
                        }

                        kstep = 2;
                    }

                    //
                    //  Swap
                    //
                    ks = Math.Abs(ipvt[k - 1]);

                    if (ks != k)
                    {
                        iks = (ks * (ks - 1)) / 2;

                        BLAS1Z.zswap(ks, ref ap, 1, ref ap, 1, xIndex: +iks, yIndex: +ik);
                        ksj = ik + ks;

                        for (jb = ks; jb <= k; jb++)
                        {
                            j = k + ks - jb;
                            jk = ik + j;

                            t2 = Complex.Conjugate(ap[jk - 1]);
                            ap[jk - 1] = Complex.Conjugate(ap[ksj - 1]);
                            ap[ksj - 1] = t2;

                            ksj = ksj - (j - 1);
                        }

                        if (kstep != 1)
                        {
                            kskp1 = ikp1 + ks;

                            t2 = ap[kskp1 - 1];
                            ap[kskp1 - 1] = ap[kkp1 - 1];
                            ap[kkp1 - 1] = t2;
                        }
                    }

                    ik = ik + k;

                    if (kstep == 2)
                    {
                        ik = ik + k + 1;
                    }

                    k = k + kstep;
                }
            }
        }

    }
}