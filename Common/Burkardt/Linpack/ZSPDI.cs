using System;
using System.Numerics;
using Burkardt.BLAS;
using Burkardt.Types;

namespace Burkardt.Linpack;

public static class ZSPDI
{
    public static void zspdi(ref Complex[] ap, int n, int[] ipvt, ref Complex[] det,
            int job)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZSPDI sets the determinant and inverse of a complex symmetric packed matrix.
        //
        //  Discussion:
        //
        //    ZSPDI uses the factors from ZSPFA.
        //
        //    The matrix is stored in packed form.
        //
        //    A division by zero will occur if the inverse is requested and ZSPCO has
        //    set RCOND to 0.0 or ZSPFA has set INFO nonzero.
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
        //    Input/output, Complex AP[N*(N+1)/2]; on input, the matrix factors
        //    from ZSPFA.  On output, if the inverse was requested, the upper
        //    triangle of the inverse of the original matrix, stored in packed
        //    form.  The columns of the upper triangle are stored sequentially
        //    in a one-dimensional array.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int IPVT[N], the pivot vector from ZSPFA.
        //
        //    Output, Complex DET[2], the determinant of the original matrix.
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
        int kkp1 = 0;
        int km1;
        int ks;
        int ksj;
        int kskp1;
        int kstep;
        bool nodet;
        bool noinv;
        Complex t;
        Complex[] work;

        noinv = job % 10 == 0;
        nodet = job % 100 / 10 == 0;

        switch (nodet)
        {
            case false:
            {
                det[0] = new Complex(1.0, 0.0);
                det[1] = new Complex(0.0, 0.0);
                t = new Complex(0.0, 0.0);
                ik = 0;

                for (k = 1; k <= n; k++)
                {
                    kk = ik + k;
                    d = ap[kk - 1];
                    switch (ipvt[k - 1])
                    {
                        //
                        //  2 by 2 block
                        //  Use det (D  T)  =  ( D / T * C - T ) * T
                        //          (T  C)
                        //  to avoid underflow/overflow troubles.
                        //  Take two passes through scaling.  Use T for flag.
                        //
                        case <= 0 when typeMethods.zabs1(t) == 0.0:
                            ikp1 = ik + k;
                            kkp1 = ikp1 + k;
                            t = ap[kkp1 - 1];
                            d = d / t * ap[kkp1] - t;
                            break;
                        case <= 0:
                            d = t;
                            t = new Complex(0.0, 0.0);
                            break;
                    }

                    switch (nodet)
                    {
                        case false:
                        {
                            det[0] *= d;

                            if (typeMethods.zabs1(det[0]) != 0.0)
                            {
                                while (typeMethods.zabs1(det[0]) < 1.0)
                                {
                                    det[0] *= new Complex(10.0, 0.0);
                                    det[1] -= new Complex(1.0, 0.0);
                                }

                                while (10.0 <= typeMethods.zabs1(det[0]))
                                {
                                    det[0] /= new Complex(10.0, 0.0);
                                    det[1] += new Complex(1.0, 0.0);
                                }
                            }

                            break;
                        }
                    }

                    ik += k;
                }

                break;
            }
        }

        switch (noinv)
        {
            //
            //  Compute inverse ( A ).
            //
            case false:
            {
                work = new Complex[n];
                k = 1;
                ik = 0;

                while (k <= n)
                {
                    km1 = k - 1;
                    kk = ik + k;
                    ikp1 = ik + k;

                    switch (ipvt[k - 1])
                    {
                        case >= 0:
                        {
                            //
                            //  1 by 1
                            //
                            ap[kk - 1] = new Complex(1.0, 0.0) / ap[kk - 1];

                            switch (km1)
                            {
                                case >= 1:
                                {
                                    for (i = 1; i <= km1; i++)
                                    {
                                        work[i - 1] = ap[ik + i - 1];
                                    }

                                    ij = 0;

                                    for (j = 1; j <= km1; j++)
                                    {
                                        jk = ik + j;
                                        ap[jk - 1] = BLAS1Z.zdotu(j, ap, 1, work, 1, xIndex: +ij);
                                        BLAS1Z.zaxpy(j - 1, work[j - 1], ap, 1, ref ap, 1, xIndex: +ij, yIndex: +ik);
                                        ij += j;
                                    }

                                    ap[kk - 1] += BLAS1Z.zdotu(km1, work, 1, ap, 1, yIndex: +ik);
                                    break;
                                }
                            }

                            kstep = 1;
                            break;
                        }
                        //
                        default:
                        {
                            kkp1 = ikp1 + k;
                            t = ap[kkp1 - 1];
                            ak = ap[kk - 1] / t;
                            akp1 = ap[kkp1] / t;
                            akkp1 = ap[kkp1 - 1] / t;
                            d = t * (ak * akp1 - new Complex(1.0, 0.0));
                            ap[kk - 1] = akp1 / d;
                            ap[kkp1] = ak / d;
                            ap[kkp1 - 1] = -akkp1 / d;

                            switch (km1)
                            {
                                case >= 1:
                                {
                                    for (i = 1; i <= km1; i++)
                                    {
                                        work[i - 1] = ap[ikp1 - 1];
                                    }

                                    ij = 0;

                                    for (j = 1; j <= km1; j++)
                                    {
                                        jkp1 = ikp1 + j;
                                        ap[jkp1 - 1] = BLAS1Z.zdotu(j, ap, 1, work, 1, xIndex: +ij);
                                        BLAS1Z.zaxpy(j - 1, work[j - 1], ap, 1, ref ap, 1, xIndex: +ij, yIndex: +ikp1);
                                        ij += j;
                                    }

                                    ap[kkp1] += BLAS1Z.zdotu(km1, work, 1, ap, 1, yIndex: +ikp1);
                                    ap[kkp1 - 1] += BLAS1Z.zdotu(km1, ap, 1, ap, 1, xIndex: +ik, yIndex: +ikp1);

                                    for (i = 1; i <= km1; i++)
                                    {
                                        work[i - 1] = ap[ik + i - 1];
                                    }

                                    ij = 0;

                                    for (j = 1; j <= km1; j++)
                                    {
                                        jk = ik + j;
                                        ap[jk - 1] = BLAS1Z.zdotu(j, ap, 1, work, 1, xIndex: +ij);
                                        BLAS1Z.zaxpy(j - 1, work[j - 1], ap, 1, ref ap, 1, xIndex: +ij, yIndex: +ik);
                                        ij += j;
                                    }

                                    ap[kk - 1] += BLAS1Z.zdotu(km1, work, 1, ap, 1, yIndex: +ik);
                                    break;
                                }
                            }

                            kstep = 2;
                            break;
                        }
                    }

                    //
                    //  Swap.
                    //
                    ks = Math.Abs(ipvt[k - 1]);

                    if (ks != k)
                    {
                        iks = ks * (ks - 1) / 2;
                        BLAS1Z.zswap(ks, ref ap, 1, ref ap, 1, xIndex: +iks, yIndex: +ik);
                        ksj = ik + ks;

                        for (jb = ks; jb <= k; jb++)
                        {
                            j = k + ks - jb;
                            jk = ik + j;

                            t = ap[jk - 1];
                            ap[jk - 1] = ap[ksj - 1];
                            ap[ksj - 1] = t;

                            ksj -= (j - 1);
                        }

                        if (kstep != 1)
                        {
                            kskp1 = ikp1 + ks;

                            t = ap[kskp1 - 1];
                            ap[kskp1 - 1] = ap[kkp1 - 1];
                            ap[kkp1 - 1] = t;
                        }
                    }

                    ik += k;

                    ik = kstep switch
                    {
                        2 => ik + k + 1,
                        _ => ik
                    };

                    k += kstep;
                }

                break;
            }
        }
    }

}