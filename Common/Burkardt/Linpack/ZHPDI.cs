using System;
using System.Numerics;
using Burkardt.BLAS;

namespace Burkardt.Linpack;

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
        double d;
        int ik;
        int ikp1;
        int k;
        int kk;
        int kkp1;
        double t;

        bool noinv = job % 10 == 0;
        bool nodet = job % 100 / 10 == 0;
        bool noert = job % 1000 / 100 == 0;

        if (!nodet || !noert)
        {
            switch (noert)
            {
                case false:
                    inert[0] = 0;
                    inert[1] = 0;
                    inert[2] = 0;
                    break;
            }

            switch (nodet)
            {
                case false:
                    det[0] = 1.0;
                    det[1] = 0.0;
                    break;
            }

            t = 0.0;
            ik = 0;

            for (k = 1; k <= n; k++)
            {
                kk = ik + k;
                d = ap[kk - 1].Real;
                switch (ipvt[k - 1])
                {
                    //
                    //  Check if 1 by 1
                    //
                    //
                    //  2 by 2 block
                    //  Use DET (D  S; S  C)  =  ( D / T * C - T ) * T, T = abs ( S )
                    //  to avoid underflow/overflow troubles.
                    //  Take two passes through scaling.  Use T for flag.
                    //
                    case <= 0 when t == 0.0:
                        ikp1 = ik + k;
                        kkp1 = ikp1 + k;
                        t = Complex.Abs(ap[kkp1 - 1]);
                        d = d / t * ap[kkp1].Real - t;
                        break;
                    case <= 0:
                        d = t;
                        t = 0.0;
                        break;
                }

                switch (noert)
                {
                    case false:
                        switch (d)
                        {
                            case > 0.0:
                                inert[0] += 1;
                                break;
                            case < 0.0:
                                inert[1] += 1;
                                break;
                            case 0.0:
                                inert[2] += 1;
                                break;
                        }

                        break;
                }

                switch (nodet)
                {
                    case false:
                    {
                        det[0] *= d;

                        if (det[0] != 0.0)
                        {
                            while (Math.Abs(det[0]) < 1.0)
                            {
                                det[0] *= 10.0;
                                det[1] -= 1.0;
                            }

                            while (10.0 <= Math.Abs(det[0]))
                            {
                                det[0] /= 10.0;
                                det[1] += 1.0;
                            }
                        }

                        break;
                    }
                }

                ik += k;
            }
        }

        switch (noinv)
        {
            //
            //  Compute inverse(A).
            //
            case false:
            {
                Complex[] work = new Complex [n];

                k = 1;
                ik = 0;

                while (k <= n)
                {
                    int km1 = k - 1;
                    kk = ik + k;
                    ikp1 = ik + k;
                    kkp1 = ikp1 + k;
                    int kstep;
                    int jk;
                    int j;
                    int ij;
                    switch (ipvt[k - 1])
                    {
                        //
                        //  1 by 1
                        //
                        case >= 0:
                        {
                            ap[kk - 1] = new Complex(1.0 / ap[kk - 1].Real, 0.0);

                            switch (km1)
                            {
                                case >= 1:
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
                                        ij += j;
                                    }

                                    ap[kk - 1] += new Complex
                                        (BLAS1Z.zdotc(km1, work, 1, ap, 1, yIndex: +ik).Real, 0.0);
                                    break;
                                }
                            }

                            kstep = 1;
                            break;
                        }
                        //
                        default:
                        {
                            t = Complex.Abs(ap[kkp1 - 1]);
                            double ak = ap[kk - 1].Real / t;
                            double akp1 = ap[kkp1].Real / t;
                            Complex akkp1 = ap[kkp1 - 1] / t;
                            d = t * (ak * akp1 - 1.0);
                            ap[kk - 1] = new Complex(akp1 / d, 0.0);
                            ap[kkp1] = new Complex(ak / d, 0.0);
                            ap[kkp1 - 1] = -akkp1 / d;

                            switch (km1)
                            {
                                case >= 1:
                                {
                                    for (j = 1; j <= km1; j++)
                                    {
                                        work[j - 1] = ap[ikp1 + j - 1];
                                    }

                                    ij = 0;
                                    for (j = 1; j <= km1; j++)
                                    {
                                        int jkp1 = ikp1 + j;
                                        ap[jkp1 - 1] = BLAS1Z.zdotc(j, ap, 1, work, 1, xIndex: +ij);
                                        BLAS1Z.zaxpy(j - 1, work[j - 1], ap, 1, ref ap, 1, xIndex: +ij, yIndex: +ikp1);
                                        ij += j;
                                    }

                                    ap[kkp1] += new Complex
                                        (BLAS1Z.zdotc(km1, work, 1, ap, 1, xIndex: +ikp1).Real, 0.0);

                                    ap[kkp1 - 1] += BLAS1Z.zdotc(km1, ap, 1, ap, 1, xIndex: +ik, yIndex: +ikp1);
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
                                        ij += j;
                                    }

                                    ap[kk - 1] += new Complex
                                        (BLAS1Z.zdotc(km1, work, 1, ap, 1, yIndex: +ik).Real, 0.0);
                                    break;
                                }
                            }

                            kstep = 2;
                            break;
                        }
                    }

                    //
                    //  Swap
                    //
                    int ks = Math.Abs(ipvt[k - 1]);

                    if (ks != k)
                    {
                        int iks = ks * (ks - 1) / 2;

                        BLAS1Z.zswap(ks, ref ap, 1, ref ap, 1, xIndex: +iks, yIndex: +ik);
                        int ksj = ik + ks;

                        Complex t2;
                        int jb;
                        for (jb = ks; jb <= k; jb++)
                        {
                            j = k + ks - jb;
                            jk = ik + j;

                            t2 = Complex.Conjugate(ap[jk - 1]);
                            ap[jk - 1] = Complex.Conjugate(ap[ksj - 1]);
                            ap[ksj - 1] = t2;

                            ksj -= (j - 1);
                        }

                        if (kstep != 1)
                        {
                            int kskp1 = ikp1 + ks;

                            t2 = ap[kskp1 - 1];
                            ap[kskp1 - 1] = ap[kkp1 - 1];
                            ap[kkp1 - 1] = t2;
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