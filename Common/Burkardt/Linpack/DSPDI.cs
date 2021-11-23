﻿using System;
using Burkardt.BLAS;

namespace Burkardt.Linpack;

public static class DSPDI
{
    public static void dspdi(ref double[] ap, int n, int[] kpvt, ref double[] det, ref int[] inert,
            double[] work, int job)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DSPDI computes the determinant, inertia and inverse of a real symmetric matrix.
        //
        //  Discussion:
        //
        //    DSPDI uses the factors from DSPFA, where the matrix is stored in
        //    packed form.
        //
        //    A division by zero will occur if the inverse is requested
        //    and DSPCO has set RCOND == 0.0D+00 or DSPFA has set INFO /= 0.
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
        //    Input/output, double AP[(N*(N+1))/2].  On input, the output from
        //    DSPFA.  On output, the upper triangle of the inverse of the original
        //    matrix, stored in packed form, if requested.  The columns of the upper 
        //    triangle are stored sequentially in a one-dimensional array.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int KPVT[N], the pivot vector from DSPFA.
        //
        //    Output, double DET[2], the determinant of the original matrix,
        //    if requested.
        //      determinant = DET[0] * 10.0**DET[1]
        //    with 1.0D+00 <= abs ( DET[0] ) < 10.0D+00 or DET[0] = 0.0.
        //
        //    Output, int INERT[3], the inertia of the original matrix, if requested.
        //    INERT(1) = number of positive eigenvalues.
        //    INERT(2) = number of negative eigenvalues.
        //    INERT(3) = number of zero eigenvalues.
        //
        //    Workspace, double WORK[N].
        //
        //    Input, int JOB, has the decimal expansion ABC where:
        //      if A /= 0, the inertia is computed,
        //      if B /= 0, the determinant is computed,
        //      if C /= 0, the inverse is computed.
        //    For example, JOB = 111  gives all three.
        //
    {
        double d;
        int ik;
        int ikp1;
        int k;
        int kk;
        int kkp1;
        double t;

        bool doinv = job % 10 != 0;
        bool dodet = job % 100 / 10 != 0;
        bool doert = job % 1000 / 100 != 0;

        if (dodet || doert)
        {
            switch (doert)
            {
                case true:
                    inert[0] = 0;
                    inert[1] = 0;
                    inert[2] = 0;
                    break;
            }

            switch (dodet)
            {
                case true:
                    det[0] = 1.0;
                    det[1] = 0.0;
                    break;
            }

            t = 0.0;
            ik = 0;

            for (k = 1; k <= n; k++)
            {
                kk = ik + k;
                d = ap[kk - 1];
                switch (kpvt[k - 1])
                {
                    //
                    //  2 by 2 block
                    //  use det (d  s)  =  (d/t * c - t) * t,  t = abs ( s )
                    //          (s  c)
                    //  to avoid underflow/overflow troubles.
                    //
                    //  Take two passes through scaling.  Use T for flag.
                    //
                    case <= 0 when t == 0.0:
                        ikp1 = ik + k;
                        kkp1 = ikp1 + k;
                        t = Math.Abs(ap[kkp1 - 1]);
                        d = d / t * ap[kkp1] - t;
                        break;
                    case <= 0:
                        d = t;
                        t = 0.0;
                        break;
                }

                switch (doert)
                {
                    case true:
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

                switch (dodet)
                {
                    case true:
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

        switch (doinv)
        {
            //
            //  Compute inverse(A).
            //
            case true:
            {
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
                    switch (kpvt[k - 1])
                    {
                        case >= 0:
                        {
                            //
                            //  1 by 1.
                            //
                            ap[kk - 1] = 1.0 / ap[kk - 1];

                            switch (k)
                            {
                                case >= 2:
                                {
                                    BLAS1D.dcopy(k - 1, ap, 1, ref work, 1, xIndex: +ik);
                                    ij = 0;

                                    for (j = 1; j <= k - 1; j++)
                                    {
                                        jk = ik + j;
                                        ap[jk - 1] = BLAS1D.ddot(j, ap, 1, work, 1, xIndex: +ij);
                                        BLAS1D.daxpy(j - 1, work[j - 1], ap, 1, ref ap, 1, xIndex: +ij, yIndex: +ik);
                                        ij += j;
                                    }

                                    ap[kk - 1] += BLAS1D.ddot(k - 1, work, 1, ap, 1, yIndex: +ik);
                                    break;
                                }
                            }

                            kstep = 1;
                            break;
                        }
                        default:
                        {
                            //
                            //  2 by 2.
                            //
                            t = Math.Abs(ap[kkp1 - 1]);
                            double ak = ap[kk - 1] / t;
                            double akp1 = ap[kkp1] / t;
                            double akkp1 = ap[kkp1 - 1] / t;
                            d = t * (ak * akp1 - 1.0);
                            ap[kk - 1] = akp1 / d;
                            ap[kkp1] = ak / d;
                            ap[kkp1 - 1] = -akkp1 / d;

                            switch (km1)
                            {
                                case >= 1:
                                {
                                    BLAS1D.dcopy(km1, ap, 1, ref work, 1, xIndex: +ikp1);
                                    ij = 0;

                                    for (j = 1; j <= km1; j++)
                                    {
                                        int jkp1 = ikp1 + j;
                                        ap[jkp1 - 1] = BLAS1D.ddot(j, ap, 1, work, 1, xIndex: +ij);
                                        BLAS1D.daxpy(j - 1, work[j - 1], ap, 1, ref ap, 1, xIndex: +ij, yIndex: +ikp1);
                                        ij += j;
                                    }

                                    ap[kkp1] += BLAS1D.ddot(km1, work, 1, ap, 1, yIndex: +ikp1);
                                    ap[kkp1 - 1] += BLAS1D.ddot(km1, ap, 1, ap, 1, xIndex: +ik, yIndex: +ikp1);
                                    BLAS1D.dcopy(km1, ap, 1, ref work, 1, xIndex: +ik);
                                    ij = 0;

                                    for (j = 1; j <= km1; j++)
                                    {
                                        jk = ik + j;
                                        ap[jk - 1] = BLAS1D.ddot(j, ap, 1, work, 1, xIndex: +ij);
                                        BLAS1D.daxpy(j - 1, work[j - 1], ap, 1, ref ap, 1, xIndex: +ij, yIndex: +ik);
                                        ij += j;
                                    }

                                    ap[kk - 1] += BLAS1D.ddot(km1, work, 1, ap, 1, yIndex: +ik);
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
                    int ks = Math.Abs(kpvt[k - 1]);

                    if (ks != k)
                    {
                        int iks = ks * (ks - 1) / 2;
                        BLAS1D.dswap(ks, ref ap, 1, ref ap, 1, xIndex: +iks, yIndex: +ik);
                        int ksj = ik + ks;

                        double temp;
                        int jb;
                        for (jb = ks; jb <= k; jb++)
                        {
                            j = k + ks - jb;
                            jk = ik + j;
                            temp = ap[jk - 1];
                            ap[jk - 1] = ap[ksj - 1];
                            ap[ksj - 1] = temp;
                            ksj -= j - 1;
                        }

                        if (kstep != 1)
                        {
                            int kskp1 = ikp1 + ks;
                            temp = ap[kskp1 - 1];
                            ap[kskp1 - 1] = ap[kkp1 - 1];
                            ap[kkp1 - 1] = temp;
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