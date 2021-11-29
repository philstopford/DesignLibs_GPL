using System;
using System.Numerics;
using Burkardt.BLAS;
using Burkardt.Types;

namespace Burkardt.Linpack;

public static class ZSPCO
{
    public static double zspco(ref Complex[] ap, int n, ref int[] ipvt)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZSPCO factors a Complex symmetric matrix stored in packed form.
        //
        //  Discussion:
        //
        //    The routine also estimates the condition of the matrix.
        //
        //    If RCOND is not needed, ZSPFA is slightly faster.
        //
        //    To solve A*X = B, follow ZSPCO by ZSPSL.
        //
        //    To compute inverse(A)*C, follow ZSPCO by ZSPSL.
        //
        //    To compute inverse(A), follow ZSPCO by ZSPDI.
        //
        //    To compute determinant(A), follow ZSPCO by ZSPDI.
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
        //    Input/output, Complex AP[N*(N+1)/2]; on input, the packed form of a
        //    symmetric matrix.  The columns of the upper triangle are stored
        //    sequentially in a one-dimensional array.  On output, a block diagonal
        //    matrix and the multipliers which were used to obtain it, stored in packed
        //    form.  The factorization can be written A = U*D*U' where U is a product
        //    of permutation and unit upper triangular matrices, U' is the transpose
        //    of U, and D is block diagonal with 1 by 1 and 2 by 2 blocks.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Output, int IPVT[N], the pivot indices.
        //
        //    Output, double ZSPCO, an estimate of RCOND, the reciprocal condition of
        //    the matrix.  For the system A*X = B, relative perturbations in A and B
        //    of size EPSILON may cause relative perturbations in X of size
        //    (EPSILON/RCOND).  If RCOND is so small that the logical expression
        //      1.0 + RCOND == 1.0
        //    is true, then A may be singular to working precision.  In particular,
        //    RCOND is zero if exact singularity is detected or the estimate underflows.
        //
        //  Local Parameters:
        //
        //    Local, Complex Z[N], a work vector whose contents are usually
        //    unimportant.  If A is close to a singular matrix, then Z is an
        //    approximate null vector in the sense that
        //      norm(A*Z) = RCOND * norm(A) * norm(Z).
        //
    {
        Complex ak;
        Complex akm1;
        Complex bk;
        Complex bkm1;
        Complex denom;
        int ikm1;
        int ikp1;
        int j;
        int kk;
        int km1k;
        int km1km1;
        int kp;
        int kps;
        int ks;
        double rcond;
        double s;
        Complex t;

        Complex[] z = new Complex[n];
        //
        //  Find norm of A using only upper half.
        //
        int j1 = 1;

        for (j = 1; j <= n; j++)
        {
            z[j - 1] = new Complex(BLAS1Z.dzasum(j, ap, 1, index: +j1 - 1), 0.0);
            int ij = j1;
            j1 += j;

            int i;
            for (i = 1; i <= j - 1; i++)
            {
                z[i - 1] = new Complex(z[i - 1].Real + typeMethods.zabs1(ap[ij - 1]), 0.0);
                ij += 1;
            }
        }

        double anorm = 0.0;
        for (j = 0; j < n; j++)
        {
            anorm = Math.Max(anorm, z[j].Real);
        }

        //
        //  Factor.
        //
        ZSPFA.zspfa(ref ap, n, ref ipvt);
        //
        //  RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))).
        //
        //  Estimate = norm(Z)/norm(Y) where A*Z = Y and A*Y = E.
        //
        //  The components of E are chosen to cause maximum local
        //  growth in the elements of W where U*D*W = E.
        //
        //  The vectors are frequently rescaled to avoid overflow.
        //
        //  Solve U*D*W = E.
        //
        Complex ek = new(1.0, 0.0);
        for (j = 0; j < n; j++)
        {
            z[j] = new Complex(0.0, 0.0);
        }

        int k = n;
        int ik = n * (n - 1) / 2;

        while (0 < k)
        {
            kk = ik + k;
            ikm1 = ik - (k - 1);

            ks = ipvt[k - 1] switch
            {
                < 0 => 2,
                _ => 1
            };

            kp = Math.Abs(ipvt[k - 1]);
            kps = k + 1 - ks;

            if (kp != kps)
            {
                t = z[kps - 1];
                z[kps - 1] = z[kp - 1];
                z[kp - 1] = t;
            }

            if (typeMethods.zabs1(z[k - 1]) != 0.0)
            {
                ek = typeMethods.zsign1(ek, z[k - 1]);
            }

            z[k - 1] += ek;
            BLAS1Z.zaxpy(k - ks, z[k - 1], ap, 1, ref z, 1, xIndex: +ik);

            if (ks != 1)
            {
                if (typeMethods.zabs1(z[k - 2]) != 0.0)
                {
                    ek = typeMethods.zsign1(ek, z[k - 2]);
                }

                z[k - 2] += ek;
                BLAS1Z.zaxpy(k - ks, z[k - 2], ap, 1, ref z, 1, xIndex: +ikm1);
            }

            if (ks != 2)
            {
                if (typeMethods.zabs1(ap[kk - 1]) < typeMethods.zabs1(z[k - 1]))
                {
                    s = typeMethods.zabs1(ap[kk - 1]) / typeMethods.zabs1(z[k - 1]);
                    BLAS1Z.zdscal(n, s, ref z, 1);
                    ek = new Complex(s, 0.0) * ek;
                }

                if (typeMethods.zabs1(ap[kk - 1]) != 0.0)
                {
                    z[k - 1] /= ap[kk - 1];
                }
                else
                {
                    z[k - 1] = new Complex(1.0, 0.0);
                }
            }
            else
            {
                km1k = ik + k - 1;
                km1km1 = ikm1 + k - 1;
                ak = ap[kk - 1] / ap[km1k - 1];
                akm1 = ap[km1km1 - 1] / ap[km1k - 1];
                bk = z[k - 1] / ap[km1k - 1];
                bkm1 = z[k - 2] / ap[km1k - 1];
                denom = ak * akm1 - new Complex(1.0, 0.0);
                z[k - 1] = (akm1 * bk - bkm1) / denom;
                z[k - 2] = (ak * bkm1 - bk) / denom;
            }

            k -= ks;
            ik -= k;
            switch (ks)
            {
                case 2:
                    ik -= k + 1;
                    break;
            }
        }

        s = 1.0 / BLAS1Z.dzasum(n, z, 1);
        BLAS1Z.zdscal(n, s, ref z, 1);
        //
        //  Solve trans(U) * Y = W.
        //
        k = 1;
        ik = 0;

        while (k <= n)
        {
            ks = ipvt[k - 1] switch
            {
                < 0 => 2,
                _ => 1
            };

            if (k != 1)
            {
                z[k - 1] += BLAS1Z.zdotu(k - 1, ap, 1, z, 1, xIndex: +ik);
                ikp1 = ik + k;

                switch (ks)
                {
                    case 2:
                        z[k] += BLAS1Z.zdotu(k - 1, ap, 1, z, 1, xIndex: +ikp1);
                        break;
                }

                kp = Math.Abs(ipvt[k - 1]);

                if (kp != k)
                {
                    t = z[k - 1];
                    z[k - 1] = z[kp - 1];
                    z[kp - 1] = t;
                }
            }

            ik += k;
            ik = ks switch
            {
                2 => ik + k + 1,
                _ => ik
            };

            k += ks;
        }

        s = 1.0 / BLAS1Z.dzasum(n, z, 1);
        BLAS1Z.zdscal(n, s, ref z, 1);
        double ynorm = 1.0;
        //
        //  Solve U*D*V = Y.
        //
        k = n;
        ik = n * (n - 1) / 2;

        while (0 < k)
        {
            kk = ik + k;
            ikm1 = ik - (k - 1);

            ks = ipvt[k - 1] switch
            {
                < 0 => 2,
                _ => 1
            };

            if (k != ks)
            {
                kp = Math.Abs(ipvt[k - 1]);
                kps = k + 1 - ks;

                if (kp != kps)
                {
                    t = z[kps - 1];
                    z[kps - 1] = z[kp - 1];
                    z[kp - 1] = t;
                }

                BLAS1Z.zaxpy(k - ks, z[k - 1], ap, 1, ref z, 1, xIndex: +ik);

                switch (ks)
                {
                    case 2:
                        BLAS1Z.zaxpy(k - ks, z[k - 2], ap, 1, ref z, 1, xIndex: +ikm1);
                        break;
                }
            }

            if (ks != 2)
            {
                if (typeMethods.zabs1(ap[kk - 1]) < typeMethods.zabs1(z[k - 1]))
                {
                    s = typeMethods.zabs1(ap[kk - 1]) / typeMethods.zabs1(z[k - 1]);
                    BLAS1Z.zdscal(n, s, ref z, 1);
                    ynorm = s * ynorm;
                }

                if (typeMethods.zabs1(ap[kk - 1]) != 0.0)
                {
                    z[k - 1] /= ap[kk - 1];
                }
                else
                {
                    z[k - 1] = new Complex(1.0, 0.0);
                }
            }
            else
            {
                km1k = ik + k - 1;
                km1km1 = ikm1 + k - 1;
                ak = ap[kk - 1] / ap[km1k - 1];
                akm1 = ap[km1km1 - 1] / ap[km1k - 1];
                bk = z[k - 1] / ap[km1k - 1];
                bkm1 = z[k - 2] / ap[km1k - 1];
                denom = ak * akm1 - new Complex(1.0, 0.0);
                z[k - 1] = (akm1 * bk - bkm1) / denom;
                z[k - 2] = (ak * bkm1 - bk) / denom;
            }

            k -= ks;
            ik -= k;

            switch (ks)
            {
                case 2:
                    ik -= k + 1;
                    break;
            }
        }

        s = 1.0 / BLAS1Z.dzasum(n, z, 1);
        BLAS1Z.zdscal(n, s, ref z, 1);
        ynorm = s * ynorm;
        //
        //  Solve U' * Z = V.
        //
        k = 1;
        ik = 0;

        while (k <= n)
        {
            ks = ipvt[k - 1] switch
            {
                < 0 => 2,
                _ => 1
            };

            if (k != 1)
            {
                z[k - 1] += BLAS1Z.zdotu(k - 1, ap, 1, z, 1, xIndex: +ik);
                ikp1 = ik + k;

                switch (ks)
                {
                    case 2:
                        z[k] += BLAS1Z.zdotu(k - 1, ap, 1, z, 1, xIndex: +ikp1);
                        break;
                }

                kp = Math.Abs(ipvt[k - 1]);

                if (kp != k)
                {
                    t = z[k - 1];
                    z[k - 1] = z[kp - 1];
                    z[kp - 1] = t;
                }
            }

            ik += k;

            ik = ks switch
            {
                2 => ik + k + 1,
                _ => ik
            };

            k += ks;
        }

        //
        //  Make ZNORM = 1.
        //
        s = 1.0 / BLAS1Z.dzasum(n, z, 1);
        BLAS1Z.zdscal(n, s, ref z, 1);
        ynorm = s * ynorm;

        if (anorm != 0.0)
        {
            rcond = ynorm / anorm;
        }
        else
        {
            rcond = 0.0;
        }

        return rcond;
    }

}