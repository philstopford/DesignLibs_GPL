using System;
using Burkardt.BLAS;
using Burkardt.Types;

namespace Burkardt.Linpack
{
    public static class DGBCO
    {
        public static double dgbco(ref double[] abd, int lda, int n, int ml, int mu, ref int[] ipvt,
        double[] z )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DGBCO factors a real band matrix and estimates its condition.
        //
        //  Discussion:
        //
        //    If RCOND is not needed, DGBFA is slightly faster.
        //
        //    To solve A*X = B, follow DGBCO by DGBSL.
        //
        //    To compute inverse(A)*C, follow DGBCO by DGBSL.
        //
        //    To compute determinant(A), follow DGBCO by DGBDI.
        //
        //  Example:
        //
        //    If the original matrix is
        //
        //      11 12 13  0  0  0
        //      21 22 23 24  0  0
        //       0 32 33 34 35  0
        //       0  0 43 44 45 46
        //       0  0  0 54 55 56
        //       0  0  0  0 65 66
        //
        //    then for proper band storage,
        //
        //      N = 6, ML = 1, MU = 2, 5 <= LDA and ABD should contain
        //
        //       *  *  *  +  +  +      * = not used
        //       *  * 13 24 35 46      + = used for pivoting
        //       * 12 23 34 45 56
        //      11 22 33 44 55 66
        //      21 32 43 54 65  *
        //
        //  Band storage:
        //
        //    If A is a band matrix, the following program segment
        //    will set up the input.
        //
        //      ml = (band width below the diagonal)
        //      mu = (band width above the diagonal)
        //      m = ml + mu + 1
        //
        //      do j = 1, n
        //        i1 = max ( 1, j-mu )
        //        i2 = min ( n, j+ml )
        //        do i = i1, i2
        //          k = i - j + m
        //          abd(k,j) = a(i,j)
        //        }
        //      }
        //
        //    This uses rows ML+1 through 2*ML+MU+1 of ABD.  In addition, the first
        //    ML rows in ABD are used for elements generated during the
        //    triangularization.  The total number of rows needed in ABD is
        //    2*ML+MU+1.  The ML+MU by ML+MU upper left triangle and the ML by ML
        //    lower right triangle are not referenced.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 June 2005
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
        //    Input/output, double ABD[LDA*N].  On input, the matrix in band
        //    storage.  The columns of the matrix are stored in the columns of ABD and
        //    the diagonals of the matrix are stored in rows ML+1 through 2*ML+MU+1
        //    of ABD.  On output, an upper triangular matrix in band storage and
        //    the multipliers which were used to obtain it.  The factorization can
        //    be written A = L*U where L is a product of permutation and unit lower
        //    triangular matrices and U is upper triangular.
        //
        //    Input, int LDA, the leading dimension of the array ABD.
        //    2*ML + MU + 1 <= LDA is required.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int ML, MU, the number of diagonals below and above the
        //    main diagonal.  0 <= ML < N, 0 <= MU < N.
        //
        //    Output, int IPVT[N], the pivot indices.
        //
        //    Workspace, double Z[N], a work vector whose contents are
        //    usually unimportant.  If A is close to a singular matrix, then Z is an
        //    approximate null vector in the sense that
        //      norm(A*Z) = RCOND * norm(A) * norm(Z).
        //
        //    Output, double DGBCO, an estimate of the reciprocal condition number RCOND
        //    of A.  For the system A*X = B, relative perturbations in A and B of size
        //    EPSILON may cause relative perturbations in X of size EPSILON/RCOND.
        //    If RCOND is so small that the logical expression
        //      1.0 + RCOND == 1.0D+00
        //    is true, then A may be singular to working precision.  In particular,
        //    RCOND is zero if exact singularity is detected or the estimate underflows.
        //
        {
            double anorm;
            double ek;
            int i;
            int is_;
            int j;
            int ju;
            int k;
            int l;
            int la;
            int lm;
            int lz;
            int m;
            int mm;
            double rcond;
            double s;
            double sm;
            double t;
            double wk;
            double wkm;
            double ynorm;
            //
            //  Compute the 1-norm of A.
            //
            anorm = 0.0;
            l = ml + 1;
                is_ = l + mu;

            for (j = 1; j <= n; j++)
            {
                anorm = Math.Max(anorm, BLAS1D.dasum(l, abd, 1, index: + is_ -1 + (j - 1) * lda));
                if (ml + 1 < is_)
                {
                    is_ = is_ -1;
                }

                if (j <= mu)
                {
                    l = l + 1;
                }

                if (n - ml <= j)
                {
                    l = l - 1;
                }
            }

            //
            //  Factor.
            //
            DGBFA.dgbfa(ref abd, lda, n, ml, mu, ref ipvt);
            //
            //  RCOND = 1/(norm(A)*(estimate of norm(inverse(A)))).
            //
            //  Estimate = norm(Z)/norm(Y) where  a*z = y  and A'*Y = E.
            //
            //  A' is the transpose of A.  The components of E are
            //  chosen to cause maximum local growth in the elements of W where
            //  U'*W = E.  The vectors are frequently rescaled to avoid
            //  overflow.
            //
            //  Solve U' * W = E.
            //
            ek = 1.0;
            for (i = 1; i <= n; i++)
            {
                z[i - 1] = 0.0;
            }

            m = ml + mu + 1;
            ju = 0;

            for (k = 1; k <= n; k++)
            {
                if (z[k - 1] != 0.0)
                {
                    ek = ek * typeMethods.r8_sign(-z[k - 1]);
                }

                if (Math.Abs(abd[m - 1 + (k - 1) * lda]) < Math.Abs(ek - z[k - 1]))
                {
                    s = Math.Abs(abd[m - 1 + (k - 1) * lda]) / Math.Abs(ek - z[k - 1]);
                    for (i = 1; i <= n; i++)
                    {
                        z[i - 1] = s * z[i - 1];
                    }

                    ek = s * ek;
                }

                wk = ek - z[k - 1];
                wkm = -ek - z[k - 1];
                s = Math.Abs(wk);
                sm = Math.Abs(wkm);

                if (abd[m - 1 + (k - 1) * lda] != 0.0)
                {
                    wk = wk / abd[m - 1 + (k - 1) * lda];
                    wkm = wkm / abd[m - 1 + (k - 1) * lda];
                }
                else
                {
                    wk = 1.0;
                    wkm = 1.0;
                }

                ju = Math.Min(Math.Max(ju, mu + ipvt[k - 1]), n);
                mm = m;

                if (k + 1 <= ju)
                {
                    for (j = k + 1; j <= ju; j++)
                    {
                        mm = mm - 1;
                        sm = sm + Math.Abs(z[j - 1] + wkm * abd[mm - 1 + (j - 1) * lda]);
                        z[j - 1] = z[j - 1] + wk * abd[mm - 1 + (j - 1) * lda];
                        s = s + Math.Abs(z[j - 1]);
                    }

                    if (s < sm)
                    {
                        t = wkm - wk;
                        wk = wkm;
                        mm = m;
                        for (j = k + 1; j <= ju; ju++)
                        {
                            mm = mm - 1;
                            z[j - 1] = z[j - 1] + t * abd[mm - 1 + (j - 1) * lda];
                        }
                    }
                }

                z[k - 1] = wk;
            }

            s = BLAS1D.dasum(n, z, 1);

            for (i = 1; i <= n; i++)
            {
                z[i - 1] = z[i - 1] / s;
            }

            //
            //  Solve L' * Y = W.
            //
            for (k = n; 1 <= k; k--)
            {
                lm = Math.Min(ml, n - k);

                if (k < m)
                {
                    z[k - 1] = z[k - 1] + BLAS1D.ddot(lm, abd, 1, z, 1, xIndex: + m + (k - 1) * lda, yIndex: + k);
                }

                if (1.0 < Math.Abs(z[k - 1]))
                {
                    s = 1.0 / Math.Abs(z[k - 1]);
                    for (i = 1; i <= n; i++)
                    {
                        z[i - 1] = s * z[i - 1];
                    }
                }

                l = ipvt[k - 1];
                t = z[l - 1];
                z[l - 1] = z[k - 1];
                z[k - 1] = t;
            }

            s = BLAS1D.dasum(n, z, 1);
            for (i = 1; i <= n; i++)
            {
                z[i - 1] = z[i - 1] / s;
            }

            ynorm = 1.0;
            //
            //  Solve L * V = Y.
            //
            for (k = 1; k <= n; k++)
            {
                l = ipvt[k - 1];
                t = z[l - 1];
                z[l - 1] = z[k - 1];
                z[k - 1] = t;
                lm = Math.Min(ml, n - k);

                if (k < n)
                {
                    BLAS1D.daxpy(lm, t, abd, 1, ref z, 1, xIndex: + m + (k - 1) * lda, yIndex: + k);
                }

                if (1.0 < Math.Abs(z[k - 1]))
                {
                    s = 1.0 / Math.Abs(z[k - 1]);
                    for (i = 1; i <= n; i++)
                    {
                        z[i - 1] = s * z[i - 1];
                    }

                    ynorm = s * ynorm;
                }
            }

            s = 1.0 / BLAS1D.dasum(n, z, 1);
            for (i = 1; i <= n; i++)
            {
                z[i - 1] = s * z[i - 1];
            }

            ynorm = s * ynorm;
            //
            //  Solve U * Z = W.
            //
            for (k = n; 1 <= k; k--)
            {
                if (Math.Abs(abd[m - 1 + (k - 1) * lda]) < Math.Abs(z[k - 1]))
                {
                    s = Math.Abs(abd[m - 1 + (k - 1) * lda]) / Math.Abs(z[k - 1]);
                    for (i = 1; i <= n; i++)
                    {
                        z[i - 1] = s * z[i - 1];
                    }

                    ynorm = s * ynorm;
                }

                if (abd[m - 1 + (k - 1) * lda] != 0.0)
                {
                    z[k - 1] = z[k - 1] / abd[m - 1 + (k - 1) * lda];
                }
                else
                {
                    z[k - 1] = 1.0;
                }

                lm = Math.Min(k, m) - 1;
                la = m - lm;
                lz = k - lm;
                t = -z[k - 1];
                BLAS1D.daxpy(lm, t, abd, 1, ref z, 1, xIndex:  + la - 1 + (k - 1) * lda, yIndex:  + lz - 1);
            }

            //
            //  Make ZNORM = 1.0.
            //
            s = 1.0 / BLAS1D.dasum(n, z, 1);
            for (i = 1; i <= n; i++)
            {
                z[i - 1] = s * z[i - 1];
            }

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
}