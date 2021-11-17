using System;
using Burkardt.BLAS;

namespace Burkardt.Linpack;

public static class DGBSL
{
    public static void dgbsl(double[] abd, int lda, int n, int ml, int mu, int[] ipvt,
            ref double[] b, int job )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DGBSL solves a real banded system factored by DGBCO or DGBFA.
        //
        //  Discussion:
        //
        //    DGBSL can solve either A * X = B  or  A' * X = B.
        //
        //    A division by zero will occur if the input factor contains a
        //    zero on the diagonal.  Technically this indicates singularity
        //    but it is often caused by improper arguments or improper
        //    setting of LDA.  It will not occur if the subroutines are
        //    called correctly and if DGBCO has set 0.0 < RCOND
        //    or DGBFA has set INFO == 0.
        //
        //    To compute inverse(A) * C  where C is a matrix with P columns:
        //
        //      call dgbco ( abd, lda, n, ml, mu, ipvt, rcond, z )
        //
        //      if ( rcond is too small ) then
        //        exit
        //      end if
        //
        //      do j = 1, p
        //        call dgbsl ( abd, lda, n, ml, mu, ipvt, c(1,j), 0 )
        //      end do
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 May 2005
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
        //    Input, double ABD[LDA*N], the output from DGBCO or DGBFA.
        //
        //    Input, integer LDA, the leading dimension of the array ABD.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int ML, MU, the number of diagonals below and above the
        //    main diagonal.  0 <= ML < N, 0 <= MU < N.
        //
        //    Input, int IPVT[N], the pivot vector from DGBCO or DGBFA.
        //
        //    Input/output, double B[N].  On input, the right hand side.
        //    On output, the solution.
        //
        //    Input, int JOB, job choice.
        //    0, solve A*X=B.
        //    nonzero, solve A'*X=B.
        //
    {
        int k;
        int l;
        int la;
        int lb;
        int lm;
        int m;
        double t;

        m = mu + ml + 1;
        switch (job)
        {
            //
            //  JOB = 0, Solve A * x = b.
            //
            //  First solve L * y = b.
            //
            case 0:
            {
                switch (ml)
                {
                    case > 0:
                    {
                        for (k = 1; k <= n - 1; k++)
                        {
                            lm = Math.Min(ml, n - k);
                            l = ipvt[k - 1];
                            t = b[l - 1];
                            if (l != k)
                            {
                                b[l - 1] = b[k - 1];
                                b[k - 1] = t;
                            }

                            BLAS1D.daxpy(lm, t, abd, 1, ref b, 1, xIndex:  + m + (k - 1) * lda, yIndex:  + k);
                        }

                        break;
                    }
                }

                //
                //  Now solve U * x = y.
                //
                for (k = n; 1 <= k; k--)
                {
                    b[k - 1] /= abd[m - 1 + (k - 1) * lda];
                    lm = Math.Min(k, m) - 1;
                    la = m - lm;
                    lb = k - lm;
                    t = -b[k - 1];
                    BLAS1D.daxpy(lm, t, abd, 1, ref b, 1, xIndex:  + la - 1 + (k - 1) * lda, yIndex:  + lb - 1);
                }

                break;
            }
            //
            default:
            {
                for (k = 1; k <= n; k++)
                {
                    lm = Math.Min(k, m) - 1;
                    la = m - lm;
                    lb = k - lm;
                    t = BLAS1D.ddot(lm, abd, 1, b, 1, xIndex:  + la - 1 + (k - 1) * lda, yIndex:  + lb - 1);
                    b[k - 1] = (b[k - 1] - t) / abd[m - 1 + (k - 1) * lda];
                }

                switch (ml)
                {
                    //
                    //  Now solve L' * x = y.
                    //
                    case > 0:
                    {
                        for (k = n - 1; 1 <= k; k--)
                        {
                            lm = Math.Min(ml, n - k);
                            b[k - 1] += BLAS1D.ddot(lm, abd, 1, b, 1, xIndex:  + m + (k - 1) * lda, yIndex: + k);
                            l = ipvt[k - 1];
                            if (l != k)
                            {
                                t = b[l - 1];
                                b[l - 1] = b[k - 1];
                                b[k - 1] = t;
                            }
                        }

                        break;
                    }
                }

                break;
            }
        }
    }
}