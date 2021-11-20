﻿using Burkardt.BLAS;

namespace Burkardt.Linpack;

public static class DGESL
{
    public static void dgesl(double[] a, int lda, int n, int[] ipvt, ref double[] b, int job )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DGESL solves a real general linear system A * X = B.
        //
        //  Discussion:
        //
        //    DGESL can solve either of the systems A * X = B or A' * X = B.
        //
        //    The system matrix must have been factored by DGECO or DGEFA.
        //
        //    A division by zero will occur if the input factor contains a
        //    zero on the diagonal.  Technically this indicates singularity
        //    but it is often caused by improper arguments or improper
        //    setting of LDA.  It will not occur if the subroutines are
        //    called correctly and if DGECO has set 0.0 < RCOND
        //    or DGEFA has set INFO == 0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 May 2005
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
        //    Input, double A[LDA*N], the output from DGECO or DGEFA.
        //
        //    Input, int LDA, the leading dimension of A.
        //
        //    Input, int N, the order of the matrix A.
        //
        //    Input, int IPVT[N], the pivot vector from DGECO or DGEFA.
        //
        //    Input/output, double B[N].
        //    On input, the right hand side vector.
        //    On output, the solution vector.
        //
        //    Input, int JOB.
        //    0, solve A * X = B;
        //    nonzero, solve A' * X = B.
        //
    {
        int k;
        int l;
        double t;
        switch (job)
        {
            //
            //  Solve A * X = B.
            //
            case 0:
            {
                for (k = 1; k <= n - 1; k++)
                {
                    l = ipvt[k - 1];
                    t = b[l - 1];

                    if (l != k)
                    {
                        b[l - 1] = b[k - 1];
                        b[k - 1] = t;
                    }

                    BLAS1D.daxpy(n - k, t, a, 1, ref b, 1, xIndex: + k + (k - 1) * lda, yIndex: + k);

                }

                for (k = n; 1 <= k; k--)
                {
                    b[k - 1] /= a[k - 1 + (k - 1) * lda];
                    t = -b[k - 1];
                    BLAS1D.daxpy(k - 1, t, a, 1, ref b, 1, xIndex: + 0 + (k - 1) * lda);
                }

                break;
            }
            //
            default:
            {
                for (k = 1; k <= n; k++)
                {
                    t = BLAS1D.ddot(k - 1, a, 1, b, 1, xIndex: + 0 + (k - 1) * lda);
                    b[k - 1] = (b[k - 1] - t) / a[k - 1 + (k - 1) * lda];
                }

                for (k = n - 1; 1 <= k; k--)
                {
                    b[k - 1] += BLAS1D.ddot(n - k, a, 1, b, 1, xIndex:  + k + (k - 1) * lda, yIndex: + k);
                    l = ipvt[k - 1];

                    if (l == k)
                    {
                        continue;
                    }

                    t = b[l - 1];
                    b[l - 1] = b[k - 1];
                    b[k - 1] = t;
                }

                break;
            }
        }
    }
}