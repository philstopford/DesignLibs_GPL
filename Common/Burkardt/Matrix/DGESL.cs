using Burkardt.BLAS;

namespace Burkardt.MatrixNS;

public static partial class Matrix
{
    public static void dge_sl ( int n, double[] a, int[] pivot, ref double[] b, int job )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DGE_SL solves a system factored by SGE_FA.
        //
        //  Discussion:
        //
        //    DGE_SL is a simplified version of the LINPACK routine SGESL.
        //
        //    The doubly dimensioned array A is treated as a one dimensional vector,
        //    stored by COLUMNS:
        //
        //      A(0,0), A(1,0), A(2,0), ..., A(N-1,0) // A(1,0), A(1,1), ... A(N-1,1)
        //
        //    Entry A(I,J) is stored as A[I+J*N]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be positive.
        //
        //    Input, double A[N*N], the LU factors from DGE_FA.
        //
        //    Input, int PIVOT[N], the pivot vector from DGE_FA.
        //
        //    Input/output, double B[N].
        //    On input, the right hand side vector.
        //    On output, the solution vector.
        //
        //    Input, int JOB, specifies the operation.
        //    0, solve A * x = b.
        //    nonzero, solve A' * x = b.
        //
    {
        int i;
        int k;
        int l;
        double t;
        switch (job)
        {
            //
            //  Solve A * x = b.
            //
            case 0:
            {
                //
                //  Solve PL * Y = B.
                //
                for (k = 1; k <= n - 1; k++)
                {
                    l = pivot[k - 1];

                    if (l != k)
                    {
                        t = b[l - 1];
                        b[l - 1] = b[k - 1];
                        b[k - 1] = t;
                    }

                    for (i = k + 1; i <= n; i++)
                    {
                        b[i - 1] += a[i - 1 + (k - 1) * n] * b[k - 1];
                    }
                }

                //
                //  Solve U * X = Y.
                //
                for (k = n; 1 <= k; k--)
                {
                    b[k - 1] /= a[k - 1 + (k - 1) * n];
                    for (i = 1; i <= k - 1; i++)
                    {
                        b[i - 1] -= a[i - 1 + (k - 1) * n] * b[k - 1];
                    }
                }
                //
                //  Solve A' * X = B.
                //
                break;
            }
            default:
            {
                //
                //  Solve U' * Y = B.
                //
                for (k = 1; k <= n; k++)
                {
                    t = 0.0;
                    for (i = 1; i <= k - 1; i++)
                    {
                        t += b[i - 1] * a[i - 1 + (k - 1) * n];
                    }

                    b[k - 1] = (b[k - 1] - t) / a[k - 1 + (k - 1) * n];
                }

                //
                //  Solve ( PL )' * X = Y.
                //
                for (k = n - 1; 1 <= k; k--)
                {
                    t = 0.0;
                    for (i = k + 1; i <= n; i++)
                    {
                        t += b[i - 1] * a[i - 1 + (k - 1) * n];
                    }

                    b[k - 1] += t;

                    l = pivot[k - 1];

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
    }

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

                    BLAS1D.daxpy(n - k, t, a, 1, ref b, 1,  + k + (k - 1) * lda, + k);

                }

                for (k = n; 1 <= k; k--)
                {
                    b[k - 1] /= a[k - 1 + (k - 1) * lda];
                    t = -b[k - 1];
                    BLAS1D.daxpy(k - 1, t, a, 1, ref b, 1,  + 0 + (k - 1) * lda);
                }

                break;
            }
            //
            default:
            {
                for (k = 1; k <= n; k++)
                {
                    t = BLAS1D.ddot(k - 1, a, 1, b, 1,  + 0 + (k - 1) * lda);
                    b[k - 1] = (b[k - 1] - t) / a[k - 1 + (k - 1) * lda];
                }

                for (k = n - 1; 1 <= k; k--)
                {
                    b[k - 1] += BLAS1D.ddot(n - k, a, 1, b, 1,  + k + (k - 1) * lda, + k);
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
    }
}