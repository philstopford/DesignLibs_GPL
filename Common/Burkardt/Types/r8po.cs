using System;
using Burkardt.BLAS;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static int r8po_fa ( ref double[] a, int lda, int n )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PO_FA factors a real symmetric positive definite matrix.
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
        //    FORTRAN77 original version by Dongarra, Moler, Bunch, Stewart.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Dongarra, Moler, Bunch and Stewart,
        //    LINPACK User's Guide,
        //    SIAM, (Society for Industrial and Applied Mathematics),
        //    3600 University City Science Center,
        //    Philadelphia, PA, 19104-2688.
        //    ISBN 0-89871-172-X
        //
        //  Parameters:
        //
        //    Input/output, double A[LDA*N].  On input, the symmetric matrix
        //    to be  factored.  Only the diagonal and upper triangle are used.
        //    On output, an upper triangular matrix R so that A = R'*R
        //    where R' is the transpose.  The strict lower triangle is unaltered.
        //    If INFO /= 0, the factorization is not complete.
        //
        //    Input, int LDA, the leading dimension of the array A.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Output, int R8PO_FA, error flag.
        //    0, for normal return.
        //    K, signals an error condition.  The leading minor of order K is not
        //    positive definite.
        //
        {
            int info;
            int j;
            int k;
            double s;
            double t;

            for ( j = 1; j <= n; j++ )
            {
                s = 0.0;

                for ( k = 1; k <= j-1; k++ )
                {
                    t = a[k-1+(j-1)*lda] - BLAS1D.ddot ( k-1, a, 1, a, 1, +0+(k-1)*lda,+0+(j-1)*lda );
                    t = t / a[k-1+(k-1)*lda];
                    a[k-1+(j-1)*lda] = t;
                    s = s + t * t;
                }

                s = a[j-1+(j-1)*lda] - s;

                if ( s <= 0.0 )
                {
                    info = j;
                    return info;
                }

                a[j-1+(j-1)*lda] = Math.Sqrt ( s );
            }

            info = 0;

            return info;
        }
        
        public static double[] r8po_fa(int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8PO_FA factors an R8PO matrix.
            //
            //  Discussion:
            //
            //    The R8PO storage format is appropriate for a symmetric positive definite 
            //    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
            //    upper triangular matrix, so it will be in R8GE storage format.)
            //
            //    Only the diagonal and upper triangle of the square array are used.
            //    This same storage format is used when the matrix is factored by
            //    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
            //    is set to zero.
            //
            //    The positive definite symmetric matrix A has a Cholesky factorization
            //    of the form:
            //
            //      A = R' * R
            //
            //    where R is an upper triangular matrix with positive elements on
            //    its diagonal.  
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    16 February 2013
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
            //    LINPACK User's Guide,
            //    SIAM, 1979,
            //    ISBN13: 978-0-898711-72-1,
            //    LC: QA214.L56.
            //
            //  Parameters:
            //
            //    Input, int N, the order of the matrix.
            //
            //    Input, double A[N*N], the matrix in R8PO storage.
            //
            //    Output, double R8PO_FA[N*N], the Cholesky factor in SGE
            //    storage, or NULL if there was an error.
            //
        {
            double s;

            double[] b = new double[n * n];

            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i < n; i++)
                {
                    b[i + j * n] = a[i + j * n];
                }
            }

            for (int j = 0; j < n; j++)
            {
                for (int k = 0; k <= j - 1; k++)
                {
                    for (int i = 0; i <= k - 1; i++)
                    {
                        b[k + j * n] = b[k + j * n] - b[i + k * n] * b[i + j * n];
                    }

                    b[k + j * n] = b[k + j * n] / b[k + k * n];
                }

                s = b[j + j * n];
                for (int i = 0; i <= j - 1; i++)
                {
                    s = s - b[i + j * n] * b[i + j * n];
                }

                if (s <= 0.0)
                {
                    return null;
                }

                b[j + j * n] = Math.Sqrt(s);
            }

            //
            //  Since the Cholesky factor is in R8GE format, zero out the lower triangle.
            //
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    b[i + j * n] = 0.0;
                }
            }

            return b;
        }

        public static void r8po_sl(double[] a, int lda, int n, double[] b, int aIndex = 0, int bIndex = 0 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PO_SL solves a linear system factored by R8PO_FA.
        //
        //  Discussion:
        //
        //    A division by zero will occur if the input factor contains
        //    a zero on the diagonal.  Technically this indicates
        //    singularity but it is usually caused by improper subroutine
        //    arguments.  It will not occur if the subroutines are called
        //    correctly and INFO == 0.
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
        //    FORTRAN77 original version by Dongarra, Moler, Bunch, Stewart.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Dongarra, Moler, Bunch and Stewart,
        //    LINPACK User's Guide,
        //    SIAM, (Society for Industrial and Applied Mathematics),
        //    3600 University City Science Center,
        //    Philadelphia, PA, 19104-2688.
        //    ISBN 0-89871-172-X
        //
        //  Parameters:
        //
        //    Input, double A[LDA*N], the output from R8PO_FA.
        //
        //    Input, int LDA, the leading dimension of the array A.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input/output, double B[N].  On input, the right hand side.
        //    On output, the solution.
        //
        {
            int k;
            double t;
            //
            //  Solve R' * Y = B.
            //
            for (k = 1; k <= n; k++)
            {
                t = BLAS1D.ddot(k - 1, a, 1, b, 1,  aIndex + 0 + (k - 1) * lda, bIndex);
                b[bIndex + k - 1] = (b[bIndex + k - 1] - t) / a[aIndex + k - 1 + (k - 1) * lda];
            }

            //
            //  Solve R * X = Y.
            //
            for (k = n; 1 <= k; k--)
            {
                b[bIndex + k - 1] = b[bIndex + k - 1] / a[aIndex + k - 1 + (k - 1) * lda];
                t = -b[bIndex + k - 1];
                BLAS1D.daxpy(k - 1, t, a, 1, ref b, 1, aIndex + 0 + (k - 1) * lda, bIndex);
            }
        }

    }
}