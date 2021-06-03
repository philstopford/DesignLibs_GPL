using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
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

    }
}
