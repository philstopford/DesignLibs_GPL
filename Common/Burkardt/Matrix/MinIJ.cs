using System;

namespace Burkardt.MatrixNS
{
    public static partial class Matrix
    {
        public static double[] minij(int m, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MINIJ returns the MINIJ matrix.
            //
            //  Discussion:
            //
            //    A(I,J) = min ( I, J )
            //
            //  Example:
            //
            //    N = 5
            //
            //    1 1 1 1 1
            //    1 2 2 2 2
            //    1 2 3 3 3
            //    1 2 3 4 4
            //    1 2 3 4 5
            //
            //  Properties:
            //
            //    A is integral, therefore det ( A ) is integral, and 
            //    det ( A ) * inverse ( A ) is integral.
            //
            //    A is positive definite.
            //
            //    A is symmetric: A' = A.
            //
            //    Because A is symmetric, it is normal.
            //
            //    Because A is normal, it is diagonalizable.
            //
            //    The inverse of A is tridiagonal.
            //
            //    The eigenvalues of A are
            //
            //      LAMBDA[i] = 0.5 / ( 1 - cos ( ( 2 * I - 1 ) * Math.PI / ( 2 * N + 1 ) ) ),
            //
            //    (N+1)*ONES[N] - A also has a tridiagonal inverse.
            //
            //    Gregory and Karney consider the matrix defined by
            //
            //      B(I,J) = N + 1 - MAX(I,J)
            //
            //    which is equal to the MINIJ matrix, but with the rows and
            //    columns reversed.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 October 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Robert Gregory, David Karney,
            //    Example 3.12, Example 4.14,
            //    A Collection of Matrices for Testing Computational Algorithms,
            //    Wiley, 1969, page 41, page 74, 
            //    LC: QA263.G68.
            //
            //    Daniel Rutherford,
            //    Some continuant determinants arising in physics and chemistry II,
            //    Proceedings of the Royal Society Edinburgh,
            //    Volume 63, A, 1952, pages 232-241.
            //
            //    John Todd,
            //    Basic Numerical Mathematics, Vol. 2: Numerical Algebra,
            //    Academic Press, 1977, page 158.
            //
            //    Joan Westlake,
            //    A Handbook of Numerical Matrix Inversion and Solution of 
            //    Linear Equations,
            //    John Wiley, 1968,
            //    ISBN13: 978-0471936756,
            //    LC: QA263.W47.
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns 
            //    of the matrix.
            //
            //    Output, double MINIJ[M*N], the matrix.
            //
        {
            double[] a;
            int i;
            int j;

            a = new double[m * n];

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    a[i + j * m] = (double) (Math.Min(i + 1, j + 1));
                }
            }

            return a;
        }
    }
}