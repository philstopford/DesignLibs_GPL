using System;
using System.Numerics;

namespace Burkardt.MatrixNS;

public static class Tridiagonal
{
    public static double[] tris(int m, int n, double x, double y, double z)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIS returns the TRIS matrix.
        //
        //  Discussion:
        //
        //    The matrix is a tridiagonal matrix defined by three scalars.
        //
        //    See page 155 of the Todd reference.
        //
        //  Formula:
        //
        //    if ( J = I-1 )
        //      A(I,J) = X
        //    else if ( J = I )
        //      A(I,J) = Y
        //    else if ( J = I + 1 )
        //      A(I,J) = Z
        //    else
        //      A(I,J) = 0
        //
        //  Example:
        //
        //    M = 5, N = 5, X = 1, Y = 2, Z = 3
        //
        //    2 3 0 0 0
        //    1 2 3 0 0
        //    0 1 2 3 0
        //    0 0 1 2 3
        //    0 0 0 1 2
        //
        //  Properties:
        //
        //    A is generally not symmetric: A' /= A.
        //
        //    A is tridiagonal.
        //
        //    Because A is tridiagonal, it has property A (bipartite).
        //
        //    A is banded, with bandwidth 3.
        //
        //    A is Toeplitz: constant along diagonals.
        //
        //    If Y is not zero, then for A to be singular, it must be the case that
        //
        //      0.5 * Y / sqrt ( X * Z ) < 1
        //
        //    and
        //
        //      cos (K*PI/(N+1)) = - 0.5 * Y / sqrt ( X * Z ) for some 1 <= K <= N.
        //
        //    If Y is zero, then A is singular when N is odd, or if X or Z is zero.
        //
        //    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).
        //
        //    A has eigenvalues
        //
        //      LAMBDA(I) = Y + 2 * sqrt(X*Z) * COS(I*PI/(N+1))
        //
        //    The eigenvalues will be complex if X * Z < 0.
        //
        //    If X = Z, the matrix is symmetric.
        //
        //    As long as X and Z are nonzero, the matrix is irreducible.
        //
        //    If X = Z = -1, and Y = 2, the matrix is a symmetric, positive
        //    definite M matrix, the negative of the second difference matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 July 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    John Todd,
        //    Basic Numerical Mathematics,
        //    Volume 2: Numerical Algebra,
        //    Birkhauser, 1980,
        //    ISBN: 0817608117,
        //    LC: QA297.T58.
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Input, double X, Y, Z, the scalars that define A.
        //
        //    Output, double TRIS[M*N], the matrix.
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
                if (j == i - 1)
                {
                    a[i + j * m] = x;
                }
                else if (j == i)
                {
                    a[i + j * m] = y;
                }
                else if (j == i + 1)
                {
                    a[i + j * m] = z;
                }
                else
                {
                    a[i + j * m] = 0.0;
                }
            }
        }

        return a;
    }

    public static Complex[] tris_eigenvalues(int n, double x, double y, double z)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIS_EIGENVALUES returns the eigenvalues of the TRIS matrix.
        //
        //  Discussion:
        //
        //    The eigenvalues will be complex if X * Z < 0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 July 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double X, Y, Z, the scalars that define A.
        //
        //    Output, complex <double> TRIS_EIGENVALUES[N], the eigenvalues.
        //
    {
        double angle;
        Complex arg;
        int i;
        Complex[] lambda;
            

        lambda = new Complex[n];

        for (i = 0; i < n; i++)
        {
            angle = (i + 1) * Math.PI / (n + 1);
            arg = new Complex(x * z, 0.0);
            lambda[i] = y + 2.0 * Complex.Sqrt(arg) * Math.Cos(angle);
        }

        return lambda;
    }
}