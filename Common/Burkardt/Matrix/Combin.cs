using System;

namespace Burkardt.MatrixNS;

public static partial class Matrix
{
    public static double[] combin(double alpha, double beta, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COMBIN returns the COMBIN matrix.
        //
        //  Discussion:
        //
        //    This matrix is known as the combinatorial matrix.
        //
        //  Formula:
        //
        //    If ( I = J ) then
        //      A(I,J) = ALPHA + BETA
        //    else
        //      A(I,J) = BETA
        //
        //  Example:
        //
        //    N = 5, ALPHA = 2, BETA = 3
        //
        //    5 3 3 3 3
        //    3 5 3 3 3
        //    3 3 5 3 3
        //    3 3 3 5 3
        //    3 3 3 3 5
        //
        //  Properties:
        //
        //    A is symmetric: A' = A.
        //
        //    Because A is symmetric, it is normal.
        //
        //    Because A is normal, it is diagonalizable.
        //
        //    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).
        //
        //    A is a circulant matrix: each row is shifted once to get the next row.
        //
        //    det ( A ) = ALPHA^(N-1) * ( ALPHA + N * BETA ).
        //
        //    A has constant row sums.
        //
        //    Because A has constant row sums,
        //    it has an eigenvalue with this value,
        //    and a (right) eigenvector of ( 1, 1, 1, ..., 1 ).
        //
        //    A has constant column sums.
        //
        //    Because A has constant column sums,
        //    it has an eigenvalue with this value,
        //    and a (left) eigenvector of ( 1, 1, 1, ..., 1 ).
        //
        //    LAMBDA(1:N-1) = ALPHA,
        //    LAMBDA(N) = ALPHA + N * BETA.
        //
        //    The eigenvector associated with LAMBDA(N) is (1,1,1,...,1)/sqrt(N).
        //
        //    The other N-1 eigenvectors are simply any (orthonormal) basis
        //    for the space perpendicular to (1,1,1,...,1).
        //
        //    A is nonsingular if ALPHA /= 0 and ALPHA + N * BETA /= 0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 June 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Robert Gregory, David Karney,
        //    A Collection of Matrices for Testing Computational Algorithms,
        //    Wiley, 1969,
        //    ISBN: 0882756494,
        //    LC: QA263.68
        //
        //    Donald Knuth,
        //    The Art of Computer Programming,
        //    Volume 1, Fundamental Algorithms, Second Edition,
        //    Addison-Wesley, Reading, Massachusetts, 1973, page 36.
        //
        //  Parameters:
        //
        //    Input, double ALPHA, BETA, scalars that define A.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Output, double COMBIN[N*N], the matrix.
        //
    {
        int j;

        double[] a = new double[n * n];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < n; i++)
            {
                if (i == j)
                {
                    a[i + j * n] = alpha + beta;
                }
                else
                {
                    a[i + j * n] = beta;
                }
            }
        }

        return a;
    }

    public static double[] combin_inverse(double alpha, double beta, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COMBIN_INVERSE returns the inverse of the COMBIN matrix.
        //
        //  Formula:
        //
        //    if ( I = J )
        //      A(I,J) = (ALPHA+(N-1)*BETA) / (ALPHA*(ALPHA+N*BETA))
        //    else
        //      A(I,J) =             - BETA / (ALPHA*(ALPHA+N*BETA))
        //
        //  Example:
        //
        //    N = 5, ALPHA = 2, BETA = 3
        //
        //           14 -3 -3 -3 -3
        //           -3 14 -3 -3 -3
        //   1/34 *  -3 -3 14 -3 -3
        //           -3 -3 -3 14 -3
        //           -3 -3 -3 -3 14
        //
        //  Properties:
        //
        //    A is symmetric: A' = A.
        //
        //    Because A is symmetric, it is normal.
        //
        //    Because A is normal, it is diagonalizable.
        //
        //    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).
        //
        //    A is a circulant matrix: each row is shifted once to get the next row.
        //
        //    A is Toeplitz: constant along diagonals.
        //
        //    det ( A ) = 1 / (ALPHA^(N-1) * (ALPHA+N*BETA)).
        //
        //    A is well defined if ALPHA /= 0 and ALPHA+N*BETA /= 0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 June 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Donald Knuth,
        //    The Art of Computer Programming,
        //    Volume 1, Fundamental Algorithms, Second Edition,
        //    Addison-Wesley, Reading, Massachusetts, 1973, page 36.
        //
        //  Parameters:
        //
        //    Input, double ALPHA, BETA, scalars that define the matrix.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Output, double COMBIN_INVERSE[N*N], the matrix.
        //
    {
        int j;

        switch (alpha)
        {
            case 0.0:
                Console.WriteLine("");
                Console.WriteLine("COMBIN_INVERSE - Fatal error!");
                Console.WriteLine("  The entries of the matrix are undefined");
                Console.WriteLine("  because ALPHA = 0.");
                return null;
        }

        switch (alpha + n * beta)
        {
            case 0.0:
                Console.WriteLine("");
                Console.WriteLine("COMBIN_INVERSE - Fatal error!");
                Console.WriteLine("  The entries of the matrix are undefined");
                Console.WriteLine("  because ALPHA+N*BETA is zero.");
                return null;
        }

        double[] a = new double[n * n];

        double bot = alpha * (alpha + n * beta);

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < n; i++)
            {
                if (i == j)
                {
                    a[i + j * n] = (alpha + (n - 1) * beta) / bot;
                }
                else
                {
                    a[i + j * n] = -beta / bot;
                }
            }
        }

        return a;
    }
}