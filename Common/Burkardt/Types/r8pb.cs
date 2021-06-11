namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static void r8pbu_cg(int n, int mu, double[] a, double[] b, double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PBU_CG uses the conjugate gradient method on a R8PBU system.
        //
        //  Discussion:
        //
        //    The R8PBU storage format is used for a symmetric positive definite band matrix.
        //
        //    To save storage, only the diagonal and upper triangle of A is stored,
        //    in a compact diagonal format that preserves columns.
        //
        //    The diagonal is stored in row MU+1 of the array.
        //    The first superdiagonal in row MU, columns 2 through N.
        //    The second superdiagonal in row MU-1, columns 3 through N.
        //    The MU-th superdiagonal in row 1, columns MU+1 through N.
        //
        //    The matrix A must be a positive definite symmetric band matrix.
        //
        //    The method is designed to reach the solution after N computational
        //    steps.  However, roundoff may introduce unacceptably large errors for
        //    some problems.  In such a case, calling the routine again, using
        //    the computed solution as the new starting estimate, should improve
        //    the results.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 February 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Frank Beckman,
        //    The Solution of Linear Equations by the Conjugate Gradient Method,
        //    in Mathematical Methods for Digital Computers,
        //    edited by John Ralston, Herbert Wilf,
        //    Wiley, 1967,
        //    ISBN: 0471706892,
        //    LC: QA76.5.R3.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be positive.
        //
        //    Input, int MU, the number of superdiagonals.
        //    MU must be at least 0, and no more than N-1.
        //
        //    Input, double A[(MU+1)*N], the R8PBU matrix.
        //
        //    Input, double B[N], the right hand side vector.
        //
        //    Input/output, double X[N].
        //    On input, an estimate for the solution.
        //    On output, the approximate solution vector.
        //
        {
            double alpha;
            double[] ap;
            double beta;
            int i;
            int it;
            double[] p;
            double pap;
            double pr;
            double[] r;
            double rap;
            //
            //  Initialize
            //    AP = A * x,
            //    R  = b - A * x,
            //    P  = b - A * x.
            //
            ap = r8pbu_mv(n, n, mu, a, x);

            r = new double[n];
            for (i = 0; i < n; i++)
            {
                r[i] = b[i] - ap[i];
            }

            p = new double[n];
            for (i = 0; i < n; i++)
            {
                p[i] = b[i] - ap[i];
            }

            //
            //  Do the N steps of the conjugate gradient method.
            //
            for (it = 1; it <= n; it++)
            {
                //
                //  Compute the matrix*vector product AP=A*P.
                //
                ap = r8pbu_mv(n, n, mu, a, p);
                //
                //  Compute the dot products
                //    PAP = P*AP,
                //    PR  = P*R
                //  Set
                //    ALPHA = PR / PAP.
                //
                pap = 0.0;
                for (i = 0; i < n; i++)
                {
                    pap = pap + p[i] * ap[i];
                }

                if (pap == 0.0)
                {
                    break;
                }

                pr = 0.0;
                for (i = 0; i < n; i++)
                {
                    pr = pr + p[i] * r[i];
                }

                alpha = pr / pap;
                //
                //  Set
                //    X = X + ALPHA * P
                //    R = R - ALPHA * AP.
                //
                for (i = 0; i < n; i++)
                {
                    x[i] = x[i] + alpha * p[i];
                }

                for (i = 0; i < n; i++)
                {
                    r[i] = r[i] - alpha * ap[i];
                }

                //
                //  Compute the vector dot product
                //    RAP = R*AP
                //  Set
                //    BETA = - RAP / PAP.
                //
                rap = 0.0;
                for (i = 0; i < n; i++)
                {
                    rap = rap + r[i] * ap[i];
                }

                beta = -rap / pap;
                //
                //  Update the perturbation vector
                //    P = R + BETA * P.
                //
                for (i = 0; i < n; i++)
                {
                    p[i] = r[i] + beta * p[i];
                }

            }
        }

        public static double[] r8pbu_dif2(int m, int n, int mu)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8PBU_DIF2 returns the DIF2 matrix in R8PBU format.
            //
            //  Example:
            //
            //    N = 5
            //
            //    2 -1  .  .  .
            //   -1  2 -1  .  .
            //    . -1  2 -1  .
            //    .  . -1  2 -1
            //    .  .  . -1  2
            //
            //  Properties:
            //
            //    A is banded, with bandwidth 3.
            //
            //    A is tridiagonal.
            //
            //    Because A is tridiagonal, it has property A (bipartite).
            //
            //    A is a special case of the TRIS or tridiagonal scalar matrix.
            //
            //    A is integral, therefore det ( A ) is integral, and 
            //    det ( A ) * inverse ( A ) is integral.
            //
            //    A is Toeplitz: constant along diagonals.
            //
            //    A is symmetric: A' = A.
            //
            //    Because A is symmetric, it is normal.
            //
            //    Because A is normal, it is diagonalizable.
            //
            //    A is persymmetric: A(I,J) = A(N+1-J,N+1-I.
            //
            //    A is positive definite.
            //
            //    A is an M matrix.
            //
            //    A is weakly diagonally dominant, but not strictly diagonally dominant.
            //
            //    A has an LU factorization A = L * U, without pivoting.
            //
            //      The matrix L is lower bidiagonal with subdiagonal elements:
            //
            //        L(I+1,I) = -I/(I+1)
            //
            //      The matrix U is upper bidiagonal, with diagonal elements
            //
            //        U(I,I) = (I+1)/I
            //
            //      and superdiagonal elements which are all -1.
            //
            //    A has a Cholesky factorization A = L * L', with L lower bidiagonal.
            //
            //      L(I,I) =    sqrt ( (I+1) / I )
            //      L(I,I-1) = -sqrt ( (I-1) / I )
            //
            //    The eigenvalues are
            //
            //      LAMBDA(I) = 2 + 2 * COS(I*PI/(N+1))
            //                = 4 SIN^2(I*PI/(2*N+2))
            //
            //    The corresponding eigenvector X(I) has entries
            //
            //       X(I)(J) = sqrt(2/(N+1)) * sin ( I*J*PI/(N+1) ).
            //
            //    Simple linear systems:
            //
            //      x = (1,1,1,...,1,1),   A*x=(1,0,0,...,0,1)
            //
            //      x = (1,2,3,...,n-1,n), A*x=(0,0,0,...,0,n+1)
            //
            //    det ( A ) = N + 1.
            //
            //    The value of the determinant can be seen by induction,
            //    and expanding the determinant across the first row:
            //
            //      det ( A(N) ) = 2 * det ( A(N-1) ) - (-1) * (-1) * det ( A(N-2) )
            //                = 2 * N - (N-1)
            //                = N + 1
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 June 2014
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
            //    Morris Newman, John Todd,
            //    Example A8,
            //    The evaluation of matrix inversion programs,
            //    Journal of the Society for Industrial and Applied Mathematics,
            //    Volume 6, Number 4, pages 466-476, 1958.
            //
            //    John Todd,
            //    Basic Numerical Mathematics,
            //    Volume 2: Numerical Algebra,
            //    Birkhauser, 1980,
            //    ISBN: 0817608117,
            //    LC: QA297.T58.
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
            //    Input, int M, N, the number of rows and columns.
            //
            //    Input, int MU, the number of superdiagonals.
            //    MU must be at least 0, and no more than N-1.
            //
            //    Output, double R8PBU_DIF2[(MU+1)*N], the matrix.
            //
        {
            double[] a = new double[(mu + 1) * n];

            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i < mu + 1; i++)
                {
                    a[i + j * (mu + 1)] = 0.0;
                }
            }

            for (int j = 1; j < n; j++)
            {
                int i = mu - 1;
                a[i + j * (mu + 1)] = -1.0;
            }

            for (int j = 0; j < n; j++)
            {
                int i = mu;
                a[i + j * (mu + 1)] = 2.0;
            }

            return a;
        }

        public static double[] r8pbu_mv(int m, int n, int mu, double[] a, double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PBU_MV multiplies a R8PBU matrix times a vector.
        //
        //  Discussion:
        //
        //    The R8PBU storage format is used for a symmetric positive definite band matrix.
        //
        //    To save storage, only the diagonal and upper triangle of A is stored,
        //    in a compact diagonal format that preserves columns.
        //
        //    The diagonal is stored in row MU+1 of the array.
        //    The first superdiagonal in row MU, columns 2 through N.
        //    The second superdiagonal in row MU-1, columns 3 through N.
        //    The MU-th superdiagonal in row 1, columns MU+1 through N.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 February 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, int MU, the number of superdiagonals in the matrix.
        //    MU must be at least 0 and no more than N-1.
        //
        //    Input, double A[(MU+1)*N], the matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double R8PBU_MV[M], the result vector A * x.
        //
        {
            double[] b = new double[m];
            //
            //  Multiply X by the diagonal of the matrix.
            //
            for (int j = 0; j < n; j++)
            {
                b[j] = a[mu + j * (mu + 1)] * x[j];
            }

            //
            //  Multiply X by the superdiagonals of the matrix.
            //
            for (int i = mu; 1 <= i; i--)
            {
                for (int j = mu + 2 - i; j <= n; j++)
                {
                    int ieqn = i + j - mu - 1;
                    b[ieqn - 1] = b[ieqn - 1] + a[i - 1 + (j - 1) * (mu + 1)] * x[j - 1];
                    b[j - 1] = b[j - 1] + a[i - 1 + (j - 1) * (mu + 1)] * x[ieqn - 1];
                }
            }

            return b;
        }

        public static double[] r8pbu_res(int m, int n, int mu, double[] a, double[] x, double[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PBU_RES computes the residual R = B-A*X for R8PBU matrices.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrix.
        //    M must be positive.
        //
        //    Input, int N, the number of columns of the matrix.
        //    N must be positive.
        //
        //    Input, int MU, the number of superdiagonals in the matrix.
        //    MU must be at least 0 and no more than N-1.
        //
        //    Input, double A[(MU+1)*N], the matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Input, double B[M], the desired result A * x.
        //
        //    Output, double R8PBU_RES[M], the residual R = B - A * X.
        //
        {
            double[] r = r8pbu_mv(m, n, mu, a, x);
            for (int i = 0; i < m; i++)
            {
                r[i] = b[i] - r[i];
            }

            return r;
        }
    }
}