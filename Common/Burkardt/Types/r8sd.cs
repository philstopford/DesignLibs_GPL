namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static void r8sd_cg(int n, int ndiag, int[] offset, double[] a, double[] b,
        ref double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SD_CG uses the conjugate gradient method on a R8SD linear system.
        //
        //  Discussion:
        //
        //    The R8SD storage format is used for symmetric matrices whose only nonzero entries
        //    occur along a few diagonals, but for which these diagonals are not all
        //    close enough to the main diagonal for band storage to be efficient.
        //
        //    In that case, we assign the main diagonal the offset value 0, and 
        //    each successive superdiagonal gets an offset value 1 higher, until
        //    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
        //
        //    Assuming there are NDIAG nonzero diagonals (ignoring subdiagonals!),
        //    we then create an array B that has N rows and NDIAG columns, and simply
        //    "collapse" the matrix A to the left:
        //
        //    For the conjugate gradient method to be applicable, the matrix A must 
        //    be a positive definite symmetric matrix.
        //
        //    The method is designed to reach the solution to the linear system
        //      A * x = b
        //    after N computational steps.  However, roundoff may introduce
        //    unacceptably large errors for some problems.  In such a case,
        //    calling the routine a second time, using the current solution estimate
        //    as the new starting guess, should result in improved results.
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
        //    Input, int NDIAG, the number of diagonals that are stored.
        //    NDIAG must be at least 1 and no more than N.
        //
        //    Input, int OFFSET[NDIAG], the offsets for the diagonal storage.
        //
        //    Input, double A[N*NDIAG], the matrix.
        //
        //    Input, double B[N], the right hand side vector.
        //
        //    Input/output, double X[N].
        //    On input, an estimate for the solution, which may be 0.
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
            ap = r8sd_mv(n, n, ndiag, offset, a, x);

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

            for (it = 1; it < n; it++)
            {
                //
                //  Compute the matrix*vector product AP = A*P.
                //
                ap = r8sd_mv(n, n, ndiag, offset, a, p);
                //
                //  Compute the dot products
                //    PAP = P*AP,
                //    PR  = P*R
                //  Set
                //    ALPHA = PR / PAP.
                //
                pap = r8vec_dot_product(n, p, ap);

                if (pap == 0.0)
                {
                    break;
                }

                pr = r8vec_dot_product(n, p, r);

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
                rap = r8vec_dot_product(n, r, ap);

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

        public static double[] r8sd_dif2(int m, int n, int ndiag, int[] offset)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8SD_DIF2 returns the DIF2 matrix in R8SD format.
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
            //    Input, int NDIAG, the number of diagonals available for storage.
            //
            //    Input, int OFFSET[NDIAG], the indices of the diagonals.  It is
            //    presumed that OFFSET[0] = 0 and OFFSET[1] = 1.
            //
            //    Output, double R8SD_DIF2[N*NDIAG], the matrix.
            //
        {
            double[] a = new double[n * ndiag];

            for (int j = 0; j < ndiag; j++)
            {
                for (int i = 0; i < n; i++)
                {
                    a[i + j * n] = 0.0;
                }
            }

            for (int i = 0; i < n; i++)
            {
                int j = 0;
                a[i + j * ndiag] = 2.0;
            }

            for (int i = 0; i < n - 1; i++)
            {
                int j = 1;
                a[i + j * ndiag] = -1.0;
            }

            return a;
        }

        public static double[] r8sd_mv(int m, int n, int ndiag, int[] offset, double[] a,
        double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SD_MV multiplies a R8SD matrix times a vector.
        //
        //  Discussion:
        //
        //    The R8SD storage format is used for symmetric matrices whose only nonzero entries
        //    occur along a few diagonals, but for which these diagonals are not all
        //    close enough to the main diagonal for band storage to be efficient.
        //
        //    In that case, we assign the main diagonal the offset value 0, and 
        //    each successive superdiagonal gets an offset value 1 higher, until
        //    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
        //
        //    Assuming there are NDIAG nonzero diagonals (ignoring subdiagonals!),
        //    we then create an array B that has N rows and NDIAG columns, and simply
        //    "collapse" the matrix A to the left.
        //
        //  Example:
        //
        //    The "offset" value is printed above each column.
        //
        //    Original matrix               New Matrix
        //
        //       0   1   2   3   4   5       0   1   3   5
        //
        //      11  12   0  14   0  16      11  12  14  16
        //      21  22  23   0  25   0      22  23  25  --
        //       0  32  33  34   0  36      33  34  36  --
        //      41   0  43  44  45   0      44  45  --  --
        //       0  52   0  54  55  56      55  56  --  --
        //      61   0  63   0  65  66      66  --  --  --
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
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, int NDIAG, the number of diagonals that are stored.
        //    NDIAG must be at least 1 and no more than N.
        //
        //    Input, int OFFSET[NDIAG], the offsets for the diagonal storage.
        //
        //    Input, double A[N*NDIAG], the matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double R8SD_MV[N], the product A * x.
        //
        {
            double[] b = new double[m];

            for (int i = 0; i < m; i++)
            {
                b[i] = 0.0;
            }

            for (int i = 0; i < n; i++)
            {
                for (int jdiag = 0; jdiag < ndiag; jdiag++)
                {
                    int j = i + offset[jdiag];
                    if (0 <= j && j < n)
                    {
                        b[i] = b[i] + a[i + jdiag * n] * x[j];
                        if (offset[jdiag] != 0)
                        {
                            b[j] = b[j] + a[i + jdiag * n] * x[i];
                        }
                    }
                }
            }

            return b;
        }

        public static double[] r8sd_res(int m, int n, int ndiag, int[] offset, double[] a,
        double[] x, double[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SD_RES computes the residual R = B-A*X for R8SD matrices.
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
        //    Input, int NDIAG, the number of diagonals that are stored.
        //    NDIAG must be at least 1 and no more than N.
        //
        //    Input, int OFFSET[NDIAG], the offsets for the diagonal storage.
        //
        //    Input, double A[N*NDIAG], the matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Input, double B[M], the desired result A * x.
        //
        //    Output, double R8SD_RES[M], the residual R = B - A * X.
        //
        {
            double[] r = r8sd_mv(m, n, ndiag, offset, a, x);

            for (int i = 0; i < m; i++)
            {
                r[i] = b[i] - r[i];
            }

            return r;
        }
    }
}