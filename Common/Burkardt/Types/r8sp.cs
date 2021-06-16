using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static void r8sp_cg(int n, int nz_num, int[] row, int[] col, double[] a,
        double[] b, ref double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SP_CG uses the conjugate gradient method on a R8SP linear system.
        //
        //  Discussion:
        //
        //    The R8SP storage format stores the row, column and value of each nonzero
        //
        //    It is possible that a pair of indices (I,J) may occur more than
        //    once.  Presumably, in this case, the intent is that the actual value
        //    of A(I,J) is the sum of all such entries.  This is not a good thing
        //    to do, but I seem to have come across this in MATLAB.
        //
        //    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP
        //    (nonsymmetric case), by MATLAB, and by SPARSEKIT ("COO" format).
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
        //    05 June 2014
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
        //    Input, int NZ_NUM, the number of nonzero elements in the matrix.
        //
        //    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and column indices
        //    of the nonzero elements.
        //
        //    Input, double A[NZ_NUM], the nonzero elements of the matrix.
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
            ap = r8sp_mv(n, n, nz_num, row, col, a, x);

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
                //  Compute the matrix*vector product AP = A*P.
                //
                ap = r8sp_mv(n, n, nz_num, row, col, a, p);
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
            return;
        }

        public static double[] r8sp_dif2(int m, int n, int nz_num, int[] row, int[] col )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SP_DIF2 returns the DIF2 matrix in R8SP format.
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
        //    Input, int NZ_NUM, the number of nonzeros.
        //
        //    Input, int ROW[NZ_NUM], COL[NZ_NUM], space in which the rows and columns
        //    of nonzero entries will be stored.
        //
        //    Output, double R8SP_DIF2[NZ_NUM], the matrix.
        //
        {
            double[] a = new double[nz_num];

            for (int i = 0; i < nz_num; i++)
            {
                row[i] = 0;
                col[i] = 0;
                a[i] = 0.0;
            }

            int mn = Math.Min(m, n);

            int k = 0;
            for (int i = 0; i < mn; i++)
            {
                if (0 < i)
                {
                    row[k] = i;
                    col[k] = i - 1;
                    a[k] = -1.0;
                    k = k + 1;
                }

                row[k] = i;
                col[k] = i;
                a[k] = 2.0;
                k = k + 1;

                if (i < n - 1)
                {
                    row[k] = i;
                    col[k] = i + 1;
                    a[k] = -1.0;
                    k = k + 1;
                }
            }

            return a;
        }

        public static double[] r8sp_mv(int m, int n, int nz_num, int[] row, int[] col,
        double[] a, double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SP_MV multiplies a R8SP matrix times a vector.
        //
        //  Discussion:
        //
        //    The R8SP storage format stores the row, column and value of each nonzero
        //
        //    It is possible that a pair of indices (I,J) may occur more than
        //    once.  Presumably, in this case, the intent is that the actual value
        //    of A(I,J) is the sum of all such entries.  This is not a good thing
        //    to do, but I seem to have come across this in MATLAB.
        //
        //    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP
        //    (nonsymmetric case), by MATLAB, and by SPARSEKIT ("COO" format).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns of the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero elements in the matrix.
        //
        //    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and column indices
        //    of the nonzero elements.
        //
        //    Input, double A[NZ_NUM], the nonzero elements of the matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double R8SP_MV[M], the product vector A*X.
        //
        {
            double[] b = new double[m];

            for (int i = 0; i < m; i++)
            {
                b[i] = 0.0;
            }

            for (int k = 0; k < nz_num; k++)
            {
                int i = row[k];
                int j = col[k];
                b[i] = b[i] + a[k] * x[j];
            }

            return b;
        }

        public static double[] r8sp_res(int m, int n, int nz_num, int[] row, int[] col, double[] a,
        double[] x, double[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SP_RES computes the residual R = B-A*X for R8SP matrices.
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
        //    Input, int NZ_NUM, the number of nonzeros.
        //
        //    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and column indices.
        //
        //    Input, double A[NZ_NUM], the values.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Input, double B[M], the desired result A * x.
        //
        //    Output, double R8SP_RES[M], the residual R = B - A * X.
        //
        {
            double[] r = r8sp_mv(m, n, nz_num, row, col, a, x);

            for (int i = 0; i < m; i++)
            {
                r[i] = b[i] - r[i];
            }

            return r;
        }
    }
}