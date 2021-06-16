using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static void r83t_cg(int n, double[] a, double[] b, ref double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83T_CG uses the conjugate gradient method on an R83T system.
        //
        //  Discussion:
        //
        //    The R83T storage format is used for a tridiagonal matrix.
        //    The superdiagonal is stored in entries (1:N-1,3), the diagonal in
        //    entries (1:N,2), and the subdiagonal in (2:N,1).  Thus, the
        //    original matrix is "collapsed" horizontally into the array.
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
        //    18 June 2014
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
        //    Input, double A[N*3], the matrix.
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
            ap = r83t_mv(n, n, a, x);

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
                ap = r83t_mv(n, n, a, p);
                //
                //  Compute the dot products
                //    PAP = P*AP,
                //    PR  = P*R
                //  Set
                //    ALPHA = PR / PAP.
                //
                pap = r8vec_dot_product(n, p, ap);
                pr = r8vec_dot_product(n, p, r);

                if (pap == 0.0)
                {
                    break;
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

        public static double[] r83t_dif2(int m, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R83T_DIF2 returns the DIF2 matrix in R83T format.
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
            //    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).
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
            //    18 June 2014
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
            //    Input, int M, N, the order of the matrix.
            //
            //    Output, double A[M*3], the matrix.
            //
        {
            double[] a;
            int i;
            int j;
            int mn;

            a = new double[m * 3];

            for (j = 0; j < 3; j++)
            {
                for (i = 0; i < m; i++)
                {
                    a[i + j * m] = 0.0;
                }
            }

            mn = Math.Min(m, n);

            j = 0;
            for (i = 1; i < mn; i++)
            {
                a[i + j * m] = -1.0;
            }

            j = 1;
            for (i = 0; i < mn; i++)
            {
                a[i + j * m] = 2.0;
            }

            j = 2;
            for (i = 0; i < mn - 1; i++)
            {
                a[i + j * m] = -1.0;
            }

            if (m < n)
            {
                i = mn - 1;
                j = 2;
                a[i + j * m] = -1.0;
            }
            else if (n < m)
            {
                i = mn;
                j = 0;
                a[i + j * m] = -1.0;
            }

            return a;
        }

        public static double[] r83t_mv(int m, int n, double[] a, double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83T_MV multiplies a R83T matrix times a vector.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A[M*3], the matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double R83T_MV[M], the product A * x.
        //
        {
            double[] b;
            int i;
            int j;
            int mn;

            b = new double[m];

            for (i = 0; i < m; i++)
            {
                b[i] = 0.0;
            }

            if (n == 1)
            {
                i = 0;
                j = 1;
                b[0] = a[i + j * m] * x[0];
                if (1 < m)
                {
                    i = 1;
                    j = 0;
                    b[1] = a[i + j * m] * x[0];
                }

                return b;
            }

            mn = Math.Min(m, n);

            b[0] = a[0 + 1 * m] * x[0]
                   + a[0 + 2 * m] * x[1];

            for (i = 1; i < mn - 1; i++)
            {
                b[i] = a[i + 0 * m] * x[i - 1]
                       + a[i + 1 * m] * x[i]
                       + a[i + 2 * m] * x[i + 1];
            }

            b[mn - 1] = a[mn - 1 + 0 * m] * x[mn - 2]
                        + a[mn - 1 + 1 * m] * x[mn - 1];

            if (n < m)
            {
                b[mn] = b[mn] + a[mn + 0 * m] * x[mn - 1];
            }
            else if (m < n)
            {
                b[mn - 1] = b[mn - 1] + a[mn - 1 + 2 * m] * x[mn];
            }

            return b;
        }

        public static double[] r83t_res(int m, int n, double[] a, double[] x, double[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83T_RES computes the residual R = B-A*X for R83T matrices.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 June 2014
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
        //    Input, double A[M*3], the matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Input, double B[M], the desired result A * x.
        //
        //    Output, double R83T_RES[M], the residual R = B - A * X.
        //
        {
            int i;
            double[] r;

            r = r83t_mv(m, n, a, x);

            for (i = 0; i < m; i++)
            {
                r[i] = b[i] - r[i];
            }

            return r;
        }
    }
}