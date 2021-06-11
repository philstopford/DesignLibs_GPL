using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static void r8ge_cg(int n, double[] a, double[] b, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8GE_CG uses the conjugate gradient method on an R8GE system.
            //
            //  Discussion:
            //
            //    The R8GE storage format is used for a general M by N matrix.  A storage 
            //    space is made for each entry.  The two dimensional logical
            //    array can be thought of as a vector of M*N entries, starting with
            //    the M entries in the column 1, then the M entries in column 2
            //    and so on.  Considered as a vector, the entry A(I,J) is then stored
            //    in vector location I+(J-1)*M.
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
            //    Input, double A[N*N], the matrix.
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
            ap = r8ge_mv(n, n, a, x);

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
                ap = r8ge_mv(n, n, a, p);
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
        }

        public static double[] r8ge_dif2(int m, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8GE_DIF2 returns the DIF2 matrix in R8GE format.
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
            //    05 July 2000
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
            //    Output, double R8GE_DIF2[M*N], the matrix.
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
                        a[i + j * m] = -1.0;
                    }
                    else if (j == i)
                    {
                        a[i * j * m] = 2.0;
                    }
                    else if (j == i + 1)
                    {
                        a[i + j * m] = -1.0;
                    }
                    else
                    {
                        a[i + j * m] = 0.0;
                    }
                }
            }

            return a;
        }

        public static double[] r8ge_mv(int m, int n, double[] a, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8GE_MV multiplies an R8GE matrix by an R8VEC.
            //
            //  Discussion:
            //
            //    The R8GE storage format is used for a general M by N matrix.  A storage 
            //    space is made for each entry.  The two dimensional logical
            //    array can be thought of as a vector of M*N entries, starting with
            //    the M entries in the column 1, then the M entries in column 2
            //    and so on.  Considered as a vector, the entry A(I,J) is then stored
            //    in vector location I+(J-1)*M.
            //
            //    R8GE storage is used by LINPACK and LAPACK.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 July 2014
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
            //    Input, double A[M*N], the matrix.
            //
            //    Input, double X[N], the vector to be multiplied by A.
            //
            //    Output, double R8GE_MV[M], the product A * x.
            //
        {
            int i;
            int j;
            double[] b;

            b = new double[m];

            for (i = 0; i < m; i++)
            {
                b[i] = 0.0;
                for (j = 0; j < n; j++)
                {
                    b[i] = b[i] + a[i + j * m] * x[j];
                }
            }

            return b;
        }

        public static double[] r8ge_res(int m, int n, double[] a, double[] x, double[] b)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8GE_RES computes the residual R = B-A*X for R8GE matrices.
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
            //    Input, double A[M*N], the matrix.
            //
            //    Input, double X[N], the vector to be multiplied by A.
            //
            //    Input, double B[M], the desired result A * x.
            //
            //    Output, double R8GE_RES[M], the residual R = B - A * X.
            //
        {
            int i;
            double[] r;

            r = r8ge_mv(m, n, a, x);
            for (i = 0; i < m; i++)
            {
                r[i] = b[i] - r[i];
            }

            return r;
        }

        public static double[] r8ge_mtm(int n, double[] a, double[] b)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8GE_MTM computes C=A'*B for R8GE matrices.
            //
            //  Discussion:
            //
            //    The R8GE storage format is used for a "general" M by N matrix.  
            //    A physical storage space is made for each logical entry.  The two 
            //    dimensional logical array is mapped to a vector, in which storage is 
            //    by columns.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    04 August 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the matrices.
            //    N must be positive.
            //
            //    Input, double A[N*N], B[N*N], the factors.
            //
            //    Output, double C[N*N], the product.
            //
        {
            double[] c = new double[n * n];

            for (int j = 0; j < n; j++)
            {
                for (int i = 0; i < n; i++)
                {
                    c[i + j * n] = 0.0;
                    for (int k = 0; k < n; k++)
                    {
                        c[i + j * n] = c[i + j * n] + a[k + i * n] * b[k + j * n];
                    }
                }
            }

            return c;
        }

        public static void r8ge_print(int m, int n, double[] a, string title)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8GE_PRINT prints an R8GE matrix.
            //
            //  Discussion:
            //
            //    The R8GE storage format is used for a "general" M by N matrix.  
            //    A physical storage space is made for each logical entry.  The two 
            //    dimensional logical array is mapped to a vector, in which storage is 
            //    by columns.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    06 April 2006
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
            //    Input, double A[M*N], the R8GE matrix.
            //
            //    Input, string TITLE, a title.
            //
        {
            r8ge_print_some(m, n, a, 1, 1, m, n, title);
        }

        public static void r8ge_print_some(int m, int n, double[] a, int ilo, int jlo, int ihi,
                int jhi, string title)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8GE_PRINT_SOME prints some of an R8GE matrix.
            //
            //  Discussion:
            //
            //    The R8GE storage format is used for a "general" M by N matrix.  
            //    A physical storage space is made for each logical entry.  The two 
            //    dimensional logical array is mapped to a vector, in which storage is 
            //    by columns.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    06 April 2006
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
            //    Input, double A[M*N], the R8GE matrix.
            //
            //    Input, int ILO, JLO, IHI, JHI, designate the first row and
            //    column, and the last row and column to be printed.
            //
            //    Input, string TITLE, a title.
            //
        {
            int INCX = 5;

            int i2hi;
            int i2lo;
            int j2hi;
            int j2lo;

            Console.WriteLine("");
            Console.WriteLine(title + "");
            //
            //  Print the columns of the matrix, in strips of 5.
            //
            for (j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX)
            {
                j2hi = j2lo + INCX - 1;
                j2hi = Math.Min(j2hi, n);
                j2hi = Math.Min(j2hi, jhi);

                Console.WriteLine("");
                //
                //  For each column J in the current range...
                //
                //  Write the header.
                //
                string cout = "  Col:    ";
                for (int j = j2lo; j <= j2hi; j++)
                {
                    cout += j.ToString().PadLeft(7) + "       ";
                }

                Console.WriteLine(cout);
                Console.WriteLine("  Row");
                Console.WriteLine("  ---");
                //
                //  Determine the range of the rows in this strip.
                //
                i2lo = Math.Max(ilo, 1);
                i2hi = Math.Min(ihi, m);

                for (int i = i2lo; i <= i2hi; i++)
                {
                    //
                    //  Print out (up to) 5 entries in row I, that lie in the current strip.
                    //
                    cout = i.ToString().PadLeft(5) + "  ";
                    for (int j = j2lo; j <= j2hi; j++)
                    {
                        cout += a[i - 1 + (j - 1) * m].ToString().PadLeft(12) + "  ";
                    }

                    Console.WriteLine(cout);
                }
            }
        }
    }
}