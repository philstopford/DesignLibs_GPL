using System;
using Burkardt.Uniform;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r83_print(int n, double[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83_PRINT prints a R83 matrix.
        //
        //  Discussion:
        //
        //    The R83 storage format is used for a tridiagonal matrix.
        //    The superdiagonal is stored in entries (1,2:N), the diagonal in
        //    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
        //    original matrix is "collapsed" vertically into the array.
        //
        //  Example:
        //
        //    Here is how a R83 matrix of order 5 would be stored:
        //
        //       *  A12 A23 A34 A45
        //      A11 A22 A33 A44 A55
        //      A21 A32 A43 A54  *
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
        //    Input, int N, the order of the matrix.
        //    N must be positive.
        //
        //    Input, double A[3*N], the R83 matrix.
        //
        //    Input, string TITLE, a title.
        //
    {
        r83_print_some(n, a, 1, 1, n, n, title);
    }

    public static void r83_print_some(int n, double[] a, int ilo, int jlo, int ihi, int jhi,
            string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83_PRINT_SOME prints some of a R83 matrix.
        //
        //  Discussion:
        //
        //    The R83 storage format is used for a tridiagonal matrix.
        //    The superdiagonal is stored in entries (1,2:N), the diagonal in
        //    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
        //    original matrix is "collapsed" vertically into the array.
        //
        //  Example:
        //
        //    Here is how a R83 matrix of order 5 would be stored:
        //
        //       *  A12 A23 A34 A45
        //      A11 A22 A33 A44 A55
        //      A21 A32 A43 A54  *
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
        //    Input, int N, the order of the matrix.
        //    N must be positive.
        //
        //    Input, double A[3*N], the R83 matrix.
        //
        //    Input, int ILO, JLO, IHI, JHI, designate the first row and
        //    column, and the last row and column, to be printed.
        //
        //    Input, string TITLE, a title.
        //
    {
        const int INCX = 5;

        int j2lo;

        Console.WriteLine("");
        Console.WriteLine(title + "");
        //
        //  Print the columns of the matrix, in strips of 5.
        //
        for (j2lo = jlo; j2lo <= jhi; j2lo += INCX)
        {
            int j2hi = j2lo + INCX - 1;
            j2hi = Math.Min(j2hi, n);
            j2hi = Math.Min(j2hi, jhi);

            int inc = j2hi + 1 - j2lo;

            Console.WriteLine("");
            string cout = "  Col: ";
            int j2;
            int j;
            for (j = j2lo; j <= j2hi; j++)
            {
                j2 = j + 1 - j2lo;
                cout += j.ToString().PadLeft(7) + "       ";
            }

            Console.WriteLine(cout);
            Console.WriteLine("  Row");
            Console.WriteLine("  ---");
            //
            //  Determine the range of the rows in this strip.
            //
            int i2lo = Math.Max(ilo, 1);
            i2lo = Math.Max(i2lo, j2lo - 1);

            int i2hi = Math.Min(ihi, n);
            i2hi = Math.Min(i2hi, j2hi + 1);

            int i;
            for (i = i2lo; i <= i2hi; i++)
            {
                //
                //  Print out (up to) 5 entries in row I, that lie in the current strip.
                //
                string cout2 = i.ToString().PadLeft(6) + ": ";

                for (j2 = 1; j2 <= inc; j2++)
                {
                    j = j2lo - 1 + j2;

                    if (1 < i - j || 1 < j - i)
                    {
                        cout2 += "              ";
                    }
                    else if (j == i + 1)
                    {
                        cout2 += a[0 + (j - 1) * 3].ToString().PadLeft(12) + "  ";
                    }
                    else if (j == i)
                    {
                        cout2 += a[1 + (j - 1) * 3].ToString().PadLeft(12) + "  ";
                    }
                    else if (j == i - 1)
                    {
                        cout2 += a[2 + (j - 1) * 3].ToString().PadLeft(12) + "  ";
                    }

                }

                Console.WriteLine(cout2);
            }
        }
    }

    public static void r83_cg(int n, double[] a, double[] b, ref double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83_CG uses the conjugate gradient method on an R83 system.
        //
        //  Discussion:
        //
        //    The R83 storage format is used for a tridiagonal matrix.
        //    The superdiagonal is stored in entries (1,2:N), the diagonal in
        //    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
        //    original matrix is "collapsed" vertically into the array.
        //
        //    The matrix A must be a positive definite symmetric band matrix.
        //
        //    The method is designed to reach the solution after N computational
        //    steps.  However, roundoff may introduce unacceptably large errors for
        //    some problems.  In such a case, calling the routine again, using
        //    the computed solution as the new starting estimate, should improve
        //    the results.
        //
        //  Example:
        //
        //    Here is how an R83 matrix of order 5 would be stored:
        //
        //       *  A12 A23 A34 A45
        //      A11 A22 A33 A44 A55
        //      A21 A32 A43 A54  *
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 June 2014
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
        //    Input, double A[3*N], the matrix.
        //
        //    Input, double B[N], the right hand side vector.
        //
        //    Input/output, double X[N].
        //    On input, an estimate for the solution, which may be 0.
        //    On output, the approximate solution vector.
        //
    {
        int i;
        int it;
        //
        //  Initialize
        //    AP = A * x,
        //    R  = b - A * x,
        //    P  = b - A * x.
        //
        double[] ap = r83_mv(n, n, a, x);

        double[] r = new double[n];
        for (i = 0; i < n; i++)
        {
            r[i] = b[i] - ap[i];
        }

        double[] p = new double[n];
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
            ap = r83_mv(n, n, a, p);
            //
            //  Compute the dot products
            //    PAP = P*AP,
            //    PR  = P*R
            //  Set
            //    ALPHA = PR / PAP.
            //
            double pap = r8vec_dot_product(n, p, ap);
            double pr = r8vec_dot_product(n, p, r);

            if (pap == 0.0)
            {
                break;
            }

            double alpha = pr / pap;
            //
            //  Set
            //    X = X + ALPHA * P
            //    R = R - ALPHA * AP.
            //
            for (i = 0; i < n; i++)
            {
                x[i] += alpha * p[i];
            }

            for (i = 0; i < n; i++)
            {
                r[i] -= alpha * ap[i];
            }

            //
            //  Compute the vector dot product
            //    RAP = R*AP
            //  Set
            //    BETA = - RAP / PAP.
            //
            double rap = r8vec_dot_product(n, r, ap);

            double beta = -rap / pap;
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
        
    public static double[] r83_dif2(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83_DIF2 returns the DIF2 matrix in R83 format.
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
        //    04 June 2014
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
        //    Output, double A[3*N], the matrix.
        //
    {
        int j;

        double[] a = new double[3 * n];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < 3; i++)
            {
                a[i + j * 3] = 0.0;
            }
        }

        int mn = Math.Min(m, n);

        for (j = 1; j < mn; j++)
        {
            a[0 + j * 3] = -1.0;
        }

        for (j = 0; j < mn; j++)
        {
            a[1 + j * 3] = 2.0;
        }

        for (j = 0; j < mn - 1; j++)
        {
            a[2 + j * 3] = -1.0;
        }

        if (n < m)
        {
            a[2 + (mn - 1) * 3] = -1.0;
        }

        return a;
    }

    public static double[] r83_indicator(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83_INDICATOR sets up an R83 indicator matrix.
        //
        //  Discussion:
        //
        //    The R83 storage format is used for a tridiagonal matrix.
        //    The superdiagonal is stored in entries (1,2:min(M+1,N)).
        //    The diagonal in entries (2,1:min(M,N)).
        //    The subdiagonal in (3,min(M-1,N)).
        //    R8GE A(I,J) = R83 A[I-J+1+J*3] (0 based indexing).
        //
        //  Example:
        //
        //    An R83 matrix of order 3x5 would be stored:
        //
        //       *  A12 A23 A34  *
        //      A11 A22 A33  *   *
        //      A21 A32  *   *   *
        //
        //    An R83 matrix of order 5x5 would be stored:
        //
        //       *  A12 A23 A34 A45
        //      A11 A22 A33 A44 A55
        //      A21 A32 A43 A54  *
        //
        //    An R83 matrix of order 5x3 would be stored:
        //
        //       *  A12 A23
        //      A11 A22 A33
        //      A21 A32 A43
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 September 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Output, double R83_INDICATOR[3*N], the R83 indicator matrix.
        //
    {
        int j;

        int fac = (int) Math.Pow(10, (int) Math.Log10(n) + 1);

        double[] a = r8ge_zeros_new(3, n);

        for (j = 0; j < n; j++)
        {
            int i_lo = Math.Max(0, j - 1);
            int i_hi = Math.Min(m - 1, j + 1);
            int i;
            for (i = i_lo; i <= i_hi; i++)
            {
                a[i - j + 1 + j * 3] = fac * (i + 1) + j + 1;
            }
        }

        return a;
    }

    public static void r83_jac_sl(int n, double[] a, double[] b, ref double[] x, int it_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83_JAC_SL solves an R83 system using Jacobi iteration.
        //
        //  Discussion:
        //
        //    The R83 storage format is used for a tridiagonal matrix.
        //    The superdiagonal is stored in entries (1,2:min(M+1,N)).
        //    The diagonal in entries (2,1:min(M,N)).
        //    The subdiagonal in (3,min(M-1,N)).
        //    R8GE A(I,J) = R83 A[I-J+1+J*3] (0 based indexing).
        //
        //    This routine simply applies a given number of steps of the
        //    iteration to an input approximate solution.  On first call, you can
        //    simply pass in the zero vector as an approximate solution.  If
        //    the returned value is not acceptable, you may call again, using
        //    it as the starting point for additional iterations.
        //
        //  Example:
        //
        //    An R83 matrix of order 3x5 would be stored:
        //
        //       *  A12 A23 A34  *
        //      A11 A22 A33  *   *
        //      A21 A32  *   *   *
        //
        //    An R83 matrix of order 5x5 would be stored:
        //
        //       *  A12 A23 A34 A45
        //      A11 A22 A33 A44 A55
        //      A21 A32 A43 A54  *
        //
        //    An R83 matrix of order 5x3 would be stored:
        //
        //       *  A12 A23
        //      A11 A22 A33
        //      A21 A32 A43
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 September 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be at least 2.
        //
        //    Input, double A[3*N], the R83 matrix.
        //
        //    Input, double B[N], the right hand side of the linear system.
        //
        //    Input/output, double X[N], an approximate solution to the system.
        //
        //    Input, int IT_MAX, the maximum number of iterations to take.
        //
    {
        int i;
        int it_num;

        double[] xnew = new double[n];
        //
        //  No diagonal matrix entry can be zero.
        //
        for (i = 0; i < n; i++)
        {
            switch (a[1 + i * 3])
            {
                case 0.0:
                    Console.WriteLine("");
                    Console.WriteLine("R83_JAC_SL - Fatal error!");
                    Console.WriteLine("  Zero diagonal entry, index = " + i + "");
                    return;
            }
        }

        for (it_num = 1; it_num <= it_max; it_num++)
        {
            //
            //  Solve A*x=b:
            //
            xnew[0] = b[0] - a[0 + 1 * 3] * x[1];
            for (i = 1; i < n - 1; i++)
            {
                xnew[i] = b[i] - a[2 + (i - 1) * 3] * x[i - 1] - a[0 + (i + 1) * 3] * x[i + 1];
            }

            xnew[n - 1] = b[n - 1] - a[2 + (n - 2) * 3] * x[n - 2];
            //
            //  Divide by the diagonal term, and overwrite X.
            //
            for (i = 0; i < n; i++)
            {
                x[i] = xnew[i] / a[1 + i * 3];
            }
        }

    }

    public static double[] r83_mtv(int m, int n, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83_MTV multiplies a vector times an R83 matrix.
        //
        //  Discussion:
        //
        //    The R83 storage format is used for a tridiagonal matrix.
        //    The superdiagonal is stored in entries (1,2:min(M+1,N)).
        //    The diagonal in entries (2,1:min(M,N)).
        //    The subdiagonal in (3,min(M-1,N)).
        //    R8GE A(I,J) = R83 A[I-J+1+J*3] (0 based indexing).
        //
        //  Example:
        //
        //    An R83 matrix of order 3x5 would be stored:
        //
        //       *  A12 A23 A34  *
        //      A11 A22 A33  *   *
        //      A21 A32  *   *   *
        //
        //    An R83 matrix of order 5x5 would be stored:
        //
        //       *  A12 A23 A34 A45
        //      A11 A22 A33 A44 A55
        //      A21 A32 A43 A54  *
        //
        //    An R83 matrix of order 5x3 would be stored:
        //
        //       *  A12 A23
        //      A11 A22 A33
        //      A21 A32 A43
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 September 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the linear system.
        //
        //    Input, double A[3*N], the R83 matrix.
        //
        //    Input, double X[M], the vector to be multiplied by A'.
        //
        //    Output, double R83_MTV[N], the product A' * x.
        //
    {
        int i;
        int j;

        double[] b = new double[n];

        for (i = 0; i < n; i++)
        {
            b[i] = 0.0;
        }

        for (j = 0; j < n; j++)
        {
            int i_lo = Math.Max(0, j - 1);
            int i_hi = Math.Min(m - 1, j + 1);
            for (i = i_lo; i <= i_hi; i++)
            {
                b[j] += x[i] * a[i - j + 1 + j * 3];
            }
        }

        return b;
    }

    public static double[] r83_mv(int m, int n, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83_MV multiplies a R83 matrix times a vector.
        //
        //  Discussion:
        //
        //    The R83 storage format is used for a tridiagonal matrix.
        //    The superdiagonal is stored in entries (1,2:N), the diagonal in
        //    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
        //    original matrix is "collapsed" vertically into the array.
        //
        //  Example:
        //
        //    Here is how a R83 matrix of order 5 would be stored:
        //
        //       *  A12 A23 A34 A45
        //      A11 A22 A33 A44 A55
        //      A21 A32 A43 A54  *
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A[3*N], the R83 matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double R83_MV[M], the product A * x.
        //
    {
        int i;

        double[] b = new double[m];

        for (i = 0; i < m; i++)
        {
            b[i] = 0.0;
        }

        int mn = Math.Min(m, n);

        for (i = 0; i < mn; i++)
        {
            b[i] += a[1 + i * 3] * x[i];
        }

        for (i = 0; i < mn - 1; i++)
        {
            b[i] += a[0 + (i + 1) * 3] * x[i + 1];
        }

        for (i = 1; i < mn; i++)
        {
            b[i] += a[2 + (i - 1) * 3] * x[i - 1];
        }

        if (n < m)
        {
            b[n] += a[2 + (n - 1) * 3] * x[n - 1];
        }
        else if (m < n)
        {
            b[m - 1] += a[0 + m * 3] * x[m];
        }

        return b;
    }

    public static double[] r83_random(int m, int n, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83_RANDOM randomizes an R83 matrix.
        //
        //  Discussion:
        //
        //    The R83 storage format is used for a tridiagonal matrix.
        //    The superdiagonal is stored in entries (1,2:min(M+1,N)).
        //    The diagonal in entries (2,1:min(M,N)).
        //    The subdiagonal in (3,min(M-1,N)).
        //    R8GE A(I,J) = R83 A[I-J+1+J*3] (0 based indexing).
        //
        //  Example:
        //
        //    An R83 matrix of order 3x5 would be stored:
        //
        //       *  A12 A23 A34  *
        //      A11 A22 A33  *   *
        //      A21 A32  *   *   *
        //
        //    An R83 matrix of order 5x5 would be stored:
        //
        //       *  A12 A23 A34 A45
        //      A11 A22 A33 A44 A55
        //      A21 A32 A43 A54  *
        //
        //    An R83 matrix of order 5x3 would be stored:
        //
        //       *  A12 A23
        //      A11 A22 A33
        //      A21 A32 A43
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 September 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the linear system.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double R83_RANDOM[3*N], the R83 matrix.
        //
    {
        int j;

        double[] a = r8ge_zeros_new(3, n);

        for (j = 0; j < n; j++)
        {
            int i_lo = Math.Max(0, j - 1);
            int i_hi = Math.Min(m - 1, j + 1);
            int i;
            for (i = i_lo; i <= i_hi; i++)
            {
                a[i - j + 1 + j * 3] = UniformRNG.r8_uniform_01(ref seed);
            }
        }

        return a;
    }

    public static double[] r83_res(int m, int n, double[] a, double[] x, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83_RES computes the residual R = B-A*X for R83 matrices.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 June 2014
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
        //    Input, double A[3*N], the matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Input, double B[M], the desired result A * x.
        //
        //    Output, double R83_RES[M], the residual R = B - A * X.
        //
    {
        int i;

        double[] r = r83_mv(m, n, a, x);

        for (i = 0; i < m; i++)
        {
            r[i] = b[i] - r[i];
        }

        return r;
    }

    public static double r83_norm(double x, double y, double z)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83_NORM returns the Euclidean norm of an R83.
        //
        //  Discussion:
        //
        //    An R83 is a vector of 3 R8's.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 November 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, Y, Z, the vector.
        //
        //    Output, double R83_NORM, the norm of the vector.
        //
    {
        double value = Math.Sqrt(x * x + y * y + z * z);

        return value;
    }

    public static double[] r83_to_r8ge(int m, int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83_TO_R8GE copies an R83 matrix to an R8GE matrix.
        //
        //  Discussion:
        //
        //    The R83 storage format is used for a tridiagonal matrix.
        //    The superdiagonal is stored in entries (1,2:min(M+1,N)).
        //    The diagonal in entries (2,1:min(M,N)).
        //    The subdiagonal in (3,min(M-1,N)).
        //    R8GE A(I,J) = R83 A[I-J+1+J*3] (0 based indexing).
        //
        //  Example:
        //
        //    An R83 matrix of order 3x5 would be stored:
        //
        //       *  A12 A23 A34  *
        //      A11 A22 A33  *   *
        //      A21 A32  *   *   *
        //
        //    An R83 matrix of order 5x5 would be stored:
        //
        //       *  A12 A23 A34 A45
        //      A11 A22 A33 A44 A55
        //      A21 A32 A43 A54  *
        //
        //    An R83 matrix of order 5x3 would be stored:
        //
        //       *  A12 A23
        //      A11 A22 A33
        //      A21 A32 A43
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 September 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Input, double A[3*N], the R83 matrix.
        //
        //    Output, double R83_TO_R8GE[M*N], the R8GE matrix.
        //
    {
        int j;

        double[] b = r8ge_zeros_new(m, n);

        for (j = 0; j < n; j++)
        {
            int i_lo = Math.Max(0, j - 1);
            int i_hi = Math.Min(m - 1, j + 1);
            int i;
            for (i = i_lo; i <= i_hi; i++)
            {
                b[i + j * m] = a[i - j + 1 + j * 3];
            }
        }

        return b;
    }

    public static double[] r83_zeros(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83_ZEROS zeros an R83 matrix.
        //
        //  Discussion:
        //
        //    The R83 storage format is used for a tridiagonal matrix.
        //    The superdiagonal is stored in entries (1,2:min(M+1,N)).
        //    The diagonal in entries (2,1:min(M,N)).
        //    The subdiagonal in (3,min(M-1,N)).
        //    R8GE A(I,J) = R83 A[I-J+1+J*3] (0 based indexing).
        //
        //  Example:
        //
        //    An R83 matrix of order 3x5 would be stored:
        //
        //       *  A12 A23 A34  *
        //      A11 A22 A33  *   *
        //      A21 A32  *   *   *
        //
        //    An R83 matrix of order 5x5 would be stored:
        //
        //       *  A12 A23 A34 A45
        //      A11 A22 A33 A44 A55
        //      A21 A32 A43 A54  *
        //
        //    An R83 matrix of order 5x3 would be stored:
        //
        //       *  A12 A23
        //      A11 A22 A33
        //      A21 A32 A43
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 September 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Output, double R83_ZERO[3*N], the R83 matrix.
        //
    {
        double[] a = r8ge_zeros_new(3, n);

        return a;
    }

}