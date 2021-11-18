using System;
using System.Globalization;
using Burkardt.Uniform;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r83t_cg(int n, double[] a, double[] b, ref double[] x)

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
        int i;
        int it;
        //
        //  Initialize
        //    AP = A * x,
        //    R  = b - A * x,
        //    P  = b - A * x.
        //
        double[] ap = r83t_mv(n, n, a, x);

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
            ap = r83t_mv(n, n, a, p);
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
        int i;
        int j;

        double[] a = new double[m * 3];

        for (j = 0; j < 3; j++)
        {
            for (i = 0; i < m; i++)
            {
                a[i + j * m] = 0.0;
            }
        }

        int mn = Math.Min(m, n);

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
            a[i] = -1.0;
        }

        return a;
    }

    public static void r83t_gs_sl(int n, double[] a, double[] b, ref double[] x, int it_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83T_GS_SL solves an R83T system using Gauss-Seidel iteration.
        //
        //  Discussion:
        //
        //    The R83T storage format is used for an MxN tridiagonal matrix.
        //    The superdiagonal is stored in entries (1:M-1,3), the diagonal in
        //    entries (1:M,2), and the subdiagonal in (2:M,1).  Thus, the
        //    the rows of the original matrix slide horizontally to form an
        //    Mx3 stack of data.
        //
        //    An R83T matrix of order 3x5 would be stored:
        //
        //       *  A11 A12
        //      A21 A22 A23
        //      A32 A33 A34
        //
        //    An R83T matrix of order 5x5 would be stored:
        //
        //       *  A11 A12
        //      A21 A22 A23
        //      A32 A33 A34
        //      A43 A44 A45
        //      A54 A55  *
        //
        //    An R83T matrix of order 5x3 would be stored:
        //
        //       *  A11 A12
        //      A21 A22 A23
        //      A32 A33  *
        //      A43  *   *
        //       *   *   *
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 May 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[N*3], the R83T matrix.
        //
        //    Input, double B[N], the right hand side of the linear system.
        //
        //    Input/output, double X[N], an approximate solution to 
        //    the system.
        //
        //    Input, int IT_MAX, the maximum number of iterations.
        //
    {
        int i;
        int it_num;
        //
        //  No diagonal matrix entry can be zero.
        //
        for (i = 0; i < n; i++)
        {
            switch (a[i + 1 * n])
            {
                case 0.0:
                    Console.WriteLine("");
                    Console.WriteLine("R83_GS_SL - Fatal error!");
                    Console.WriteLine("  Zero diagonal entry, index = " + i + "");
                    return;
            }
        }

        for (it_num = 1; it_num <= it_max; it_num++)
        {
            x[0] = (b[0] - a[0 + 2 * n] * x[1]) / a[0 + 1 * n];

            for (i = 1; i < n - 1; i++)
            {
                x[i] = (b[i] - a[i + 0 * n] * x[i - 1] - a[i + 2 * n] * x[i + 1]) / a[i + 1 * n];
            }

            x[n - 1] = (b[n - 1] - a[n - 1 + 0 * n] * x[n - 2]) / a[n - 1 + 1 * n];
        }

    }

    public static double[] r83t_indicator(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83T_INDICATOR sets the indicator matrix in R83T format.
        //
        //  Discussion:
        //
        //    The R83T storage format is used for an MxN tridiagonal matrix.
        //    The superdiagonal is stored in entries (1:M-1,3), the diagonal in
        //    entries (1:M,2), and the subdiagonal in (2:M,1).  Thus, the
        //    the rows of the original matrix slide horizontally to form an
        //    Mx3 stack of data.
        //
        //    An R83T matrix of order 3x5 would be stored:
        //
        //       *  A11 A12
        //      A21 A22 A23
        //      A32 A33 A34
        //
        //    An R83T matrix of order 5x5 would be stored:
        //
        //       *  A11 A12
        //      A21 A22 A23
        //      A32 A33 A34
        //      A43 A44 A45
        //      A54 A55  *
        //
        //    An R83T matrix of order 5x3 would be stored:
        //
        //       *  A11 A12
        //      A21 A22 A23
        //      A32 A33  *
        //      A43  *   *
        //       *   *   *
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 May 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrix.
        //
        //    Input, int N, the number of columns of the matrix.
        //
        //    Output, double R83T_INDICATOR[M*3], the matrix.
        //
    {
        int i;

        int fac = (int) Math.Pow(10, (int) Math.Log10(n) + 1);

        double[] a = new double[m * 3];

        for (i = 0; i < m; i++)
        {
            int k;
            for (k = 0; k < 3; k++)
            {
                int j = i + k - 1;
                a[i + k * m] = j switch
                {
                    >= 0 when j <= n - 1 => fac * (i + 1) + j + 1,
                    _ => 0.0
                };
            }
        }

        return a;
    }

    public static void r83t_jac_sl(int n, double[] a, double[] b, ref double[] x, int it_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83T_JAC_SL solves an R83T system using Jacobi iteration.
        //
        //  Discussion:
        //
        //    The R83T storage format is used for an MxN tridiagonal matrix.
        //    The superdiagonal is stored in entries (1:M-1,3), the diagonal in
        //    entries (1:M,2), and the subdiagonal in (2:M,1).  Thus, the
        //    the rows of the original matrix slide horizontally to form an
        //    Mx3 stack of data.
        //
        //    An R83T matrix of order 3x5 would be stored:
        //
        //       *  A11 A12
        //      A21 A22 A23
        //      A32 A33 A34
        //
        //    An R83T matrix of order 5x5 would be stored:
        //
        //       *  A11 A12
        //      A21 A22 A23
        //      A32 A33 A34
        //      A43 A44 A45
        //      A54 A55  *
        //
        //    An R83T matrix of order 5x3 would be stored:
        //
        //       *  A11 A12
        //      A21 A22 A23
        //      A32 A33  *
        //      A43  *   *
        //       *   *   *
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 May 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[N*3], the R83T matrix.
        //
        //    Input, double B[N], the right hand side of the linear system.
        //
        //    Input/output, double X[N], an approximate solution 
        //    to the system.
        //
        //    Input, int IT_MAX, the maximum number of iterations.
        //
    {
        int i;
        int it_num;
        //
        //  No diagonal matrix entry can be zero.
        //
        for (i = 0; i < n; i++)
        {
            switch (a[i + 1 * n])
            {
                case 0.0:
                    Console.WriteLine("");
                    Console.WriteLine("R83_JAC_SL - Fatal error!");
                    Console.WriteLine("  Zero diagonal entry, index = " + i + "");
                    return;
            }
        }

        double[] x_new = new double[n];

        for (it_num = 1; it_num <= it_max; it_num++)
        {
            x_new[0] = b[0] - a[0 + 2 * n] * x[1];
            for (i = 1; i < n - 1; i++)
            {
                x_new[i] = b[i] - a[i + 0 * n] * x[i - 1] - a[i + 2 * n] * x[i + 1];
            }

            x_new[n - 1] = b[n - 1] - a[n - 1 + 0 * n] * x[n - 2];
            //
            //  Divide by diagonal terms.
            //
            for (i = 0; i < n; i++)
            {
                x_new[i] /= a[i + 1 * n];
            }

            //
            //  Update.
            //
            for (i = 0; i < n; i++)
            {
                x[i] = x_new[i];
            }
        }
    }

    public static double[] r83t_mtv(int m, int n, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83T_MTV multiplies an R83T matrix transposed times an R8VEC.
        //
        //  Discussion:
        //
        //    The R83T storage format is used for an MxN tridiagonal matrix.
        //    The superdiagonal is stored in entries (1:M-1,3), the diagonal in
        //    entries (1:M,2), and the subdiagonal in (2:M,1).  Thus, the
        //    the rows of the original matrix slide horizontally to form an
        //    Mx3 stack of data.
        //
        //    An R83T matrix of order 3x5 would be stored:
        //
        //       *  A11 A12
        //      A21 A22 A23
        //      A32 A33 A34
        //
        //    An R83T matrix of order 5x5 would be stored:
        //
        //       *  A11 A12
        //      A21 A22 A23
        //      A32 A33 A34
        //      A43 A44 A45
        //      A54 A55  *
        //
        //    An R83T matrix of order 5x3 would be stored:
        //
        //       *  A11 A12
        //      A21 A22 A23
        //      A32 A33  *
        //      A43  *   *
        //       *   *   *
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 May 2016
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
        //    Input, double X[M], the vector to be multiplied by A.
        //
        //    Output, double R83T_MTV[N], the product A' * x.
        //
    {
        int i;
        int j;

        double[] b = new double[n];

        for (j = 0; j < n; j++)
        {
            b[j] = 0.0;
        }

        for (i = 0; i < m; i++)
        {
            int k;
            for (k = 0; k < 3; k++)
            {
                j = i + k - 1;
                switch (j)
                {
                    case >= 0 when j <= n - 1:
                        b[j] += x[i] * a[i + k * m];
                        break;
                }
            }
        }

        return b;
    }

    public static double[] r83t_mv(int m, int n, double[] a, double[] x)

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
        int i;

        double[] b = new double[m];

        for (i = 0; i < m; i++)
        {
            b[i] = 0.0;
        }

        switch (n)
        {
            case 1:
            {
                b[0] = a[m] * x[0];
                switch (m)
                {
                    case > 1:
                        i = 1;
                        b[1] = a[i] * x[0];
                        break;
                }

                return b;
            }
        }

        int mn = Math.Min(m, n);

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
            b[mn] += a[mn] * x[mn - 1];
        }
        else if (m < n)
        {
            b[mn - 1] += a[mn - 1 + 2 * m] * x[mn];
        }

        return b;
    }

    public static void r83t_print(int m, int n, double[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83T_PRINT prints an R83T matrix.
        //
        //  Discussion:
        //
        //    The R83T storage format is used for an MxN tridiagonal matrix.
        //    The superdiagonal is stored in entries (1:M-1,3), the diagonal in
        //    entries (1:M,2), and the subdiagonal in (2:M,1).  Thus, the
        //    the rows of the original matrix slide horizontally to form an
        //    Mx3 stack of data.
        //
        //    An R83T matrix of order 3x5 would be stored:
        //
        //       *  A11 A12
        //      A21 A22 A23
        //      A32 A33 A34
        //
        //    An R83T matrix of order 5x5 would be stored:
        //
        //       *  A11 A12
        //      A21 A22 A23
        //      A32 A33 A34
        //      A43 A44 A45
        //      A54 A55  *
        //
        //    An R83T matrix of order 5x3 would be stored:
        //
        //       *  A11 A12
        //      A21 A22 A23
        //      A32 A33  *
        //      A43  *   *
        //       *   *   *
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 May 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Input, double A[3*N], the R83T matrix.
        //
        //    Input, string TITLE, a title.
        //
    {
        r83t_print_some(m, n, a, 0, 0, m - 1, n - 1, title);
    }

    public static void r83t_print_some(int m, int n, double[] a, int ilo, int jlo, int ihi,
            int jhi, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83T_PRINT_SOME prints some of an R83T matrix.
        //
        //  Discussion:
        //
        //    The R83T storage format is used for an MxN tridiagonal matrix.
        //    The superdiagonal is stored in entries (1:M-1,3), the diagonal in
        //    entries (1:M,2), and the subdiagonal in (2:M,1).  Thus, the
        //    the rows of the original matrix slide horizontally to form an
        //    Mx3 stack of data.
        //
        //    An R83T matrix of order 3x5 would be stored:
        //
        //       *  A11 A12
        //      A21 A22 A23
        //      A32 A33 A34
        //
        //    An R83T matrix of order 5x5 would be stored:
        //
        //       *  A11 A12
        //      A21 A22 A23
        //      A32 A33 A34
        //      A43 A44 A45
        //      A54 A55  *
        //
        //    An R83T matrix of order 5x3 would be stored:
        //
        //       *  A11 A12
        //      A21 A22 A23
        //      A32 A33  *
        //      A43  *   *
        //       *   *   *
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 May 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Input, double A[3*N], the R83T matrix.
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
            j2hi = Math.Min(j2hi, n - 1);
            j2hi = Math.Min(j2hi, jhi);

            int inc = j2hi + 1 - j2lo;

            Console.WriteLine("");
            string cout = "  Col: ";
            int j2;
            int j;
            for (j = j2lo; j <= j2hi; j++)
            {
                j2 = j + 1 - j2lo;
                cout += j.ToString(CultureInfo.InvariantCulture).PadLeft(7) + "       ";
            }

            Console.WriteLine(cout);
            Console.WriteLine("  Row");
            Console.WriteLine("  ---");
            //
            //  Determine the range of the rows in this strip.
            //
            int i2lo = Math.Max(ilo, 0);
            i2lo = Math.Max(i2lo, j2lo - 1);

            int i2hi = Math.Min(ihi, m - 1);
            i2hi = Math.Min(i2hi, j2hi + 1);

            int i;
            for (i = i2lo; i <= i2hi; i++)
            {
                //
                //  Print out (up to) 5 entries in row I, that lie in the current strip.
                //
                cout = i.ToString(CultureInfo.InvariantCulture).PadLeft(6);

                for (j2 = 1; j2 <= inc; j2++)
                {
                    j = j2lo - 1 + j2;
                    int k = j - i + 1;
                    switch (k)
                    {
                        case < 0:
                        case > 2:
                            cout += "              ";
                            break;
                        default:
                            cout += "  " + a[i + k * m].ToString(CultureInfo.InvariantCulture).PadLeft(12);
                            break;
                    }
                }

                Console.WriteLine(cout);
            }
        }
    }

    public static double[] r83t_random(int m, int n, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83T_RANDOM returns a random R83T matrix.
        //
        //  Discussion:
        //
        //    The R83T storage format is used for an MxN tridiagonal matrix.
        //    The superdiagonal is stored in entries (1:M-1,3), the diagonal in
        //    entries (1:M,2), and the subdiagonal in (2:M,1).  Thus, the
        //    the rows of the original matrix slide horizontally to form an
        //    Mx3 stack of data.
        //
        //    An R83T matrix of order 3x5 would be stored:
        //
        //       *  A11 A12
        //      A21 A22 A23
        //      A32 A33 A34
        //
        //    An R83T matrix of order 5x5 would be stored:
        //
        //       *  A11 A12
        //      A21 A22 A23
        //      A32 A33 A34
        //      A43 A44 A45
        //      A54 A55  *
        //
        //    An R83T matrix of order 5x3 would be stored:
        //
        //       *  A11 A12
        //      A21 A22 A23
        //      A32 A33  *
        //      A43  *   *
        //       *   *   *
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 May 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrix.
        //
        //    Input, int N, the number of columns of the matrix.
        //
        //    Input/output, int *SEED, a seed for the random number
        //    generator.
        //
        //    Output, double R83T_RANDOM[M*3], the matrix.
        //
    {
        int i;

        double[] a = new double[m * 3];

        for (i = 0; i < m; i++)
        {
            int k;
            for (k = 0; k < 3; k++)
            {
                int j = i + k - 1;
                a[i + k * m] = j switch
                {
                    >= 0 when j <= n - 1 => UniformRNG.r8_uniform_01(ref seed),
                    _ => 0.0
                };
            }
        }

        return a;
    }

    public static double[] r83t_res(int m, int n, double[] a, double[] x, double[] b)

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

        double[] r = r83t_mv(m, n, a, x);

        for (i = 0; i < m; i++)
        {
            r[i] = b[i] - r[i];
        }

        return r;
    }

    public static double[] r83t_to_r8ge(int m, int n, double[] a_r83t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83T_TO_R8GE copies an R83T matrix to an R8GE matrix.
        //
        //  Discussion:
        //
        //    The R83T storage format is used for an MxN tridiagonal matrix.
        //    The superdiagonal is stored in entries (1:M-1,3), the diagonal in
        //    entries (1:M,2), and the subdiagonal in (2:M,1).  Thus, the
        //    the rows of the original matrix slide horizontally to form an
        //    Mx3 stack of data.
        //
        //    An R83T matrix of order 3x5 would be stored:
        //
        //       *  A11 A12
        //      A21 A22 A23
        //      A32 A33 A34
        //
        //    An R83T matrix of order 5x5 would be stored:
        //
        //       *  A11 A12
        //      A21 A22 A23
        //      A32 A33 A34
        //      A43 A44 A45
        //      A54 A55  *
        //
        //    An R83T matrix of order 5x3 would be stored:
        //
        //       *  A11 A12
        //      A21 A22 A23
        //      A32 A33  *
        //      A43  *   *
        //       *   *   *
        //
        //    The R8GE storage format is used for a general M by N matrix.  A storage 
        //    space is made for each entry.  The two dimensional logical
        //    array can be thought of as a vector of M*N entries, starting with
        //    the M entries in the column 1, then the M entries in column 2
        //    and so on.  Considered as a vector, the entry A(I,J) is then stored
        //    in vector location I+(J-1)*M.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 May 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Input, double A_R83T[M*3], the R83T matrix.
        //
        //    Output, double R83T_TO_R8GE[M*N], the R8GE matrix.
        //
    {
        int i;
        int j;

        double[] a_r8ge = new double[m * n];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                a_r8ge[i + j * m] = 0.0;
            }
        }

        for (i = 0; i < m; i++)
        {
            int k;
            for (k = 0; k < 3; k++)
            {
                j = i + k - 1;
                a_r8ge[i + j * m] = j switch
                {
                    >= 0 when j <= n - 1 => a_r83t[i + k * m],
                    _ => a_r8ge[i + j * m]
                };
            }
        }

        return a_r8ge;
    }

    public static double[] r83t_zeros(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83T_ZEROS zeros an R83T matrix.
        //
        //  Discussion:
        //
        //    The R83T storage format is used for an MxN tridiagonal matrix.
        //    The superdiagonal is stored in entries (1:M-1,3), the diagonal in
        //    entries (1:M,2), and the subdiagonal in (2:M,1).  Thus, the
        //    the rows of the original matrix slide horizontally to form an
        //    Mx3 stack of data.
        //
        //    An R83T matrix of order 3x5 would be stored:
        //
        //       *  A11 A12
        //      A21 A22 A23
        //      A32 A33 A34
        //
        //    An R83T matrix of order 5x5 would be stored:
        //
        //       *  A11 A12
        //      A21 A22 A23
        //      A32 A33 A34
        //      A43 A44 A45
        //      A54 A55  *
        //
        //    An R83T matrix of order 5x3 would be stored:
        //
        //       *  A11 A12
        //      A21 A22 A23
        //      A32 A33  *
        //      A43  *   *
        //       *   *   *
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 May 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Output, double R83T_ZEROS[M*3], the matrix.
        //
    {
        int j;

        double[] a = new double[m * 3];
        for (j = 0; j < 3; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                a[i + j * m] = 0.0;
            }
        }

        return a;
    }
}