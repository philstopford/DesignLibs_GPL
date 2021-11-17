using System;
using Burkardt.Uniform;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r83s_cg(int n, double[] a, double[] b, ref double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83S_CG uses the conjugate gradient method on an R83S system.
        //
        //  Discussion:
        //
        //    The R83S storage format is used for a tridiagonal scalar matrix.
        //    The vector A(3) contains the subdiagonal, diagonal, and superdiagonal
        //    values that occur on every row.
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
        //    Here is how an R83S matrix of order 5, stored as (A1,A2,A3), would
        //    be interpreted:
        //
        //      A2  A3   0   0   0
        //      A1  A2  A3   0   0
        //       0  A1  A2  A3   0 
        //       0   0  A1  A2  A3
        //       0   0   0  A1  A2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 July 2014
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
        //    Input, double A[3], the matrix.
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
        ap = r83s_mv(n, n, a, x);

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
            ap = r83s_mv(n, n, a, p);
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

    public static double[] r83s_dif2(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83S_DIF2 returns the DIF2 matrix in R83S format.
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
        //    A is tridiagonal.
        //    Because A is tridiagonal, it has property A (bipartite).
        //    A is a special case of the TRIS or tridiagonal scalar matrix.
        //    A is integral, therefore det ( A ) is integral, and 
        //    det ( A ) * inverse ( A ) is integral.
        //    A is Toeplitz: constant along diagonals.
        //    A is symmetric: A' = A.
        //    Because A is symmetric, it is normal.
        //    Because A is normal, it is diagonalizable.
        //    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).
        //    A is positive definite.
        //    A is an M matrix.
        //    A is weakly diagonally dominant, but not strictly diagonally dominant.
        //    A has an LU factorization A = L * U, without pivoting.
        //      The matrix L is lower bidiagonal with subdiagonal elements:
        //        L(I+1,I) = -I/(I+1)
        //      The matrix U is upper bidiagonal, with diagonal elements
        //        U(I,I) = (I+1)/I
        //      and superdiagonal elements which are all -1.
        //    A has a Cholesky factorization A = L * L', with L lower bidiagonal.
        //      L(I,I) =    sqrt ( (I+1) / I )
        //      L(I,I-1) = -sqrt ( (I-1) / I )
        //    The eigenvalues are
        //      LAMBDA(I) = 2 + 2 * COS(I*PI/(N+1))
        //                = 4 SIN^2(I*PI/(2*N+2))
        //    The corresponding eigenvector X(I) has entries
        //       X(I)(J) = sqrt(2/(N+1)) * sin ( I*J*PI/(N+1) ).
        //    Simple linear systems:
        //      x = (1,1,1,...,1,1),   A*x=(1,0,0,...,0,1)
        //      x = (1,2,3,...,n-1,n), A*x=(0,0,0,...,0,n+1)
        //    det ( A ) = N + 1.
        //    The value of the determinant can be seen by induction,
        //    and expanding the determinant across the first row:
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
        //    09 July 2014
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
        //    Output, double A[3], the matrix.
        //
    {
        double[] a = new double[3];

        a[0] = -1.0;
        a[1] = 2.0;
        a[2] = -1.0;

        return a;
    }

    public static void r83s_gs_sl(int n, double[] a, double[] b, ref double[] x, int it_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83S_GS_SL solves an R83S system using Gauss-Seidel iteration.
        //
        //  Discussion:
        //
        //    The R83S storage format is used for a tridiagonal scalar matrix.
        //    The vector A(3) contains the subdiagonal, diagonal, and superdiagonal
        //    values that occur on every row.
        //    RGE A(I,J) = R83S A[I-J+1].
        //
        //    This routine simply applies a given number of steps of the
        //    iteration to an input approximate solution.  On first call, you can
        //    simply pass in the zero vector as an approximate solution.  If
        //    the returned value is not acceptable, you may call again, using
        //    it as the starting point for additional iterations.
        //
        //  Example:
        //
        //    Here is how an R83S matrix of order 5, stored as (A1,A2,A3), would
        //    be interpreted:
        //
        //      A2  A1   0   0   0
        //      A3  A2  A1   0   0
        //       0  A3  A2  A1   0 
        //       0   0  A3  A2  A1
        //       0   0   0  A3  A2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 September 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[3], the R83S matrix.
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
        switch (a[1])
        {
            //
            //  No diagonal matrix entry can be zero.
            //
            case 0.0:
                Console.WriteLine("");
                Console.WriteLine("R83S_GS_SL - Fatal error!");
                Console.WriteLine("  Zero diagonal entry");
                return;
        }

        for (it_num = 1; it_num <= it_max; it_num++)
        {
            x[0] = (b[0] - a[2] * x[1]) / a[1];
            for (i = 1; i < n - 1; i++)
            {
                x[i] = (b[i] - a[0] * x[i - 1] - a[2] * x[i + 1]) / a[1];
            }

            x[n - 1] = (b[n - 1] - a[0] * x[n - 2]) / a[1];
        }

    }

    public static double[] r83s_indicator(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83S_INDICATOR sets up an R83s indicator matrix.
        //
        //  Discussion:
        //
        //    The R83S storage format is used for a tridiagonal scalar matrix.
        //    The vector A(3) contains the subdiagonal, diagonal, and superdiagonal
        //    values that occur on every row.
        //    RGE A(I,J) = R83S A[I-J+1].
        //
        //  Example:
        //
        //    Here is how an R83S matrix of order 5, stored as (A1,A2,A3), would
        //    be interpreted:
        //
        //      A2  A1   0   0   0
        //      A3  A2  A1   0   0
        //       0  A3  A2  A1   0 
        //       0   0  A3  A2  A1
        //       0   0   0  A3  A2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 September 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Output, double R83_INDICATOR[3], the R83S indicator matrix.
        //
    {
        double[] a;

        a = new double[3];

        a[0] = 3.0;
        a[1] = 2.0;
        a[2] = 1.0;

        return a;
    }

    public static void r83s_jac_sl(int n, double[] a, double[] b, ref double[] x, int it_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83S_JAC_SL solves an R83S system using Jacobi iteration.
        //
        //  Discussion:
        //
        //    The R83S storage format is used for a tridiagonal scalar matrix.
        //    The vector A(3) contains the subdiagonal, diagonal, and superdiagonal
        //    values that occur on every row.
        //    RGE A(I,J) = R83S A[I-J+1].
        //
        //    This routine simply applies a given number of steps of the
        //    iteration to an input approximate solution.  On first call, you can
        //    simply pass in the zero vector as an approximate solution.  If
        //    the returned value is not acceptable, you may call again, using
        //    it as the starting point for additional iterations.
        //
        //  Example:
        //
        //    Here is how an R83S matrix of order 5, stored as (A1,A2,A3), would
        //    be interpreted:
        //
        //      A2  A1   0   0   0
        //      A3  A2  A1   0   0
        //       0  A3  A2  A1   0 
        //       0   0  A3  A2  A1
        //       0   0   0  A3  A2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 September 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[3], the R83S matrix.
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
        double[] xnew;

        xnew = new double[n];
        switch (a[1])
        {
            //
            //  No diagonal matrix entry can be zero.
            //
            case 0.0:
                Console.WriteLine("");
                Console.WriteLine("R83S_JAC_SL - Fatal error!");
                Console.WriteLine("  Zero diagonal entry");
                return;
        }

        for (it_num = 1; it_num <= it_max; it_num++)
        {
            //
            //  Solve A*x=b:
            //
            xnew[0] = b[0] - a[2] * x[1];
            for (i = 1; i < n - 1; i++)
            {
                xnew[i] = b[i] - a[0] * x[i - 1] - a[2] * x[i + 1];
            }

            xnew[n - 1] = b[n - 1] - a[0] * x[n - 2];
            //
            //  Divide by the diagonal term, and overwrite X.
            //
            for (i = 0; i < n; i++)
            {
                xnew[i] /= a[1];
            }

            for (i = 0; i < n; i++)
            {
                x[i] = xnew[i];
            }
        }

    }

    public static double[] r83s_mtv(int m, int n, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83S_MTV multiplies a vector times an R83S matrix.
        //
        //  Discussion:
        //
        //    The R83S storage format is used for a tridiagonal scalar matrix.
        //    The vector A(3) contains the subdiagonal, diagonal, and superdiagonal
        //    values that occur on every row.
        //    RGE A(I,J) = R83S A[I-J+1].
        //
        //  Example:
        //
        //    Here is how an R83S matrix of order 5, stored as (A1,A2,A3), would
        //    be interpreted:
        //
        //      A2  A1   0   0   0
        //      A3  A2  A1   0   0
        //       0  A3  A2  A1   0 
        //       0   0  A3  A2  A1
        //       0   0   0  A3  A2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 September 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the linear system.
        //
        //    Input, double A[3], the R83S matrix.
        //
        //    Input, double X[M], the vector to be multiplied.
        //
        //    Output, double R83S_MTV[N], the product A'*x.
        //
    {
        double[] b;
        int i;
        int i_hi;
        int i_lo;
        int j;

        b = new double[n];

        for (i = 0; i < n; i++)
        {
            b[i] = 0.0;
        }

        for (j = 0; j < n; j++)
        {
            i_lo = Math.Max(0, j - 1);
            i_hi = Math.Min(m - 1, j + 1);
            for (i = i_lo; i <= i_hi; i++)
            {
                b[j] += x[i] * a[i - j + 1];
            }
        }

        return b;
    }


    public static double[] r83s_mv(int m, int n, double[] a, double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83S_MV multiplies a R83S matrix times a vector.
        //
        //  Discussion:
        //
        //    The R83S storage format is used for a tridiagonal scalar matrix.
        //    The vector A(3) contains the subdiagonal, diagonal, and superdiagonal
        //    values that occur on every row.
        //
        //  Example:
        //
        //    Here is how an R83S matrix of order 5, stored as (A1,A2,A3), would
        //    be interpreted:
        //
        //      A2  A3   0   0   0
        //      A1  A2  A3   0   0
        //       0  A1  A2  A3   0 
        //       0   0  A1  A2  A3
        //       0   0   0  A1  A2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A[3], the matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double R83S_MV[M], the product A * x.
        //
    {
        int i;
        int ihi;

        double[] b = new double[m];

        for (i = 0; i < m; i++)
        {
            b[i] = 0.0;
        }

        ihi = Math.Min(m, n + 1);

        for (i = 1; i < ihi; i++)
        {
            b[i] += a[0] * x[i - 1];
        }

        ihi = Math.Min(m, n);
        for (i = 0; i < ihi; i++)
        {
            b[i] += a[1] * x[i];
        }

        ihi = Math.Min(m, n - 1);
        for (i = 0; i < ihi; i++)
        {
            b[i] += a[2] * x[i + 1];
        }

        return b;
    }


    public static void r83s_print(int m, int n, double[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83S_PRINT prints an R83S matrix.
        //
        //  Discussion:
        //
        //    The R83S storage format is used for a tridiagonal scalar matrix.
        //    The vector A(3) contains the subdiagonal, diagonal, and superdiagonal
        //    values that occur on every row.
        //    RGE A(I,J) = R83S A[I-J+1].
        //
        //  Example:
        //
        //    Here is how an R83S matrix of order 5, stored as (A1,A2,A3), would
        //    be interpreted:
        //
        //      A2  A1   0   0   0
        //      A3  A2  A1   0   0
        //       0  A3  A2  A1   0 
        //       0   0  A3  A2  A1
        //       0   0   0  A3  A2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 September 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[3], the R83S matrix.
        //
        //    Input, string TITLE, a title.
        //
    {
        r83s_print_some(m, n, a, 1, 1, m, n, title);
    }

    public static void r83s_print_some(int m, int n, double[] a, int ilo, int jlo, int ihi,
            int jhi, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83S_PRINT_SOME prints some of an R83S matrix.
        //
        //  Discussion:
        //
        //    The R83S storage format is used for a tridiagonal scalar matrix.
        //    The vector A(3) contains the subdiagonal, diagonal, and superdiagonal
        //    values that occur on every row.
        //    RGE A(I,J) = R83S A[I-J+1].
        //
        //  Example:
        //
        //    Here is how an R83S matrix of order 5, stored as (A1,A2,A3), would
        //    be interpreted:
        //
        //      A2  A1   0   0   0
        //      A3  A2  A1   0   0
        //       0  A3  A2  A1   0 
        //       0   0  A3  A2  A1
        //       0   0   0  A3  A2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 September 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Input, double A[3], the R83S matrix.
        //
        //    Input, int ILO, JLO, IHI, JHI, designate the first row and
        //    column, and the last row and column, to be printed.
        //
        //    Input, string TITLE, a title.
        //
    {
        int INCX = 5;

        int i;
        int i2hi;
        int i2lo;
        int inc;
        int j;
        int j2;
        int j2hi;
        int j2lo;
        string cout = "";

        Console.WriteLine("");
        Console.WriteLine(title + "");
        //
        //  Print the columns of the matrix, in strips of 5.
        //
        for (j2lo = jlo; j2lo <= jhi; j2lo += INCX)
        {
            j2hi = j2lo + INCX - 1;
            j2hi = Math.Min(j2hi, n);
            j2hi = Math.Min(j2hi, jhi);

            inc = j2hi + 1 - j2lo;

            Console.WriteLine("");
            cout = "  Col: ";
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
            i2lo = Math.Max(ilo, 1);
            i2lo = Math.Max(i2lo, j2lo - 1);

            i2hi = Math.Min(ihi, m);
            i2hi = Math.Min(i2hi, j2hi + 1);

            for (i = i2lo; i <= i2hi; i++)
            {
                //
                //  Print out (up to) 5 entries in row I, that lie in the current strip.
                //
                cout = i.ToString().PadLeft(6) + ":   ";

                for (j2 = 1; j2 <= inc; j2++)
                {
                    j = j2lo - 1 + j2;

                    switch (i - j + 1)
                    {
                        case < 0:
                        case > 2:
                            cout += "              ";
                            break;
                        default:
                            cout += "  " + a[i - j + 1].ToString().PadLeft(12);
                            break;
                    }
                }

                Console.WriteLine(cout);
            }
        }
    }

    public static double[] r83s_random(int m, int n, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83S_RANDOM randomizes an R83S matrix.
        //
        //  Discussion:
        //
        //    The R83S storage format is used for a tridiagonal scalar matrix.
        //    The vector A(3) contains the subdiagonal, diagonal, and superdiagonal
        //    values that occur on every row.
        //    RGE A(I,J) = R83S A[I-J+1].
        //
        //  Example:
        //
        //    Here is how an R83S matrix of order 5, stored as (A1,A2,A3), would
        //    be interpreted:
        //
        //      A2  A1   0   0   0
        //      A3  A2  A1   0   0
        //       0  A3  A2  A1   0 
        //       0   0  A3  A2  A1
        //       0   0   0  A3  A2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 September 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double R83S_RANDOM[3], the R83S matrix.
        //
    {
        double[] a;

        a = UniformRNG.r8vec_uniform_01_new(3, ref seed);

        return a;
    }


    public static double[] r83s_res(int m, int n, double[] a, double[] x, double[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83S_RES computes the residual R = B-A*X for R83S matrices.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 July 2014
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
        //    Input, double A[3], the matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Input, double B[M], the desired result A * x.
        //
        //    Output, double R83S_RES[M], the residual R = B - A * X.
        //
    {
        int i;
        double[] r;

        r = r83s_mv(m, n, a, x);

        for (i = 0; i < m; i++)
        {
            r[i] = b[i] - r[i];
        }

        return r;
    }

    public static double[] r83s_to_r8ge(int m, int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83S_TO_R8GE copies an R83S matrix to an R8GE matrix.
        //
        //  Discussion:
        //
        //    The R83S storage format is used for a tridiagonal scalar matrix.
        //    The vector A(3) contains the subdiagonal, diagonal, and superdiagonal
        //    values that occur on every row.
        //    RGE A(I,J) = R83S A[I-J+1].
        //
        //  Example:
        //
        //    Here is how an R83S matrix of order 5, stored as (A1,A2,A3), would
        //    be interpreted:
        //
        //      A2  A1   0   0   0
        //      A3  A2  A1   0   0
        //       0  A3  A2  A1   0 
        //       0   0  A3  A2  A1
        //       0   0   0  A3  A2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 September 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Input, double A[3], the R83S matrix.
        //
        //    Output, double R83S_TO_R8GE[M*N], the R8GE matrix.
        //
    {
        double[] b;
        int i;
        int i_hi;
        int i_lo;
        int j;

        b = new double[m * n];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                b[i + j * m] = 0.0;
            }
        }

        for (j = 0; j < n; j++)
        {
            i_lo = Math.Max(0, j - 1);
            i_hi = Math.Min(m - 1, j + 1);
            for (i = i_lo; i <= i_hi; i++)
            {
                b[i + j * m] = a[i - j + 1];
            }
        }

        return b;
    }

    public static double[] r83s_zeros(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83S_ZEROS zeros an R83S matrix.
        //
        //  Discussion:
        //
        //    The R83S storage format is used for a tridiagonal scalar matrix.
        //    The vector A(3) contains the subdiagonal, diagonal, and superdiagonal
        //    values that occur on every row.
        //    RGE A(I,J) = R83S A[I-J+1].
        //
        //  Example:
        //
        //    Here is how an R83S matrix of order 5, stored as (A1,A2,A3), would
        //    be interpreted:
        //
        //      A2  A1   0   0   0
        //      A3  A2  A1   0   0
        //       0  A3  A2  A1   0 
        //       0   0  A3  A2  A1
        //       0   0   0  A3  A2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 September 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the linear system.
        //
        //    Output, double R83S_ZEROS[3], the R83S matrix.
        //
    {
        double[] a;
        int i;

        a = new double[3];

        for (i = 0; i < 3; i++)
        {
            a[i] = 0.0;
        }

        return a;
    }

}