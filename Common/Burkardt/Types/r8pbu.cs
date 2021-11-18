using System;
using System.Globalization;
using Burkardt.Uniform;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r8pbu_cg(int n, int mu, double[] a, double[] b, ref double[] x)

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
        int i;
        int it;
        //
        //  Initialize
        //    AP = A * x,
        //    R  = b - A * x,
        //    P  = b - A * x.
        //
        double[] ap = r8pbu_mv(n, n, mu, a, x);

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
            ap = r8pbu_mv(n, n, mu, a, p);
            //
            //  Compute the dot products
            //    PAP = P*AP,
            //    PR  = P*R
            //  Set
            //    ALPHA = PR / PAP.
            //
            double pap = 0.0;
            for (i = 0; i < n; i++)
            {
                pap += p[i] * ap[i];
            }

            if (pap == 0.0)
            {
                break;
            }

            double pr = 0.0;
            for (i = 0; i < n; i++)
            {
                pr += p[i] * r[i];
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
            double rap = 0.0;
            for (i = 0; i < n; i++)
            {
                rap += r[i] * ap[i];
            }

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

    public static double r8pbu_det(int n, int mu, double[] a_lu)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PBU_DET computes the determinant of a matrix factored by R8PBU_FA.
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
        //    03 October 2003
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
        //    LINPACK User's Guide,
        //    SIAM, 1979,
        //    ISBN13: 978-0-898711-72-1,
        //    LC: QA214.L56.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be positive.
        //
        //    Input, int MU, the number of superdiagonals of the matrix.
        //    MU must be at least 0 and no more than N-1.
        //
        //    Input, double A_LU[(MU+1)*N], the LU factors from R8PBU_FA.
        //
        //    Output, double R8PBU_DET, the determinant of the matrix.
        //
    {
        int j;

        double det = 1.0;

        for (j = 0; j < n; j++)
        {
            det = det * a_lu[mu + j * (mu + 1)] * a_lu[mu + j * (mu + 1)];
        }

        return det;
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
        //      L(I,I) =    Math.Sqrt ( (I+1) / I )
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

    public static double[] r8pbu_fa(int n, int mu, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PBU_FA factors an R8PBU matrix.
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
        //    Once factored, linear systems A*x=b involving the matrix can be solved
        //    by calling R8PBU_SL.  No pivoting is performed.  Pivoting is not necessary
        //    for positive definite symmetric matrices.  If the matrix is not positive
        //    definite, the algorithm may behave correctly, but it is also possible
        //    that an illegal divide by zero will occur.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 February 2004
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
        //    LINPACK User's Guide,
        //    SIAM, 1979,
        //    ISBN13: 978-0-898711-72-1,
        //    LC: QA214.L56.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be positive.
        //
        //    Input, int MU, the number of superdiagonals of the matrix.
        //    MU must be at least 0, and no more than N-1.
        //
        //    Input, double A[(MU+1)*N], the R8PBU matrix.
        //
        //    Output, double R8PBU_FA[(MU+1)*N], information describing a factored
        //    form of the matrix.
        //
    {
        int i;
        int j;

        double[] b = new double[(mu + 1) * n];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < mu + 1; i++)
            {
                b[i + j * (mu + 1)] = a[i + j * (mu + 1)];
            }
        }

        for (j = 1; j <= n; j++)
        {

            int ik = mu + 1;
            int jk = Math.Max(j - mu, 1);
            int mm = Math.Max(mu + 2 - j, 1);

            double s = 0.0;

            int k;
            for (k = mm; k <= mu; k++)
            {
                double t = 0.0;
                for (i = 0; i <= k - mm - 1; i++)
                {
                    t += b[ik + i - 1 + (jk - 1) * (mu + 1)] * b[mm + i - 1 + (j - 1) * (mu + 1)];
                }

                b[k - 1 + (j - 1) * (mu + 1)] = (b[k - 1 + (j - 1) * (mu + 1)] - t) /
                                                b[mu + (jk - 1) * (mu + 1)];

                s += b[k - 1 + (j - 1) * (mu + 1)] * b[k - 1 + (j - 1) * (mu + 1)];
                ik -= 1;
                jk += 1;
            }

            s = b[mu + (j - 1) * (mu + 1)] - s;

            switch (s)
            {
                case <= 0.0:
                    return null;
                default:
                    b[mu + (j - 1) * (mu + 1)] = Math.Sqrt(s);
                    break;
            }
        }

        return b;
    }

    public static double[] r8pbu_indicator(int n, int mu)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PBU_INDICATOR sets up an R8PBU indicator matrix.
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
        //    15 January 2004
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
        //    Input, int MU, the number of superdiagonals in the matrix.
        //    MU must be at least 0 and no more than N-1.
        //
        //    Output, double R8PBU_INDICATOR[(MU+1)*N], the R8PBU matrix.
        //
    {
        int i;
        int j;

        double[] a = new double[(mu + 1) * n];

        int fac = (int) Math.Pow(10, (int) Math.Log10(n) + 1);
        //
        //  Zero out the "junk" entries.
        //
        for (j = 0; j < mu; j++)
        {
            for (i = 0; i <= mu - j; i++)
            {
                a[i + j * (mu + 1)] = 0.0;
            }
        }

        //
        //  Set the meaningful values.
        //
        for (i = 1; i <= n; i++)
        {
            for (j = i; j <= Math.Min(i + mu, n); j++)
            {
                a[mu + i - j + (j - 1) * (mu + 1)] = fac * i + j;
            }
        }

        return a;
    }

    public static double[] r8pbu_ml(int n, int mu, double[] a_lu, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PBU_ML multiplies a vector times a matrix that was factored by R8PBU_FA.
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
        //    06 October 2003
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
        //    Input, int MU, the number of superdiagonals of the matrix.
        //    MU must be at least 0 and no more than N-1.
        //
        //    Input, double A_LU[(MU+1)*N], the LU factors from R8PBU_FA.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double R8PBU_ML[N], the product A * x.
        //
    {
        int i;
        int k;

        double[] b = new double[n];

        for (i = 0; i < n; i++)
        {
            b[i] = x[i];
        }

        //
        //  Multiply U * X = Y.
        //
        for (k = 1; k <= n; k++)
        {
            int ilo = Math.Max(1, k - mu);
            for (i = ilo; i <= k - 1; i++)
            {
                b[i - 1] += a_lu[mu + i - k + (k - 1) * (mu + 1)] * b[k - 1];
            }

            b[k - 1] = a_lu[mu + (k - 1) * (mu + 1)] * b[k - 1];
        }

        //
        //  Multiply L * Y = B.
        //
        for (k = n; 1 <= k; k--)
        {
            int jhi = Math.Min(k + mu, n);
            int j;
            for (j = k + 1; j <= jhi; j++)
            {
                b[j - 1] += a_lu[mu + k - j + (j - 1) * (mu + 1)] * b[k - 1];
            }

            b[k - 1] = a_lu[mu + (k - 1) * (mu + 1)] * b[k - 1];
        }

        return b;
    }

    public static double[] r8pbu_mv(int m, int n, int mu, double[] a, double[] x)

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
                b[ieqn - 1] += a[i - 1 + (j - 1) * (mu + 1)] * x[j - 1];
                b[j - 1] += a[i - 1 + (j - 1) * (mu + 1)] * x[ieqn - 1];
            }
        }

        return b;
    }

    public static void r8pbu_print(int n, int mu, double[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PBU_PRINT prints an R8PBU matrix.
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
        //    Input, int MU, the upper (and lower) bandwidth.
        //    MU must be nonnegative, and no greater than N-1.
        //
        //    Input, double A[(MU+1)*N], the R8PBU matrix.
        //
        //    Input, string TITLE, a title.
        //
    {
        r8pbu_print_some(n, mu, a, 1, 1, n, n, title);
    }

    public static void r8pbu_print_some(int n, int mu, double[] a, int ilo, int jlo, int ihi,
            int jhi, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PBU_PRINT_SOME prints some of an R8PBU matrix.
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
        //    Input, int MU, the upper (and lower) bandwidth.
        //    MU must be nonnegative, and no greater than N-1.
        //
        //    Input, double A[(MU+1)*N], the R8PBU matrix.
        //
        //    Input, int ILO, JLO, IHI, JHI, designate the first row and
        //    column, and the last row and column to be printed.
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

            Console.WriteLine("");

            string cout = "  Col: ";
            int j;
            for (j = j2lo; j <= j2hi; j++)
            {
                cout += j.ToString(CultureInfo.InvariantCulture).PadLeft(7) + "       ";
            }

            Console.WriteLine(cout);
            Console.WriteLine("  Row");
            Console.WriteLine("  ---");
            //
            //  Determine the range of the rows in this strip.
            //
            int i2lo = Math.Max(ilo, 1);
            i2lo = Math.Max(i2lo, j2lo - mu);
            int i2hi = Math.Min(ihi, n);
            i2hi = Math.Min(i2hi, j2hi + mu);

            int i;
            for (i = i2lo; i <= i2hi; i++)
            {
                cout = i.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  ";
                //
                //  Print out (up to) 5 entries in row I, that lie in the current strip.
                //
                for (j = j2lo; j <= j2hi; j++)
                {
                    if (mu < i - j || mu < j - i)
                    {
                        cout += "              ";
                    }
                    else if (i <= j && j <= i + mu)
                    {
                        cout += a[mu + i - j + (j - 1) * (mu + 1)].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  ";
                    }
                    else if (i - mu <= j && j <= i)
                    {
                        cout += a[mu + j - i + (i - 1) * (mu + 1)].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  ";
                    }
                }

                Console.WriteLine(cout);
            }
        }
    }

    public static double[] r8pbu_random(int n, int mu, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PBU_RANDOM randomizes an R8PBU matrix.
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
        //    The matrix returned will be positive definite, but of limited
        //    randomness.  The off diagonal elements are random values between
        //    0 and 1, and the diagonal element of each row is selected to
        //    ensure strict diagonal dominance.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 January 2004
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
        //    Input, int MU, the number of superdiagonals in the matrix.
        //    MU must be at least 0 and no more than N-1.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double R8PBU_RANDOM[(MU+1)*N], the R8PBU matrix.
        //
    {
        int i;
        int j;

        double[] a = new double[(mu + 1) * n];
        //
        //  Zero out the "junk" entries.
        //
        for (j = 0; j < mu; j++)
        {
            for (i = 0; i <= mu - j; i++)
            {
                a[i + j * (mu + 1)] = 0.0;
            }
        }

        //
        //  Set the off diagonal values.
        //
        for (i = 0; i < n; i++)
        {
            for (j = i + 1; j <= Math.Min(i + mu, n - 1); j++)
            {
                a[mu + i - j + j * (mu + 1)] = UniformRNG.r8_uniform_01(ref seed);
            }
        }

        //
        //  Set the diagonal values.
        //
        for (i = 1; i <= n; i++)
        {
            double sum2 = 0.0;

            int jlo = Math.Max(1, i - mu);
            for (j = jlo; j <= i - 1; j++)
            {
                sum2 += Math.Abs(a[mu + j - i + (i - 1) * (mu + 1)]);
            }

            int jhi = Math.Min(i + mu, n);
            for (j = i + 1; j <= jhi; j++)
            {
                sum2 += Math.Abs(a[mu + i - j + (j - 1) * (mu + 1)]);
            }

            double r = UniformRNG.r8_uniform_01(ref seed);

            a[mu + (i - 1) * (mu + 1)] = (1.0 + r) * (sum2 + 0.01);

        }

        return a;
    }

    public static double[] r8pbu_res(int m, int n, int mu, double[] a, double[] x, double[] b)

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

    public static double[] r8pbu_sl(int n, int mu, double[] a_lu, double[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PBU_SL solves an R8PBU system factored by R8PBU_FA.
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
        //    04 June 2016
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
        //    LINPACK User's Guide,
        //    SIAM, 1979,
        //    ISBN13: 978-0-898711-72-1,
        //    LC: QA214.L56.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be positive.
        //
        //    Input, int MU, the number of superdiagonals of the matrix.
        //    MU must be at least 0 and no more than N-1.
        //
        //    Input, double A_LU[(MU+1)*N], the LU factors from R8PBU_FA.
        //
        //    Input, double B[N], the right hand side of the linear system.
        //
        //    Output, double R8PBU_SL[N], the solution vector.
        //
    {
        int i;
        int ilo;
        int k;

        double[] x = new double[n];

        for (k = 0; k < n; k++)
        {
            x[k] = b[k];
        }

        //
        //  Solve L * Y = B.
        //
        for (k = 1; k <= n; k++)
        {
            ilo = Math.Max(1, k - mu);
            double t = 0.0;
            for (i = ilo; i <= k - 1; i++)
            {
                t += x[i - 1] * a_lu[mu + i - k + (k - 1) * (mu + 1)];
            }

            x[k - 1] = (x[k - 1] - t) / a_lu[mu + (k - 1) * (mu + 1)];
        }

        //
        //  Solve U * X = Y.
        //
        for (k = n; 1 <= k; k--)
        {
            x[k - 1] /= a_lu[mu + (k - 1) * (mu + 1)];

            ilo = Math.Max(1, k - mu);
            for (i = ilo; i <= k - 1; i++)
            {
                x[i - 1] -= x[k - 1] * a_lu[mu + i - k + (k - 1) * (mu + 1)];
            }
        }

        return x;
    }

    public static double[] r8pbu_sor(int n, int mu, double[] a, double[] b, double eps, int itchk,
            int itmax, double omega, double[] x_init )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PBU_SOR uses SOR iteration to solve an R8PBU linear system.
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
        //    A relaxation factor OMEGA may be used.
        //
        //    The iteration will proceed until a convergence test is met,
        //    or the iteration limit is reached.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 October 2003
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
        //    Input, int MU, the number of superdiagonals in the matrix.
        //    MU must be at least 0, and no more than N-1.
        //
        //    Input, double A[(MU+1)*N], the R8PBU matrix.
        //
        //    Input, double B[N], the right hand side of the system.
        //
        //    Input, double EPS, convergence tolerance for the system.  The vector
        //    b - A * x is computed every ITCHK iterations, and if the maximum
        //    entry of this vector is of norm less than EPS, the program
        //    will return.
        //
        //    Input, int ITCHK, the interval between convergence checks.  ITCHK steps
        //    will be taken before any check is made on whether the iteration
        //    has converged.  ITCHK should be at least 1 and no greater
        //    than ITMAX.
        //
        //    Input, int ITMAX, the maximum number of iterations allowed.  The
        //    program will return to the user if this many iterations are taken
        //    without convergence.
        //
        //    Input, double OMEGA, the relaxation factor.  OMEGA must be strictly between
        //    0 and 2.  Use OMEGA = 1 for no relaxation, classical Jacobi iteration.
        //
        //    Input, double X_INIT[N], a starting vector for the iteration.
        //
        //    Output, double R8PBU_SOR[N], the approximation to the solution.
        //
    {
        int i;

        if (itchk <= 0 || itmax < itchk)
        {
            Console.WriteLine("");
            Console.WriteLine("R8PBU_SOR - Fatal error!");
            Console.WriteLine("  Illegal ITCHK = " + itchk + "");
            return null;
        }

        switch (omega)
        {
            case <= 0.0:
            case >= 2.0:
                Console.WriteLine("");
                Console.WriteLine("R8PBU_SOR - Fatal error!");
                Console.WriteLine("  Illegal value of OMEGA = " + omega + "");
                return null;
        }

        double[] x = new double[n];
        for (i = 0; i < n; i++)
        {
            x[i] = x_init[i];
        }

        //
        //  Take ITCHK steps of the iteration before doing a convergence check.
        //
        while (true)
        {
            double[] xtemp;
            int it;
            for (it = 1; it <= itchk; it++)
            {
                //
                //  Compute XTEMP(I) = B(I) + A(I,I) * X(I) - SUM ( J=1 to N ) A(I,J) * X(J).
                //
                xtemp = r8pbu_mv(n, n, mu, a, x);

                for (i = 0; i < n; i++)
                {
                    xtemp[i] = x[i] + (b[i] - xtemp[i]) / a[mu + i * (mu + 1)];
                }

                //
                //  Compute the next iterate as a weighted combination of the
                //  old iterate and the just computed standard Jacobi iterate.
                //
                if (Math.Abs(omega - 1.0) > double.Epsilon)
                {
                    for (i = 0; i < n; i++)
                    {
                        xtemp[i] = (1.0 - omega) * x[i] + omega * xtemp[i];
                    }
                }

                //
                //  Copy the new result into the old result vector.
                //
                for (i = 0; i < n; i++)
                {
                    x[i] = xtemp[i];
                }
            }

            //
            //  Compute the maximum residual, the greatest entry in the vector
            //  RESID(I) = B(I) - A(I,J) * X(J).
            //
            xtemp = r8pbu_mv(n, n, mu, a, x);

            double err = 0.0;
            for (i = 0; i < n; i++)
            {
                err = Math.Max(err, Math.Abs(b[i] - xtemp[i]));
            }

            //
            //  Test to see if we can quit because of convergence,
            //
            if (err <= eps)
            {
                return x;
            }

        }
    }

    public static double[] r8pbu_to_r8ge(int n, int mu, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PBU_TO_R8GE copies an R8PBU matrix to an R8GE matrix.
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
        //    19 May 2016
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
        //    Input, int MU, the upper bandwidth of A1.
        //    MU must be nonnegative, and no greater than N-1.
        //
        //    Input, double A[(MU+1)*N], the R8PBU matrix.
        //
        //    Output, double R8PBU_TO_R8GE[N*N], the R8GE matrix.
        //
    {
        int i;
        int j;

        double[] b = new double[n * n];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                b[i + j * n] = 0.0;
            }
        }

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                if (i <= j && j <= i + mu)
                {
                    b[i + j * n] = a[mu + i - j + j * (mu + 1)];
                }
                else if (i - mu <= j && j < i)
                {
                    b[i + j * n] = a[mu + j - i + i * (mu + 1)];
                }
                else
                {
                    b[i + j * n] = 0.0;
                }
            }
        }

        return b;
    }

    public static double[] r8pbu_zeros(int n, int mu)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PBU_ZEROS zeros an R8PBU matrix.
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
        //    15 September 2003
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
        //    Input, int MU, the number of superdiagonals in the matrix.
        //    MU must be at least 0 and no more than N-1.
        //
        //    Output, double R8PBU_ZERO[(MU+1)*N], the R8PBU matrix.
        //
    {
        int j;

        double[] a = new double[(mu + 1) * n];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < mu + 1; i++)
            {
                a[i + j * (mu + 1)] = 0.0;
            }
        }

        return a;
    }
}