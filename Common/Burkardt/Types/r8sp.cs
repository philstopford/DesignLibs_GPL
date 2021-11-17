using System;
using System.IO;
using Burkardt.Uniform;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r8sp_cg(int n, int nz_num, int[] row, int[] col, double[] a,
            double[] b, ref double[] x)

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

    public static bool r8sp_check(int m, int n, int nz_num, int[] row, int[] col)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SP_CHECK checks that an R8SP matrix data structure is properly sorted.
        //
        //  Discussion:
        //
        //    This routine assumes that the data structure has been sorted,
        //    so that the entries of ROW are ascending sorted, and that the
        //    entries of COL are ascending sorted, within the group of entries
        //    that have a common value of ROW.
        //
        //    The R8SP storage format stores the row, column and value of each nonzero
        //    entry of a sparse matrix.
        //
        //    The R8SP format is used by CSPARSE ("sparse triplet"), SLAP
        //    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns of
        //    the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero elements in
        //    the matrix.
        //
        //    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and
        //    column indices of the nonzero elements.
        //
        //    Output, bool R8SP_CHECK, is TRUE if the matrix is properly defined.
        //
    {
        bool check;
        int k;

        check = true;
        //
        //  Check 1 <= ROW(*) <= M.
        //
        for (k = 0; k < nz_num; k++)
        {
            if (row[k] < 0 || m <= row[k])
            {
                check = false;
                return check;
            }
        }

        //
        //  Check 1 <= COL(*) <= N.
        //
        for (k = 0; k < nz_num; k++)
        {
            if (col[k] < 0 || n <= col[k])
            {
                check = false;
                return check;
            }
        }

        //
        //  Check that ROW(K) <= ROW(K+1).
        //
        for (k = 0; k < nz_num - 1; k++)
        {
            if (row[k + 1] < row[k])
            {
                check = false;
                return check;
            }
        }

        //
        //  Check that, if ROW(K) == ROW(K+1), that COL(K) < COL(K+1).
        //
        for (k = 0; k < nz_num - 1; k++)
        {
            if (row[k] == row[k + 1])
            {
                if (col[k + 1] <= col[k])
                {
                    check = false;
                    return check;
                }
            }
        }

        return check;
    }

    public static void r8sp_diagonal(int m, int n, int nz_num, int[] row, int[] col, ref double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SP_DIAGONAL reorders an R8SP matrix so diagonal entries are first.
        //
        //  Discussion:
        //
        //    The R8SP storage format corresponds to the SLAP Triad format.
        //
        //    The R8SP storage format stores the row, column and value of each nonzero
        //    entry of a sparse matrix.  The entries may be given in any order.  No
        //    check is made for the erroneous case in which a given matrix entry is
        //    specified more than once.
        //
        //    This routine reorders the entries of A so that the first N entries
        //    are exactly the diagonal entries of the matrix, in order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 September 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero elements in 
        //    the matrix.
        //
        //    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and 
        //    column indices of the nonzero elements.
        //
        //    Input/output, double A[NZ_NUM], the nonzero elements 
        //    of the matrix.
        //
    {
        int found;
        int i;
        int j;
        int k;
        double t;

        found = 0;

        for (k = 0; k < nz_num; k++)
        {
            while (row[k] == col[k])
            {
                if (row[k] == k)
                {
                    found += 1;
                    break;
                }

                i = row[k];

                j = row[i];
                row[i] = row[k];
                row[k] = j;

                j = col[i];
                col[i] = col[k];
                col[k] = j;

                t = a[i];
                a[i] = a[k];
                a[k] = t;

                found += 1;

                if (Math.Min(m, n) <= found)
                {
                    break;
                }
            }

            if (Math.Min(m, n) <= found)
            {
                break;
            }
        }

        if (found < Math.Min(m, n))
        {
            Console.WriteLine("");
            Console.WriteLine("R8SP_DIAGONAL - Warning!");
            Console.WriteLine("  Number of diagonal entries expected: " + Math.Min(m, n) + "");
            Console.WriteLine("  Number found was " + found + "");
        }
    }

    public static double[] r8sp_dif2(int m, int n, int nz_num, int[] row, int[] col)

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
            switch (i)
            {
                case > 0:
                    row[k] = i;
                    col[k] = i - 1;
                    a[k] = -1.0;
                    k += 1;
                    break;
            }

            row[k] = i;
            col[k] = i;
            a[k] = 2.0;
            k += 1;

            if (i < n - 1)
            {
                row[k] = i;
                col[k] = i + 1;
                a[k] = -1.0;
                k += 1;
            }
        }

        return a;
    }

    public static int r8sp_ij_to_k(int nz_num, int[] row, int[] col, int i, int j)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SP_IJ_TO_K seeks the compressed index of the (I,J) entry of A.
        //
        //  Discussion:
        //
        //    If A(I,J) is nonzero, then its value is stored in location K.
        //
        //    This routine searches the R8SP storage structure for the index K
        //    corresponding to (I,J), returning -1 if no such entry was found.
        //
        //    This routine assumes that the data structure has been sorted,
        //    so that the entries of ROW are ascending sorted, and that the
        //    entries of COL are ascending sorted, within the group of entries
        //    that have a common value of ROW.
        //
        //    The R8SP storage format stores the row, column and value of each nonzero
        //    entry of a sparse matrix.
        //
        //    The R8SP format is used by CSPARSE ("sparse triplet"), SLAP
        //    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NZ_NUM, the number of nonzero elements in
        //    the matrix.
        //
        //    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and
        //    column indices of the nonzero elements.
        //
        //    Input, int I, J, the row and column indices of the
        //    matrix entry.
        //
        //    Output, int R8SP_IJ_TO_K, the R8SP index of the (I,J) entry.
        //
    {
        int hi;
        int k;
        int lo;
        int md;

        lo = 0;
        hi = nz_num - 1;

        for (;;)
        {
            if (hi < lo)
            {
                k = -1;
                break;
            }

            md = (lo + hi) / 2;

            if (row[md] < i || row[md] == i && col[md] < j)
            {
                lo = md + 1;
            }
            else if (i < row[md] || row[md] == i && j < col[md])
            {
                hi = md - 1;
            }
            else
            {
                k = md;
                break;
            }
        }

        return k;
    }

    public static double[] r8sp_indicator(int m, int n, int nz_num, int[] row, int[] col)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SP_INDICATOR sets up an R8SP indicator matrix.
        //
        //  Discussion:
        //
        //    The R8SP storage format stores the row, column and value of each nonzero
        //    entry of a sparse matrix.
        //
        //    It is possible that a pair of indices (I,J) may occur more than
        //    once.  Presumably, in this case, the intent is that the actual value
        //    of A(I,J) is the sum of all such entries.  This is not a good thing
        //    to do, but I seem to have come across this in MATLAB.
        //
        //    The R8SP format is used by CSPARSE ("sparse triplet"), SLAP
        //    (nonsymmetric case), by MATLAB, and by SPARSEKIT ("COO" format).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 September 2015
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
        //    Output, double R8SP_INDICATOR[NZ_NUM], the nonzero elements of the matrix.
        //
    {
        double[] a;
        int fac;
        int i;
        int j;
        int k;

        a = new double[nz_num];

        fac = (int) Math.Pow(10, (int) Math.Log10(n) + 1);

        for (k = 0; k < nz_num; k++)
        {
            i = row[k];
            j = col[k];
            a[k] = fac * (i + 1) + j + 1;
        }

        return a;
    }

    public static void r8sp_jac_sl(int n, int nz_num, int[] row, int[] col, double[] a,
            double[] b, ref double[] x, int it_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SP_JAC_SL solves an R8SP system using Jacobi iteration.
        //
        //  Discussion:
        //
        //    The R8SP storage format corresponds to the SLAP Triad format.
        //
        //    The R8SP storage format stores the row, column and value of each nonzero
        //    entry of a sparse matrix.  The entries may be given in any order.  No
        //    check is made for the erroneous case in which a given matrix entry is
        //    specified more than once.
        //
        //    This routine REQUIRES that the matrix be square, that the matrix
        //    have nonzero diagonal entries, and that the first N entries of
        //    the array A be exactly the diagonal entries of the matrix, in order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 September 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero elements in 
        //    the matrix.
        //
        //    Input, int ROW[NZ_NUM], COL[NZ_NUM], the row and column 
        //    indices of the nonzero elements.
        //
        //    Input, double A[NZ_NUM], the nonzero elements of the matrix.
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
        int j;
        int k;
        double[] x_new;

        r8sp_diagonal(n, n, nz_num, row, col, ref a);

        x_new = new double[n];

        for (it_num = 1; it_num <= it_max; it_num++)
        {
            //
            //  Initialize to right hand side.
            //
            for (j = 0; j < n; j++)
            {
                x_new[j] = b[j];
            }

            //
            //  Subtract off-diagonal terms.
            //
            for (k = n; k < nz_num; k++)
            {
                i = row[k];
                j = col[k];
                x_new[i] -= a[k] * x[j];
            }

            //
            //  Divide by diagonal terms and update.
            //
            for (j = 0; j < n; j++)
            {
                x[j] = x_new[j] / a[j];
            }
        }
    }

    public static double[] r8sp_mtv(int m, int n, int nz_num, int[] row, int[] col,
            double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SP_MTV multiplies a vector times an R8SP matrix.
        //
        //  Discussion:
        //
        //    The R8SP storage format stores the row, column and value of each nonzero
        //    entry of a sparse matrix.
        //
        //    It is possible that a pair of indices (I,J) may occur more than
        //    once.  Presumably, in this case, the intent is that the actual value
        //    of A(I,J) is the sum of all such entries.  This is not a good thing
        //    to do, but I seem to have come across this in MATLAB.
        //
        //    The R8SP format is used by CSPARSE ("sparse triplet"), SLAP
        //    (nonsymmetric case), by MATLAB, and by SPARSEKIT ("COO" format).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 October 2004
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
        //    Input, double X[M], the vector to be multiplied by A.
        //
        //    Output, double B[N], the product vector A'*X.
        //
    {
        double[] b;
        int i;
        int j;
        int k;

        b = r8vec_zeros_new(n);

        for (k = 0; k < nz_num; k++)
        {
            i = row[k];
            j = col[k];
            b[j] += a[k] * x[i];
        }

        return b;
    }

    public static double[] r8sp_mv(int m, int n, int nz_num, int[] row, int[] col,
            double[] a, double[] x)

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
            b[i] += a[k] * x[j];
        }

        return b;
    }

    public static void r8sp_print(int m, int n, int nz_num, int[] row, int[] col,
            double[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SP_PRINT prints an R8SP matrix.
        //
        //  Discussion:
        //
        //    This version of R8SP_PRINT has been specifically modified to allow,
        //    and correctly handle, the case in which a single matrix location
        //    A(I,J) is referenced more than once by the sparse matrix structure.
        //    In such cases, the routine prints out the sum of all the values.
        //
        //    The R8SP storage format stores the row, column and value of each nonzero
        //    entry of a sparse matrix.
        //
        //    It is possible that a pair of indices (I,J) may occur more than
        //    once.  Presumably, in this case, the intent is that the actual value
        //    of A(I,J) is the sum of all such entries.  This is not a good thing
        //    to do, but I seem to have come across this in MATLAB.
        //
        //    The R8SP format is used by CSPARSE ("sparse triplet"), SLAP
        //    (nonsymmetric case), by MATLAB, and by SPARSEKIT ("COO" format).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 September 2015
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
        //    Input, string TITLE, a title.
        //
    {
        r8sp_print_some(m, n, nz_num, row, col, a, 0, 0, m - 1, n - 1, title);

    }

    public static void r8sp_print_some(int m, int n, int nz_num, int[] row, int[] col,
            double[] a, int ilo, int jlo, int ihi, int jhi, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SP_PRINT_SOME prints some of an R8SP matrix.
        //
        //  Discussion:
        //
        //    This version of R8SP_PRINT_SOME has been specifically modified to allow,
        //    and correctly handle, the case in which a single matrix location
        //    A(I,J) is referenced more than once by the sparse matrix structure.
        //    In such cases, the routine prints out the sum of all the values.
        //
        //    The R8SP storage format stores the row, column and value of each nonzero
        //    entry of a sparse matrix.
        //
        //    It is possible that a pair of indices (I,J) may occur more than
        //    once.  Presumably, in this case, the intent is that the actual value
        //    of A(I,J) is the sum of all such entries.  This is not a good thing
        //    to do, but I seem to have come across this in MATLAB.
        //
        //    The R8SP format is used by CSPARSE ("sparse triplet"), SLAP
        //    (nonsymmetric case), by MATLAB, and by SPARSEKIT ("COO" format).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 September 2015
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
        //    Input, int ILO, JLO, IHI, JHI, the first row and
        //    column, and the last row and column to be printed.
        //
        //    Input, string TITLE, a title.
        //
    {
        int INCX = 5;

        double[] aij = new double[INCX];
        int i;
        int i2hi;
        int i2lo;
        int inc;
        int j;
        int j2;
        int j2hi;
        int j2lo;
        int k;
        bool nonzero;
        string cout = "";

        Console.WriteLine("");
        Console.WriteLine(title + "");
        //
        //  Print the columns of the matrix, in strips of 5.
        //
        for (j2lo = jlo; j2lo <= jhi; j2lo += INCX)
        {
            j2hi = j2lo + INCX - 1;
            j2hi = Math.Min(j2hi, n - 1);
            j2hi = Math.Min(j2hi, jhi);

            inc = j2hi + 1 - j2lo;

            Console.WriteLine("");

            cout = "  Col:  ";
            for (j = j2lo; j <= j2hi; j++)
            {
                cout +=  j.ToString().PadLeft(7) + "       ";
            }

            Console.WriteLine(cout);
            Console.WriteLine("  Row");
            Console.WriteLine("  ---");
            //
            //  Determine the range of the rows in this strip.
            //
            i2lo = Math.Max(ilo, 0);
            i2hi = Math.Min(ihi, m - 1);

            for (i = i2lo; i <= i2hi; i++)
            {
                //
                //  Print out (up to) 5 entries in row I, that lie in the current strip.
                //
                nonzero = false;
                for (j2 = 0; j2 < INCX; j2++)
                {
                    aij[j2] = 0.0;
                }

                for (k = 1; k <= nz_num; k++)
                {
                    if (i == row[k - 1] && j2lo <= col[k - 1] && col[k - 1] <= j2hi)
                    {
                        j2 = col[k - 1] - j2lo;

                        switch (a[k - 1])
                        {
                            case 0.0:
                                continue;
                            default:
                                nonzero = true;
                                aij[j2] += a[k - 1];
                                break;
                        }
                    }
                }

                switch (nonzero)
                {
                    case true:
                    {
                        cout = i.ToString().PadLeft(6);
                        for (j2 = 0; j2 < inc; j2++)
                        {
                            cout += aij[j2].ToString().PadLeft(12) + "  ";
                        }

                        Console.WriteLine(cout);
                        break;
                    }
                }
            }
        }

        Console.WriteLine("");

    }

    public static double[] r8sp_random(int m, int n, int nz_num, int[] row, int[] col,
            ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SP_RANDOM sets a random R8SP matrix.
        //
        //  Discussion:
        //
        //    The R8SP storage format stores the row, column and value of each nonzero
        //    entry of a sparse matrix.
        //
        //    It is possible that a pair of indices (I,J) may occur more than
        //    once.  Presumably, in this case, the intent is that the actual value
        //    of A(I,J) is the sum of all such entries.  This is not a good thing
        //    to do, but I seem to have come across this in MATLAB.
        //
        //    The R8SP format is used by CSPARSE ("sparse triplet"), SLAP
        //    (nonsymmetric case), by MATLAB, and by SPARSEKIT ("COO" format).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 January 2004
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
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double R8SP_RANDOM[NZ_NUM], the nonzero elements of the matrix.
        //
    {
        double[] r;

        r = UniformRNG.r8vec_uniform_01_new(nz_num, ref seed);

        return r;
    }

    public static void r8sp_read(string input_file, int m, int n, int nz_num, ref int[] row,
            ref int[] col, ref double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SP_READ reads an R8SP matrix from a file.
        //
        //  Discussion:
        //
        //    This routine needs the value of NZ_NUM, which can be determined
        //    by a call to R8SP_READ_SIZE.
        //
        //    The R8SP storage format stores the row, column and value of each nonzero
        //    entry of a sparse matrix.
        //
        //    It is possible that a pair of indices (I,J) may occur more than
        //    once.  Presumably, in this case, the intent is that the actual value
        //    of A(I,J) is the sum of all such entries.  This is not a good thing
        //    to do, but I seem to have come across this in MATLAB.
        //
        //    The R8SP format is used by CSPARSE ("sparse triplet"), SLAP
        //    (nonsymmetric case), by MATLAB, and by SPARSEKIT ("COO" format).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 August 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string INPUT_FILE, the name of the file to be read.
        //
        //    Unused, int M, N, the number of rows and columns of the matrix.
        //
        //    Input, int NZ_NUM, the number of nonzero elements in the matrix.
        //
        //    Output, int ROW[NZ_NUM], COL[NZ_NUM], the row and column indices
        //    of the nonzero elements.
        //
        //    Output, double A[NZ_NUM], the nonzero elements of the matrix.
        //
    {
        string[] input;
        int k = 0;

        try
        {
            input = File.ReadAllLines(input_file);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("R8SP_READ - Fatal error!");
            Console.WriteLine("  Could not open the input file: \"" + input_file + "\"");
            return;
        }

        foreach (string line in input)
        {
            string[] tokens = Helpers.splitStringByWhitespace(line);
            row[k] = Convert.ToInt32(tokens[0]);
            col[k] = Convert.ToInt32(tokens[1]);
            a[k] = Convert.ToInt32(tokens[2]);
            k++;
        }
    }

    public static void r8sp_read_size(string input_file, ref int m, ref int n, ref int nz_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SP_READ_SIZE reads the size of an R8SP matrix from a file.
        //
        //  Discussion:
        //
        //    The value of NZ_NUM is simply the number of records in the input file.
        //
        //    The values of M and N are determined as the maximum entry in the row 
        //    and column vectors.
        //
        //    The R8SP storage format stores the row, column and value of each nonzero
        //    entry of a sparse matrix.
        //
        //    It is possible that a pair of indices (I,J) may occur more than
        //    once.  Presumably, in this case, the intent is that the actual value
        //    of A(I,J) is the sum of all such entries.  This is not a good thing
        //    to do, but I seem to have come across this in MATLAB.
        //
        //    The R8SP format is used by CSPARSE ("sparse triplet"), SLAP
        //    (nonsymmetric case), by MATLAB, and by SPARSEKIT ("COO" format).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 August 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string INPUT_FILE, the name of the file to 
        //    be read.
        //
        //    Output, int &M, &N, the number of rows and columns of the matrix.
        //
        //    Output, int &NZ_NUM, the number of nonzero elements in the matrix.
        //
    {
        int col_k;
        string[] input;
        int row_k;

        m = 0;
        n = 0;
        nz_num = 0;

        try
        {
            input = File.ReadAllLines(input_file);
        }
        catch (Exception)
        {
            Console.WriteLine("");
            Console.WriteLine("R8SP_READ_SIZE - Fatal error!");
            Console.WriteLine("  Could not open the input file: \"" + input_file + "\"");
            return;
        }

        foreach (string line in input)
        {
            string[] tokens = Helpers.splitStringByWhitespace(line);
            row_k = Convert.ToInt32(tokens[0]);
            col_k = Convert.ToInt32(tokens[1]);


            m = Math.Max(m, row_k + 1);
            n = Math.Max(n, col_k + 1);
            nz_num += 1;
        }


    }

    public static double[] r8sp_res(int m, int n, int nz_num, int[] row, int[] col, double[] a,
            double[] x, double[] b)

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