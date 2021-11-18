﻿using System;
using System.Globalization;
using Burkardt.Uniform;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r8sd_cg(int n, int ndiag, int[] offset, double[] a, double[] b,
            ref double[] x)

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
        int i;
        int it;
        //
        //  Initialize
        //    AP = A * x,
        //    R  = b - A * x,
        //    P  = b - A * x.
        //
        double[] ap = r8sd_mv(n, n, ndiag, offset, a, x);

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
            double pap = r8vec_dot_product(n, p, ap);

            if (pap == 0.0)
            {
                break;
            }

            double pr = r8vec_dot_product(n, p, r);

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
            a[i] = 2.0;
        }

        for (int i = 0; i < n - 1; i++)
        {
            const int j = 1;
            a[i + j * ndiag] = -1.0;
        }

        return a;
    }

    public static double[] r8sd_indicator(int n, int ndiag, int[] offset)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SD_INDICATOR sets up an R8SD indicator matrix.
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
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 January 2004
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
        //    Input, int NDIAG, the number of diagonals that are stored.
        //    NDIAG must be at least 1 and no more than N.
        //
        //    Input, int OFFSET[NDIAG], the offsets for the diagonal storage.
        //
        //    Output, double R8SD_INDICATOR[N*NDIAG], the R8SD matrix.
        //
    {
        int i;

        double[] a = r8vec_zeros_new(n * ndiag);

        int fac = (int) Math.Pow(10, (int) Math.Log10(n) + 1);

        for (i = 1; i <= n; i++)
        {
            int diag;
            for (diag = 1; diag <= ndiag; diag++)
            {
                int j = i + offset[diag - 1];
                a[i - 1 + (diag - 1) * n] = j switch
                {
                    >= 1 when j <= n => fac * i + j,
                    _ => 0.0
                };
            }
        }

        return a;
    }

    public static double[] r8sd_mv(int m, int n, int ndiag, int[] offset, double[] a,
            double[] x)

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
                switch (j)
                {
                    case >= 0 when j < n:
                    {
                        b[i] += a[i + jdiag * n] * x[j];
                        if (offset[jdiag] != 0)
                        {
                            b[j] += a[i + jdiag * n] * x[i];
                        }

                        break;
                    }
                }
            }
        }

        return b;
    }

    public static void r8sd_print(int n, int ndiag, int[] offset, double[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SD_PRINT prints an R8SD matrix.
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
        //    06 April 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of columns of the matrix.
        //    N must be positive.
        //
        //    Input, int NDIAG, the number of diagonals of the matrix
        //    that are stored in the array.
        //    NDIAG must be at least 1, and no more than N.
        //
        //    Input, int OFFSET[NDIAG], the offsets for the diagonal storage.
        //
        //    Input, double A[N*NDIAG], the R8SD matrix.
        //
        //    Input, string TITLE, a title.
        //
    {
        r8sd_print_some(n, ndiag, offset, a, 1, 1, n, n, title);

    }

    public static void r8sd_print_some(int n, int ndiag, int[] offset, double[] a, int ilo,
            int jlo, int ihi, int jhi, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SD_PRINT_SOME prints some of an R8SD matrix.
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
        //    06 April 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of columns of the matrix.
        //    N must be positive.
        //
        //    Input, int NDIAG, the number of diagonals of the matrix
        //    that are stored in the array.
        //    NDIAG must be at least 1, and no more than N.
        //
        //    Input, int OFFSET[NDIAG], the offsets for the diagonal storage.
        //
        //    Input, double A[N*NDIAG], the R8SD matrix.
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
            int i2hi = Math.Min(ihi, n);

            int i;
            for (i = i2lo; i <= i2hi; i++)
            {
                cout = i.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  ";
                //
                //  Print out (up to) 5 entries in row I, that lie in the current strip.
                //
                for (j = j2lo; j <= j2hi; j++)
                {
                    double aij = 0.0;

                    int jdiag;
                    for (jdiag = 0; jdiag < ndiag; jdiag++)
                    {
                        if (j - i == offset[jdiag])
                        {
                            aij = a[i - 1 + jdiag * n];
                        }
                        else if (j - i == -offset[jdiag])
                        {
                            aij = a[j - 1 + jdiag * n];
                        }
                    }

                    cout += aij.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  ";
                }

                Console.WriteLine(cout);
            }
        }
    }

    public static double[] r8sd_random(int n, int ndiag, int[] offset, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SD_RANDOM randomizes an R8SD matrix.
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
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 January 2004
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
        //    Input, int NDIAG, the number of diagonals that are stored.
        //    NDIAG must be at least 1 and no more than N.
        //
        //    Input, int OFFSET[NDIAG], the offsets for the diagonal storage.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double R8SD_RANDOM[N*NDIAG], the R8SD matrix.
        //
    {
        int i;

        double[] a = r8vec_zeros_new(n * ndiag);

        for (i = 1; i <= n; i++)
        {
            int j;
            for (j = 1; j <= ndiag; j++)
            {
                int jj = i + offset[j - 1];
                a[i - 1 + (j - 1) * n] = jj switch
                {
                    >= 1 when jj <= n => UniformRNG.r8_uniform_01(ref seed),
                    _ => a[i - 1 + (j - 1) * n]
                };
            }
        }

        return a;
    }

    public static double[] r8sd_res(int m, int n, int ndiag, int[] offset, double[] a,
            double[] x, double[] b)

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

    public static double[] r8sd_to_r8ge(int n, int ndiag, int[] offset, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SD_TO_R8GE copies an R8SD matrix to an R8GE matrix.
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
        //    18 October 2003
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
        //    Input, int NDIAG, the number of diagonals that are stored.
        //    NDIAG must be at least 1 and no more than N.
        //
        //    Input, int OFFSET[NDIAG], the offsets for the diagonal storage.
        //
        //    Input, double A[N*NDIAG], the R8SD matrix.
        //
        //    Output, double R8SD_TO_R8GE[N*N], the R8GE matrix.
        //
    {
        int i;
        int j;

        double[] b = r8vec_zeros_new(n * n);

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                b[i + j * n] = 0.0;
            }
        }

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < ndiag; j++)
            {
                int jj = i + offset[j];
                switch (jj)
                {
                    case >= 0 when jj <= n - 1:
                    {
                        b[i + jj * n] = a[i + j * n];
                        if (i != jj)
                        {
                            b[jj + i * n] = a[i + j * n];
                        }

                        break;
                    }
                }
            }
        }

        return b;
    }

    public static double[] r8sd_zeros(int n, int ndiag, int[] offset)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8SD_ZEROS zeros an R8SD matrix.
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
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 July 2017
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
        //    Input, int NDIAG, the number of diagonals that are stored.
        //    NDIAG must be at least 1 and no more than N.
        //
        //    Input, int OFFSET[NDIAG], the offsets for the diagonal storage.
        //
        //    Output, double R8SD_ZERO[N*NDIAG], the R8SD matrix.
        //
    {
        double[] a = r8vec_zeros_new(n * ndiag);

        return a;
    }

}