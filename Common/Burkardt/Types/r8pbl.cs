using System;
using Burkardt.Uniform;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double[] r8pbl_dif2(int n, int ml)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PBL_DIF2 returns the DIF2 matrix in R8PBL format.
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
        //    21 July 2016
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
        //    Input, int N, the number of rows and columns.
        //
        //    Input, int ML, the number of subdiagonals.
        //    ML must be at least 0, and no more than N-1.
        //
        //    Output, double R8PBL_DIF2[(ML+1)*N], the matrix.
        //
    {
        int j;

        double[] a = r8vec_zeros_new((ml + 1) * n);

        for (j = 0; j < n; j++)
        {
            a[0 + j * (ml + 1)] = 2.0;
        }

        for (j = 0; j < n - 1; j++)
        {
            a[1 + j * (ml + 1)] = -1.0;
        }

        return a;
    }

    public static double[] r8pbl_indicator(int n, int ml)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PBL_INDICATOR sets up an R8PBL indicator matrix.
        //
        //  Discussion:
        //
        //    The R8PBL storage format is used for a symmetric positive definite band matrix.
        //
        //    To save storage, only the diagonal and lower triangle of A is stored,
        //    in a compact diagonal format that preserves columns.
        //
        //    The diagonal is stored in row 1 of the array.
        //    The first subdiagonal in row 2, columns 1 through ML.
        //    The second subdiagonal in row 3, columns 1 through ML-1.
        //    The ML-th subdiagonal in row ML+1, columns 1 through 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 January 2004
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
        //    Input, int ML, the number of subdiagonals in the matrix.
        //    ML must be at least 0 and no more than N-1.
        //
        //    Output, double R8PBL_INDICATOR[(ML+1)*N], the R8PBL matrix.
        //
    {
        int i;

        double[] a = r8vec_zeros_new((ml + 1) * n);

        int fac = (int) Math.Pow(10, (int) Math.Log10(n) + 1);

        for (i = 0; i <= n; i++)
        {
            int j;
            for (j = Math.Max(1, i - ml); j <= i; j++)
            {
                a[i - j + (j - 1) * (ml + 1)] = fac * (i + 1) + j + 1;
            }
        }

        return a;
    }

    public static double[] r8pbl_mv(int n, int ml, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PBL_MV multiplies an R8PBL matrix by an R8VEC.
        //
        //  Discussion:
        //
        //    The R8PBL storage format is for a symmetric positive definite band matrix.
        //
        //    To save storage, only the diagonal and lower triangle of A is stored,
        //    in a compact diagonal format that preserves columns.
        //
        //    The diagonal is stored in row 1 of the array.
        //    The first subdiagonal in row 2, columns 1 through ML.
        //    The second subdiagonal in row 3, columns 1 through ML-1.
        //    The ML-th subdiagonal in row ML+1, columns 1 through 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 July 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int ML, the number of subdiagonals in the matrix.
        //    ML must be at least 0 and no more than N-1.
        //
        //    Input, double A([ML+1)*N], the R8PBL matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double R8PBL_MV[M], the result vector A * x.
        //
    {
        int j;
        int k;

        double[] b = r8vec_zeros_new(n);
        //
        //  Multiply X by the diagonal of the matrix.
        //
        for (j = 0; j < n; j++)
        {
            b[j] = a[0 + j * (ml + 1)] * x[j];
        }

        //
        //  Multiply X by the subdiagonals of the matrix.
        //
        for (k = 0; k < ml; k++)
        {
            for (j = 0; j < n - k; j++)
            {
                int i = j + k;
                double aij = a[k + 1 + j * (ml + 1)];
                b[i] += aij * x[j];
                b[j] += aij * x[i];
            }
        }

        return b;
    }

    public static void r8pbl_print(int n, int ml, double[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PBL_PRINT prints an R8PBL matrix.
        //
        //  Discussion:
        //
        //    The R8PBL storage format is used for a symmetric positive definite band matrix.
        //
        //    To save storage, only the diagonal and lower triangle of A is stored,
        //    in a compact diagonal format that preserves columns.
        //
        //    The diagonal is stored in row 1 of the array.
        //    The first subdiagonal in row 2, columns 1 through ML.
        //    The second subdiagonal in row 3, columns 1 through ML-1.
        //    The ML-th subdiagonal in row ML+1, columns 1 through 1.
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
        //    Input, int ML, the upper (and lower) bandwidth.
        //    ML must be nonnegative, and no greater than N-1.
        //
        //    Input, double A[(ML+1)*N], the R8PBL matrix.
        //
        //    Input, string TITLE, a title.
        //
    {
        r8pbl_print_some(n, ml, a, 0, 0, n - 1, n - 1, title);
    }

    public static void r8pbl_print_some(int n, int ml, double[] a, int ilo, int jlo, int ihi,
            int jhi, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PBL_PRINT_SOME prints some of an R8PBL matrix.
        //
        //  Discussion:
        //
        //    The R8PBL storage format is used for a symmetric positive definite band matrix.
        //
        //    To save storage, only the diagonal and lower triangle of A is stored,
        //    in a compact diagonal format that preserves columns.
        //
        //    The diagonal is stored in row 1 of the array.
        //    The first subdiagonal in row 2, columns 1 through ML.
        //    The second subdiagonal in row 3, columns 1 through ML-1.
        //    The ML-th subdiagonal in row ML+1, columns 1 through 1.
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
        //    Input, int ML, the upper (and lower) bandwidth.
        //    ML must be nonnegative, and no greater than N-1.
        //
        //    Input, double A[(ML+1)*N], the R8PBL matrix.
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
            j2hi = Math.Min(j2hi, n - 1);
            j2hi = Math.Min(j2hi, jhi);

            Console.WriteLine("");
            string cout = "  Col: ";
            int j;
            for (j = j2lo; j <= j2hi; j++)
            {
                cout += j.ToString().PadLeft(7) + "       ";
            }

            Console.WriteLine(cout);
            Console.WriteLine("  Row");
            Console.WriteLine("  ---");
            //
            //  Determine the range of the rows in this strip.
            //
            int i2lo = Math.Max(ilo, 0);
            i2lo = Math.Max(i2lo, j2lo - ml);

            int i2hi = Math.Min(ihi, n);
            i2hi = Math.Min(i2hi, j2hi + ml);

            int i;
            for (i = i2lo; i <= i2hi; i++)
            {
                cout = i.ToString().PadLeft(4) + "  ";
                //
                //  Print out (up to) 5 entries in row I, that lie in the current strip.
                //
                for (j = j2lo; j <= j2hi; j++)
                {
                    if (i <= j && j <= i + ml)
                    {
                        cout += a[j - i + i * (ml + 1)].ToString().PadLeft(12) + "  ";
                    }
                    else if (j <= i && i <= j + ml)
                    {
                        cout += a[i - j + j * (ml + 1)].ToString().PadLeft(12) + "  ";
                    }
                    else
                    {
                        cout += "              ";
                    }
                }

                Console.WriteLine(cout);
            }
        }
    }

    public static double[] r8pbl_random(int n, int ml, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PBL_RANDOM randomizes an R8PBL matrix.
        //
        //  Discussion:
        //
        //    The R8PBL storage format is used for a symmetric positive definite band matrix.
        //
        //    To save storage, only the diagonal and lower triangle of A is stored,
        //    in a compact diagonal format that preserves columns.
        //
        //    The diagonal is stored in row 1 of the array.
        //    The first subdiagonal in row 2, columns 1 through ML.
        //    The second subdiagonal in row 3, columns 1 through ML-1.
        //    The ML-th subdiagonal in row ML+1, columns 1 through 1.
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
        //    15 October 2003
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
        //    Input, int ML, the number of subdiagonals in the matrix.
        //    ML must be at least 0 and no more than N-1.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double R8PBL_RANDOM[(ML+1)*N], the R8PBL matrix.
        //
    {
        int i;
        int j;

        double[] a = r8vec_zeros_new((ml + 1) * n);

        for (i = 0; i < n; i++)
        {
            for (j = Math.Max(0, i - ml); j <= i - 1; j++)
            {
                a[i - j + j * (ml + 1)] = UniformRNG.r8_uniform_01(ref seed);
            }
        }

        //
        //  Set the diagonal values.
        //
        for (i = 0; i < n; i++)
        {
            double sum2 = 0.0;

            for (j = Math.Max(0, i - ml); j <= i - 1; j++)
            {
                sum2 += Math.Abs(a[i - j + j * (ml + 1)]);
            }

            for (j = i + 1; j <= Math.Min(i + ml, n - 1); j++)
            {
                sum2 += Math.Abs(a[j - i + i * (ml + 1)]);
            }

            double r = UniformRNG.r8_uniform_01(ref seed);

            a[0 + i * (ml + 1)] = (1.0 + r) * (sum2 + 0.01);
        }

        return a;
    }

    public static double[] r8pbl_to_r8ge(int n, int ml, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PBL_TO_R8GE copies an R8PBL matrix to an R8GE matrix.
        //
        //  Discussion:
        //
        //    The R8PBL storage format is used for a symmetric positive definite band matrix.
        //
        //    To save storage, only the diagonal and lower triangle of A is stored,
        //    in a compact diagonal format that preserves columns.
        //
        //    The diagonal is stored in row 1 of the array.
        //    The first subdiagonal in row 2, columns 1 through ML.
        //    The second subdiagonal in row 3, columns 1 through ML-1.
        //    The ML-th subdiagonal in row ML+1, columns 1 through 1.
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
        //    Input, int ML, the upper bandwidth of A1.
        //    ML must be nonnegative, and no greater than N-1.
        //
        //    Input, double A[(ML+1)*N], the R8PBL matrix.
        //
        //    Output, double R8PBL_TO_R8GE[N*N], the R8GE matrix.
        //
    {
        int i;

        double[] b = r8vec_zeros_new(n * n);

        for (i = 0; i < n; i++)
        {
            int j;
            for (j = 0; j < n; j++)
            {
                if (i <= j && j <= i + ml)
                {
                    b[i + j * n] = a[j - i + i * (ml + 1)];
                }
                else if (i - ml <= j && j < i)
                {
                    b[i + j * n] = a[i - j + j * (ml + 1)];
                }
            }
        }

        return b;
    }

    public static double[] r8pbl_zeros(int n, int ml)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PBL_ZEROS zeros an R8PBL matrix.
        //
        //  Discussion:
        //
        //    The R8PBL storage format is used for a symmetric positive definite band matrix.
        //
        //    To save storage, only the diagonal and lower triangle of A is stored,
        //    in a compact diagonal format that preserves columns.
        //
        //    The diagonal is stored in row 1 of the array.
        //    The first subdiagonal in row 2, columns 1 through ML.
        //    The second subdiagonal in row 3, columns 1 through ML-1.
        //    The ML-th subdiagonal in row ML+1, columns 1 through 1.
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
        //    Input, int ML, the number of subdiagonals in the matrix.
        //    ML must be at least 0 and no more than N-1.
        //
        //    Output, double R8PBL_ZERO[(ML+1)*N], the R8PBL matrix.
        //
    {
        double[] a = r8vec_zeros_new((ml + 1) * n);

        return a;
    }
}