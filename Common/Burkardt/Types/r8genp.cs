using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double r8ge_np_det(int n, double[] a_lu)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8GE_NP_DET computes the determinant of a matrix factored by R8GE_NP_FA.
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
            //    14 October 2003
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
            //    Input, double A_LU[N*N], the LU factors from R8GE_NP_FA.
            //
            //    Output, double R8GE_NP_DET, the determinant of the matrix.
            //
        {
            double det;
            int i;

            det = 1.0;
            for (i = 0; i < n; i++)
            {
                det = det * a_lu[i + i * n];
            }

            return det;
        }

        public static int r8ge_np_fa(int n, ref double[] a, int aIndex = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8GE_NP_FA factors an R8GE matrix by nonpivoting Gaussian elimination.
            //
            //  Discussion:
            //
            //    The R8GE storage format is used for a "general" M by N matrix.  
            //    A physical storage space is made for each logical entry.  The two 
            //    dimensional logical array is mapped to a vector, in which storage is 
            //    by columns.
            //
            //    R8GE_NP_FA is a version of the LINPACK routine SGEFA, but uses no
            //    pivoting.  It will fail if the matrix is singular, or if any zero
            //    pivot is encountered.
            //
            //    If R8GE_NP_FA successfully factors the matrix, R8GE_NP_SL may be called
            //    to solve linear systems involving the matrix.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 October 2003
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
            //    Input/output, double A[N*N].
            //    On input, A contains the matrix to be factored.
            //    On output, A contains information about the factorization,
            //    which must be passed unchanged to R8GE_NP_SL for solutions.
            //
            //    Output, int R8GE_NP_FA, singularity flag.
            //    0, no singularity detected.
            //    nonzero, the factorization failed on the INFO-th step.
            //
        {
            int i;
            int j;
            int k;

            for (k = 1; k <= n - 1; k++)
            {
                if (a[((k - 1 + (k - 1) * n) + aIndex ) % a.Length] == 0.0)
                {
                    return k;
                }

                for (i = k + 1; i <= n; i++)
                {
                    a[((i - 1 + (k - 1) * n) + aIndex ) % a.Length] = -a[((i - 1 + (k - 1) * n) + aIndex ) % a.Length] / a[((k - 1 + (k - 1) * n) + aIndex ) % a.Length];
                }

                for (j = k + 1; j <= n; j++)
                {
                    for (i = k + 1; i <= n; i++)
                    {
                        a[((i - 1 + (j - 1) * n) + aIndex ) % a.Length] =
                            a[((i - 1 + (j - 1) * n) + aIndex ) % a.Length] + a[((i - 1 + (k - 1) * n) + aIndex ) % a.Length] * a[((k - 1 + (j - 1) * n) + aIndex ) % a.Length];
                    }
                }
            }

            if (a[((n - 1 + (n - 1) * n) + aIndex ) % a.Length] == 0.0)
            {
                return n;
            }

            return 0;
        }

        public static double[] r8ge_np_inverse(int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8GE_NP_INVERSE computes the inverse of a matrix factored by R8GE_NP_FA.
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
            //    02 November 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the matrix A.
            //
            //    Input, double A[N*N], the factor information computed by R8GE_NP_FA.
            //
            //    Output, double R8GE_NP_INVERSE[N*N], the inverse matrix.
            //
        {
            double[] b;
            int i;
            int j;
            int k;
            double temp;
            double[] work;

            b = new double[n * n];
            work = new double[n];

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < n; i++)
                {
                    b[i + j * n] = a[i + j * n];
                }
            }

            //
            //  Compute Inverse(U).
            //
            for (k = 1; k <= n; k++)
            {
                b[k - 1 + (k - 1) * n] = 1.0 / b[k - 1 + (k - 1) * n];
                for (i = 1; i <= k - 1; i++)
                {
                    b[i - 1 + (k - 1) * n] = -b[i - 1 + (k - 1) * n] * b[k - 1 + (k - 1) * n];
                }

                for (j = k + 1; j <= n; j++)
                {
                    temp = b[k - 1 + (j - 1) * n];
                    b[k - 1 + (j - 1) * n] = 0.0;
                    for (i = 1; i <= k; i++)
                    {
                        b[i - 1 + (j - 1) * n] = b[i - 1 + (j - 1) * n] + temp * b[i - 1 + (k - 1) * n];
                    }
                }
            }

            //
            //  Form Inverse(U) * Inverse(L).
            //
            for (k = n - 1; 1 <= k; k--)
            {
                for (i = k + 1; i <= n; i++)
                {
                    work[i - 1] = b[i - 1 + (k - 1) * n];
                    b[i - 1 + (k - 1) * n] = 0.0;
                }

                for (j = k + 1; j <= n; j++)
                {
                    for (i = 1; i <= n; i++)
                    {
                        b[i - 1 + (k - 1) * n] = b[i - 1 + (k - 1) * n] + b[i - 1 + (j - 1) * n] * work[j - 1];
                    }
                }
            }

            return b;
        }

        public static double[] r8ge_np_ml(int n, double[] a_lu, double[] x, int job)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8GE_NP_ML computes A * x or x * A, for a matrix factored by R8GE_NP_FA.
            //
            //  Discussion:
            //
            //    The R8GE storage format is used for a "general" M by N matrix.  
            //    A physical storage space is made for each logical entry.  The two 
            //    dimensional logical array is mapped to a vector, in which storage is 
            //    by columns.
            //
            //    The matrix A is assumed to have been factored by R8GE_NP_FA.
            //
            //    R8GE_NP_ML allows the user to check that the solution of a linear
            //    system is correct, without having to save an unfactored copy
            //    of the matrix.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    29 January 2004
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
            //    Input, double A_LU[N*N], the LU factors from R8GE_NP_FA.
            //
            //    Input, double X[N], the vector to be multiplied.
            //
            //    Input, int JOB, determines the multiplication to
            //    be carried out:
            //    JOB = 0, compute A * x.
            //    JOB nonzero, compute A' * X.
            //
            //    Output, double R8GE_NP_ML[N], the result of the multiplication.
            //
        {
            double[] b;
            int i;
            int j;
            double t;

            b = new double[n];

            for (i = 0; i < n; i++)
            {
                b[i] = x[i];
            }

            if (job == 0)
            {
                //
                //  Compute U * X = Y:
                //
                for (i = 0; i < n; i++)
                {
                    t = 0.0;
                    for (j = i; j < n; j++)
                    {
                        t = t + a_lu[i + j * n] * b[j];
                    }

                    b[i] = t;
                }

                //
                //  Compute L * Y = B:
                //
                for (j = n - 2; 0 <= j; j--)
                {
                    for (i = j + 1; i < n; i++)
                    {
                        b[i] = b[i] - a_lu[i + j * n] * b[j];
                    }
                }
            }
            else
            {
                //
                //  Compute L' * X = Y:
                //
                for (j = 0; j < n - 1; j++)
                {
                    for (i = j + 1; i < n; i++)
                    {
                        b[j] = b[j] - b[i] * a_lu[i + j * n];
                    }
                }

                //
                //  Compute U' * Y = B:
                //
                for (j = n - 1; 0 <= j; j--)
                {
                    b[j] = b[j] * a_lu[j + j * n];
                    for (i = 0; i < j; i++)
                    {
                        b[j] = b[j] + b[i] * a_lu[i + j * n];
                    }
                }
            }

            return b;
        }

        public static double[] r8ge_np_sl(int n, double[] a_lu, double[] b, int job, int aluIndex = 0, int bIndex = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8GE_NP_SL solves a system factored by R8GE_NP_FA.
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
            //    01 November 2003
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
            //    Input, double A_LU[N*N], the LU factors from R8GE_NP_FA.
            //
            //    Input, double B[N], the right hand side.
            //
            //    Input, int JOB.
            //    If JOB is zero, the routine will solve A * x = b.
            //    If JOB is nonzero, the routine will solve A' * x = b.
            //
            //    Output, double R8GE_NP_SL[N], the solution.
            //
        {
            int i;
            int k;
            double[] x;
            //
            //  Solve A * x = b.
            //
            x = new double[n];
            for (i = 0; i < n; i++)
            {
                x[i] = b[(i + bIndex) % b.Length];
            }

            if (job == 0)
            {
                for (k = 0; k < n - 1; k++)
                {
                    for (i = k + 1; i < n; i++)
                    {
                        x[i] = x[i] + a_lu[((i + k * n) + aluIndex) % a_lu.Length] * x[k];
                    }
                }

                for (k = n - 1; 0 <= k; k--)
                {
                    x[k] = x[k] / a_lu[((k + k * n) + aluIndex) % a_lu.Length];
                    for (i = 0; i <= k - 1; i++)
                    {
                        x[i] = x[i] - a_lu[((i + k * n) + aluIndex) % a_lu.Length] * x[k];
                    }
                }
            }
            //
            //  Solve A' * X = B.
            //
            else
            {
                for (k = 0; k < n; k++)
                {
                    for (i = 0; i <= k - 1; i++)
                    {
                        x[k] = x[k] - x[i] * a_lu[((i + k * n) + aluIndex) % a_lu.Length];
                    }

                    x[k] = x[k] / a_lu[((k + k * n) + aluIndex) % a_lu.Length];
                }

                for (k = n - 2; 0 <= k; k--)
                {
                    for (i = k + 1; i < n; i++)
                    {
                        x[k] = x[k] + x[i] * a_lu[((i + k * n) + aluIndex) % a_lu.Length];
                    }
                }

            }

            return x;
        }

        public static int r8ge_np_trf(int m, int n, ref double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8GE_NP_TRF computes the LU factorization of an R8GE matrix.
            //
            //  Discussion:
            //
            //    The R8GE storage format is used for a "general" M by N matrix.  
            //    A physical storage space is made for each logical entry.  The two 
            //    dimensional logical array is mapped to a vector, in which storage is 
            //    by columns.
            //
            //    R8GE_NP_TRF is a nonpivoting version of R8GE_TRF, and will fail if
            //    a zero element is encountered along the diagonal.
            //
            //    The factorization has the form
            //      A = L * U
            //    where L is lower triangular with unit diagonal elements (lower
            //    trapezoidal if N < M), and U is upper triangular (upper trapezoidal
            //    if M < N).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    10 November 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the number of rows of the matrix A.  0 <= M.
            //
            //    Input, int N, the number of columns of the matrix A.  0 <= N.
            //
            //    Input/output, double A[M*N].
            //    On entry, the M by N matrix to be factored.
            //    On exit, the factors L and U from the factorization
            //    A = L*U; the unit diagonal elements of L are not stored.
            //
            //    Output, int R8GE_NP_TRF.
            //    = 0: successful exit
            //    = -K, the K-th argument had an illegal value
            //    = K, U(K,K) is exactly zero. The factorization
            //         has been completed, but the factor U is exactly
            //         singular, and division by zero will occur if it is used
            //         to solve a system of equations.
            //
        {
            int i;
            int ii;
            int info;
            int j;
            //
            //  Test the input parameters.
            //
            info = 0;

            if (m < 0)
            {
                return (-1);
            }
            else if (n < 0)
            {
                return (-2);
            }

            if (m == 0 || n == 0)
            {
                return 0;
            }

            for (j = 1; j <= Math.Min(m, n); j++)
            {
                //
                //  Compute elements J+1:M of the J-th column.
                //
                if (a[j - 1 + (j - 1) * m] != 0.0)
                {
                    for (i = j + 1; i <= m; i++)
                    {
                        a[i - 1 + (j - 1) * m] = a[i - 1 + (j - 1) * m] / a[j - 1 + (j - 1) * m];
                    }
                }
                else if (info == 0)
                {
                    info = j;
                }

                //
                //  Update the trailing submatrix.
                //
                if (j < Math.Min(m, n))
                {
                    for (ii = j + 1; ii <= m; ii++)
                    {
                        for (i = j + 1; i <= n; i++)
                        {
                            a[ii - 1 + (i - 1) * m] = a[ii - 1 + (i - 1) * m] -
                                                      a[ii - 1 + (j - 1) * m] * a[j - 1 + (i - 1) * m];
                        }
                    }
                }
            }

            return info;
        }

        public static double[] r8ge_np_trm(int m, int n, double[] a, double[] x, int job)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8GE_NP_TRM computes A * x or A' * x, for a matrix factored by R8GE_NP_TRF.
            //
            //  Discussion:
            //
            //    The R8GE storage format is used for a "general" M by N matrix.  
            //    A physical storage space is made for each logical entry.  The two 
            //    dimensional logical array is mapped to a vector, in which storage is 
            //    by columns.
            //
            //    The matrix A is assumed to have been factored by R8GE_NP_TRF.
            //
            //    R8GE_NP_TRM allows the user to check that the solution of a linear
            //    system is correct, without having to save an unfactored copy
            //    of the matrix.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 November 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
            //    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
            //    Sven Hammarling, Alan McKenney, Danny Sorensen,
            //    LAPACK User's Guide,
            //    Second Edition,
            //    SIAM, 1995.
            //
            //  Parameters:
            //
            //    Input, int M, N, the number of rows and columns in the matrix.
            //    M and N must be positive.
            //
            //    Input, double A[M*N], the M by N matrix factors computed by R8GE_NP_TRF.
            //
            //    Input, double X[*], the vector to be multiplied.
            //    If JOB is 0, X must have dimension N.
            //    If JOB is nonzero, X must have dimension M.
            //
            //    Input, int JOB, determines the multiplication to
            //    be carried out:
            //    JOB = 0, compute A * x.
            //    JOB nonzero, compute A' * X.
            //
            //    Output, double R8GE_NP_TRM[*], the result of the multiplication.
            //    If JOB is 0, the output has dimension M.
            //    If JOB is nonzero, the output has dimension N.
            //
        {
            double[] b;
            int i;
            int j;
            double temp;

            if (job == 0)
            {
                b = r8vec_zeros_new(m);
                //
                //  Compute U * X = Y:
                //
                for (i = 0; i < Math.Min(m, n); i++)
                {
                    for (j = i; j < n; j++)
                    {
                        b[i] = b[i] + a[i + j * m] * x[j];
                    }
                }

                //
                //  Compute L * Y = B:
                //
                for (i = Math.Min(m - 1, n); 1 <= i; i--)
                {
                    for (j = 0; j < i; j++)
                    {
                        b[i] = b[i] + a[i + j * m] * b[j];
                    }
                }
            }
            else
            {
                b = r8vec_zeros_new(n);
                //
                //  Compute L' * X = Y:
                //
                for (i = 0; i < Math.Min(m, n); i++)
                {
                    b[i] = x[i];
                    for (j = i + 1; j < m; j++)
                    {
                        b[i] = b[i] + a[j + i * m] * x[j];
                    }
                }

                //
                //  Compute U' * Y = B:
                //
                for (i = Math.Min(m, n) - 1; 0 <= i; i--)
                {
                    temp = 0.0;
                    for (j = 0; j <= i; j++)
                    {
                        temp = temp + a[j + i * m] * b[j];
                    }

                    b[i] = temp;
                }

            }

            return b;
        }

        public static double[] r8ge_np_trs(int n, int nrhs, char trans, double[] a, double[] b)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8GE_NP_TRS solves a system of linear equations factored by R8GE_NP_TRF.
            //
            //  Discussion:
            //
            //    The R8GE storage format is used for a "general" M by N matrix.  
            //    A physical storage space is made for each logical entry.  The two 
            //    dimensional logical array is mapped to a vector, in which storage is 
            //    by columns.
            //
            //    R8GE_NP_TRS is a nonpivoting version of R8GE_TRS.
            //
            //    R8GE_TRS solves a system of linear equations
            //      A * x = b  or  A' * X = B
            //    with a general N by N matrix A using the LU factorization computed
            //    by R8GE_NP_TRF.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    03 November 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
            //    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
            //    Sven Hammarling, Alan McKenney, Danny Sorensen,
            //    LAPACK User's Guide,
            //    Second Edition,
            //    SIAM, 1995.
            //
            //  Parameters:
            //
            //    Input, int N, the order of the matrix A.  0 <= N.
            //
            //    Input, int NRHS, the number of right hand sides.  0 <= NRHS.
            //
            //    Input, char TRANS, pecifies the form of the system of equations:
            //    'N':  A * x = b  (No transpose)
            //    'T':  A'* X = B  (Transpose)
            //    'C':  A'* X = B  (Conjugate transpose = Transpose)
            //
            //    Input, double A[N*N], the factors L and U from the factorization
            //    A = L*U as computed by R8GE_NP_TRF.
            //
            //    Input, double B[N*NRHS], the right hand side matrix B.
            //
            //    Output, double R8GE_NP_TRS[N*NRHS], the solution matrix X.
            //
        {
            int i;
            int j;
            int k;
            double[] x;

            if (trans != 'n' && trans != 'N' &&
                trans != 't' && trans != 'T' &&
                trans != 'c' && trans != 'C')
            {
                return null;
            }

            if (n < 0)
            {
                return null;
            }

            if (nrhs < 0)
            {
                return null;
            }

            if (n == 0 || nrhs == 0)
            {
                return null;
            }

            x = new double[n * nrhs];

            for (j = 0; j < nrhs; j++)
            {
                for (i = 0; i < n; i++)
                {
                    x[i + j * n] = b[i + j * n];
                }
            }

            if (trans == 'n' || trans == 'N')
            {
                //
                //  Solve L * x = b, overwriting b with x.
                //
                for (k = 0; k < nrhs; k++)
                {
                    for (j = 1; j <= n - 1; j++)
                    {
                        for (i = j + 1; i <= n; i++)
                        {
                            x[i - 1 + k * n] = x[i - 1 + k * n] - a[i - 1 + (j - 1) * n] * x[j - 1 + k * n];
                        }
                    }
                }

                //
                //  Solve U * x = b, overwriting b with x.
                //
                for (k = 0; k < nrhs; k++)
                {
                    for (j = n; 1 <= j; j--)
                    {
                        x[j - 1 + k * n] = x[j - 1 + k * n] / a[j - 1 + (j - 1) * n];
                        for (i = 1; i <= j - 1; i++)
                        {
                            x[i - 1 + k * n] = x[i - 1 + k * n] - a[i - 1 + (j - 1) * n] * x[j - 1 + k * n];
                        }
                    }
                }
            }
            else
                //
                //  Solve U' * x = b, overwriting b with x.
                //
            {
                for (k = 0; k < nrhs; k++)
                {
                    for (j = 1; j <= n; j++)
                    {
                        x[j - 1 + k * n] = x[j - 1 + k * n] / a[j - 1 + (j - 1) * n];
                        for (i = j + 1; i <= n; i++)
                        {
                            x[i - 1 + k * n] = x[i - 1 + k * n] - a[j - 1 + (i - 1) * n] * x[j - 1 + k * n];
                        }
                    }
                }

                //
                //  Solve L' * x = b, overwriting b with x.
                //
                for (k = 0; k < nrhs; k++)
                {
                    for (j = n; 2 <= j; j--)
                    {
                        for (i = 1; i <= j - 1; i++)
                        {
                            x[i - 1 + k * n] = x[i - 1 + k * n] - a[j - 1 + (i - 1) * n] * x[j - 1 + k * n];
                        }
                    }
                }
            }

            return x;
        }

    }
}