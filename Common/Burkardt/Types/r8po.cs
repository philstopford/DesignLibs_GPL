using System;
using Burkardt.BLAS;
using Burkardt.Uniform;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static int r8po_fa(ref double[] a, int lda, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PO_FA factors a real symmetric positive definite matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 May 2005
        //
        //  Author:
        //
        //    FORTRAN77 original version by Dongarra, Moler, Bunch, Stewart.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Dongarra, Moler, Bunch and Stewart,
        //    LINPACK User's Guide,
        //    SIAM, (Society for Industrial and Applied Mathematics),
        //    3600 University City Science Center,
        //    Philadelphia, PA, 19104-2688.
        //    ISBN 0-89871-172-X
        //
        //  Parameters:
        //
        //    Input/output, double A[LDA*N].  On input, the symmetric matrix
        //    to be  factored.  Only the diagonal and upper triangle are used.
        //    On output, an upper triangular matrix R so that A = R'*R
        //    where R' is the transpose.  The strict lower triangle is unaltered.
        //    If INFO /= 0, the factorization is not complete.
        //
        //    Input, int LDA, the leading dimension of the array A.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Output, int R8PO_FA, error flag.
        //    0, for normal return.
        //    K, signals an error condition.  The leading minor of order K is not
        //    positive definite.
        //
    {
        int info;
        int j;
        int k;
        double s;
        double t;

        for (j = 1; j <= n; j++)
        {
            s = 0.0;

            for (k = 1; k <= j - 1; k++)
            {
                t = a[k - 1 + (j - 1) * lda] -
                    BLAS1D.ddot(k - 1, a, 1, a, 1, +0 + (k - 1) * lda, +0 + (j - 1) * lda);
                t /= a[k - 1 + (k - 1) * lda];
                a[k - 1 + (j - 1) * lda] = t;
                s += t * t;
            }

            s = a[j - 1 + (j - 1) * lda] - s;

            switch (s)
            {
                case <= 0.0:
                    info = j;
                    return info;
                default:
                    a[j - 1 + (j - 1) * lda] = Math.Sqrt(s);
                    break;
            }
        }

        info = 0;

        return info;
    }

    public static double r8po_det(int n, double[] a_lu)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PO_DET computes the determinant of a matrix factored by R8PO_FA.
        //
        //  Discussion:
        //
        //    The R8PO storage format is appropriate for a symmetric positive definite 
        //    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
        //    upper triangular matrix, so it will be in R8GE storage format.)
        //
        //    Only the diagonal and upper triangle of the square array are used.
        //    This same storage format is used when the matrix is factored by
        //    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
        //    is set to zero.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A_LU[N*N], the LU factors from R8PO_FA.
        //
        //    Output, double R8PO_DET, the determinant of A.
        //
    {
        double det;
        int i;

        det = 1.0;

        for (i = 0; i < n; i++)
        {
            det *= Math.Pow(a_lu[i + i * n], 2);
        }

        return det;
    }

    public static double[] r8po_dif2(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PO_DIF2 returns the second difference matrix in R8PO format.
        //
        //  Discussion:
        //
        //    The R8PO storage format is appropriate for a symmetric positive definite 
        //    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
        //    upper triangular matrix, so it will be in R8GE storage format.)
        //
        //    Only the diagonal and upper triangle of the square array are used.
        //    This same storage format is used when the matrix is factored by
        //    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
        //    is set to zero.
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
        //    Input, int N, the number of rows and columns of the matrix.
        //    N must be positive.
        //
        //    Output, double R8PO_DIF2[N*N], the matrix.
        //
    {
        double[] a;
        int i;
        int j;

        a = new double[n * n];

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                if (j == i)
                {
                    a[i + j * n] = 2.0;
                }
                else if (j == i + 1)
                {
                    a[i + j * n] = -1.0;
                }
                else
                {
                    a[i + j * n] = 0.0;
                }
            }
        }

        return a;
    }

    public static double[] r8po_fa(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PO_FA factors an R8PO matrix.
        //
        //  Discussion:
        //
        //    The R8PO storage format is appropriate for a symmetric positive definite 
        //    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
        //    upper triangular matrix, so it will be in R8GE storage format.)
        //
        //    Only the diagonal and upper triangle of the square array are used.
        //    This same storage format is used when the matrix is factored by
        //    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
        //    is set to zero.
        //
        //    The positive definite symmetric matrix A has a Cholesky factorization
        //    of the form:
        //
        //      A = R' * R
        //
        //    where R is an upper triangular matrix with positive elements on
        //    its diagonal.  
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
        //
        //    Input, double A[N*N], the matrix in R8PO storage.
        //
        //    Output, double R8PO_FA[N*N], the Cholesky factor in SGE
        //    storage, or NULL if there was an error.
        //
    {
        double s;

        double[] b = new double[n * n];

        for (int j = 0; j < n; j++)
        {
            for (int i = 0; i < n; i++)
            {
                b[i + j * n] = a[i + j * n];
            }
        }

        for (int j = 0; j < n; j++)
        {
            for (int k = 0; k <= j - 1; k++)
            {
                for (int i = 0; i <= k - 1; i++)
                {
                    b[k + j * n] -= b[i + k * n] * b[i + j * n];
                }

                b[k + j * n] /= b[k + k * n];
            }

            s = b[j + j * n];
            for (int i = 0; i <= j - 1; i++)
            {
                s -= b[i + j * n] * b[i + j * n];
            }

            switch (s)
            {
                case <= 0.0:
                    return null;
                default:
                    b[j + j * n] = Math.Sqrt(s);
                    break;
            }
        }

        //
        //  Since the Cholesky factor is in R8GE format, zero out the lower triangle.
        //
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < i; j++)
            {
                b[i + j * n] = 0.0;
            }
        }

        return b;
    }

    public static double[] r8po_sl(int n, double[] a_lu, double[] b, int bIndex = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PO_SL solves a linear system that has been factored by R8PO_FA.
        //
        //  Discussion:
        //
        //    The R8PO storage format is appropriate for a symmetric positive definite 
        //    matrix and its inverse.  (The Cholesky factor of a R8PO matrix is an
        //    upper triangular matrix, so it will be in R8GE storage format.)
        //
        //    Only the diagonal and upper triangle of the square array are used.
        //    This same storage format is used when the matrix is factored by
        //    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
        //    is set to zero.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 February 2004
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
        //
        //    Input, double A_LU[N*N], the Cholesky factor from R8PO_FA.
        //
        //    Input, double B[N], the right hand side.
        //
        //    Output, double R8PO_SL[N], the solution vector.
        //
    {
        int i;
        int k;
        double[] x;

        x = new double[n];

        for (k = 0; k < n; k++)
        {
            x[k] = b[bIndex + k];
        }

        //
        //  Solve R' * y = b.
        //
        for (k = 0; k < n; k++)
        {
            for (i = 0; i < k; i++)
            {
                x[k] -= x[i] * a_lu[i + k * n];
            }

            x[k] /= a_lu[k + k * n];
        }

        //
        //  Solve R * x = y.
        //
        for (k = n - 1; 0 <= k; k--)
        {
            x[k] /= a_lu[k + k * n];
            for (i = 0; i < k; i++)
            {
                x[i] -= a_lu[i + k * n] * x[k];
            }
        }

        return x;
    }

    public static double[] r8po_indicator(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PO_INDICATOR sets up an R8PO indicator matrix.
        //
        //  Discussion:
        //
        //    The R8PO storage format is appropriate for a symmetric positive definite 
        //    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
        //    upper triangular matrix, so it will be in R8GE storage format.)
        //
        //    Only the diagonal and upper triangle of the square array are used.
        //    This same storage format is used when the matrix is factored by
        //    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
        //    is set to zero.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 January 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of rows and columns of the matrix.
        //    N must be positive.
        //
        //    Output, double R8PO_INDICATOR[N*N], the R8PO matrix.
        //
    {
        double[] a;
        int fac;
        int i;
        int j;

        a = new double[n * n];

        fac = (int) Math.Pow(10, (int) Math.Log10(n) + 1);

        for (i = 1; i <= n; i++)
        {
            for (j = 1; j <= i - 1; j++)
            {
                a[i - 1 + (j - 1) * n] = 0.0;
            }

            for (j = i; j <= n; j++)
            {
                a[i - 1 + (j - 1) * n] = fac * i + j;
            }
        }

        return a;
    }

    public static double[] r8po_inverse(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PO_INVERSE computes the inverse of a matrix factored by R8PO_FA.
        //
        //  Discussion:
        //
        //    The R8PO storage format is appropriate for a symmetric positive definite 
        //    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
        //    upper triangular matrix, so it will be in R8GE storage format.)
        //
        //    Only the diagonal and upper triangle of the square array are used.
        //    This same storage format is used when the matrix is factored by
        //    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
        //    is set to zero.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 February 2004
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
        //
        //    Input, double A[N*N], the Cholesky factor, in R8GE storage, returned by R8PO_FA.
        //
        //    Output, double R8PO_INVERSE[N*N], the inverse, in R8PO storage.
        //
    {
        double[] b;
        int i;
        int j;
        int k;
        double t;

        b = new double[n * n];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                b[i + j * n] = a[i + j * n];
            }
        }

        //
        //  Compute Inverse ( R ).
        //
        for (k = 0; k < n; k++)
        {
            b[k + k * n] = 1.0 / b[k + k * n];
            for (i = 0; i < k; i++)
            {
                b[i + k * n] = -b[i + k * n] * b[k + k * n];
            }

            for (j = k + 1; j < n; j++)
            {
                t = b[k + j * n];
                b[k + j * n] = 0.0;
                for (i = 0; i <= k; i++)
                {
                    b[i + j * n] += t * b[i + k * n];
                }
            }
        }

        //
        //  Compute Inverse ( R ) * Transpose ( Inverse ( R ) ).
        //
        for (j = 0; j < n; j++)
        {
            for (k = 0; k < j; k++)
            {
                t = b[k + j * n];
                for (i = 0; i <= k; i++)
                {
                    b[i + k * n] += t * b[i + j * n];
                }
            }

            t = b[j + j * n];
            for (i = 0; i <= j; i++)
            {
                b[i + j * n] *= t;
            }
        }

        return b;
    }

    public static double[] r8po_ml(int n, double[] a_lu, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PO_ML computes A * x = b after A has been factored by R8PO_FA.
        //
        //  Discussion:
        //
        //    The R8PO storage format is appropriate for a symmetric positive definite 
        //    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
        //    upper triangular matrix, so it will be in R8GE storage format.)
        //
        //    Only the diagonal and upper triangle of the square array are used.
        //    This same storage format is used when the matrix is factored by
        //    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
        //    is set to zero.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 October 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A_LU[N*N], the Cholesky factor from R8PO_FA.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double R8PO_ML[N], the product A * x.
        //
    {
        double[] b;
        int i;
        int j;

        b = new double[n];
        //
        //  Compute R * x = y.
        //
        for (i = 0; i < n; i++)
        {
            b[i] = a_lu[i + i * n] * x[i];
            for (j = i + 1; j < n; j++)
            {
                b[i] += a_lu[i + j * n] * x[j];
            }
        }

        //
        //  Compute R' * y = b.
        //
        for (j = n - 1; 0 <= j; j--)
        {
            b[j] = a_lu[j + j * n] * b[j];
            for (i = 0; i < j; i++)
            {
                b[j] += b[i] * a_lu[i + j * n];
            }
        }

        return b;
    }

    public static double[] r8po_mm(int n, double[] a, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PO_MM multiplies two R8PO matrices.
        //
        //  Discussion:
        //
        //    The R8PO storage format is appropriate for a symmetric positive definite
        //    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
        //    upper triangular matrix, so it will be in R8GE storage format.)
        //
        //    Only the diagonal and upper triangle of the square array are used.
        //    This same storage format is used when the matrix is factored by
        //    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
        //    is set to zero.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 December 2003
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
        //    Input, double A[N*N], B[N*N], the R8PO factor matrices.
        //
        //    Output, double R8PO_MM[N*N], the R8PO product matrix.
        //
    {
        double aik;
        double bkj;
        double[] c;
        int i;
        int j;
        int k;

        c = new double[n * n];

        for (i = 1; i <= n; i++)
        {
            for (j = 1; j <= n; j++)
            {
                c[i - 1 + (j - 1) * n] = 0.0;
            }
        }

        for (i = 1; i <= n; i++)
        {
            for (j = i; j <= n; j++)
            {
                for (k = 1; k <= n; k++)
                {
                    if (i <= k)
                    {
                        aik = a[i - 1 + (k - 1) * n];
                    }
                    else
                    {
                        aik = a[k - 1 + (i - 1) * n];
                    }

                    if (k <= j)
                    {
                        bkj = b[k - 1 + (j - 1) * n];
                    }
                    else
                    {
                        bkj = b[j - 1 + (k - 1) * n];
                    }

                    c[i - 1 + (j - 1) * n] += aik * bkj;

                }
            }

        }

        return c;
    }

    public static double[] r8po_mv(int n, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PO_MV multiplies an R8PO matrix times a vector.
        //
        //  Discussion:
        //
        //    The R8PO storage format is appropriate for a symmetric positive definite 
        //    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
        //    upper triangular matrix, so it will be in R8GE storage format.)
        //
        //    Only the diagonal and upper triangle of the square array are used.
        //    This same storage format is used when the matrix is factored by
        //    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
        //    is set to zero.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 October 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[N*N], the R8PO matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double R8PO_MV(N), the product A * x.
        //
    {
        double[] b;
        int i;
        int j;

        b = new double[n];

        for (i = 0; i < n; i++)
        {
            b[i] = 0.0;
            for (j = 0; j < i; j++)
            {
                b[i] += a[j + i * n] * x[j];
            }

            for (j = i; j < n; j++)
            {
                b[i] += a[i + j * n] * x[j];
            }
        }

        return b;
    }

    public static void r8po_print(int n, double[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PO_PRINT prints an R8PO matrix.
        //
        //  Discussion:
        //
        //    The R8PO storage format is appropriate for a symmetric positive definite 
        //    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
        //    upper triangular matrix, so it will be in R8GE storage format.)
        //
        //    Only the diagonal and upper triangle of the square array are used.
        //    This same storage format is used when the matrix is factored by
        //    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
        //    is set to zero.
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
        //
        //    Input, double A[M*N], the R8PO matrix.
        //
        //    Input, string TITLE, a title.
        //
    {
        r8po_print_some(n, a, 1, 1, n, n, title);
    }

    public static void r8po_print_some(int n, double[] a, int ilo, int jlo, int ihi,
            int jhi, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PO_PRINT_SOME prints some of an R8PO matrix.
        //
        //  Discussion:
        //
        //    The R8PO storage format is appropriate for a symmetric positive definite 
        //    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
        //    upper triangular matrix, so it will be in R8GE storage format.)
        //
        //    Only the diagonal and upper triangle of the square array are used.
        //    This same storage format is used when the matrix is factored by
        //    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
        //    is set to zero.
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
        //    Input, double A[M*N], the R8PO matrix.
        //
        //    Input, int ILO, JLO, IHI, JHI, designate the first row and
        //    column, and the last row and column to be printed.
        //
        //    Input, string TITLE, a title.
        //
    {
        int INCX = 5;

        int i;
        int i2hi;
        int i2lo;
        int j;
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

            Console.WriteLine("");
            //
            //  For each column J in the current range...
            //
            //  Write the header.
            //
            cout = "  Col:    ";
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
            i2lo = Math.Max(ilo, 1);
            i2hi = Math.Min(ihi, n);

            for (i = i2lo; i <= i2hi; i++)
            {
                //
                //  Print out (up to) 5 entries in row I, that lie in the current strip.
                //
                cout = i.ToString().PadLeft(5) + "  ";
                for (j = j2lo; j <= j2hi; j++)
                {
                    if (i <= j)
                    {
                        cout += a[i - 1 + (j - 1) * n].ToString().PadLeft(12) + "  ";
                    }
                    else
                    {
                        cout += a[j - 1 + (i - 1) * n].ToString().PadLeft(12) + "  ";
                    }
                }

                Console.WriteLine(cout);
            }
        }

        Console.WriteLine("");
    }

    public static double[] r8po_random(int n, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PO_RANDOM randomizes an R8PO matrix.
        //
        //  Discussion:
        //
        //    The R8PO storage format is appropriate for a symmetric positive definite 
        //    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
        //    upper triangular matrix, so it will be in R8GE storage format.)
        //
        //    Only the diagonal and upper triangle of the square array are used.
        //    This same storage format is used when the matrix is factored by
        //    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
        //    is set to zero.
        //
        //    The matrix computed here is not simply a set of random numbers in
        //    the nonzero slots of the R8PO array.  It is also a positive definite
        //    matrix.  It is computed by setting a "random" upper triangular
        //    Cholesky factor R, and then computing A = R'*R.
        //    The randomness is limited by the fact that all the entries of
        //    R will be between 0 and 1.  A truly random R is only required
        //    to have positive entries on the diagonal.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 December 2003
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
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double R8PO_RANDOM[N*N], the R8PO matrix.
        //
    {
        double[] a;
        int i;
        int j;
        int k;

        a = new double[n * n];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                a[i + j * n] = 0.0;
            }
        }

        for (i = n; 1 <= i; i--)
        {
            //
            //  Set row I of R.
            //
            for (j = i; j <= n; j++)
            {
                a[i - 1 + (j - 1) * n] = UniformRNG.r8_uniform_01(ref seed);
            }

            //
            //  Consider element J of row I, last to first.
            //
            for (j = n; i <= j; j--)
            {
                //
                //  Add multiples of row I to lower elements of column J.
                //
                for (k = i + 1; k <= j; k++)
                {
                    a[k - 1 + (j - 1) * n] += a[i - 1 + (k - 1) * n] * a[i - 1 + (j - 1) * n];
                }

                //
                //  Reset element J.
                //
                a[i - 1 + (j - 1) * n] = a[i - 1 + (i - 1) * n] * a[i - 1 + (j - 1) * n];
            }
        }

        return a;
    }

    public static void r8po_sl(double[] a, int lda, int n, ref double[] b, int aIndex = 0, int bIndex = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PO_SL solves a linear system factored by R8PO_FA.
        //
        //  Discussion:
        //
        //    A division by zero will occur if the input factor contains
        //    a zero on the diagonal.  Technically this indicates
        //    singularity but it is usually caused by improper subroutine
        //    arguments.  It will not occur if the subroutines are called
        //    correctly and INFO == 0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 May 2005
        //
        //  Author:
        //
        //    FORTRAN77 original version by Dongarra, Moler, Bunch, Stewart.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Dongarra, Moler, Bunch and Stewart,
        //    LINPACK User's Guide,
        //    SIAM, (Society for Industrial and Applied Mathematics),
        //    3600 University City Science Center,
        //    Philadelphia, PA, 19104-2688.
        //    ISBN 0-89871-172-X
        //
        //  Parameters:
        //
        //    Input, double A[LDA*N], the output from R8PO_FA.
        //
        //    Input, int LDA, the leading dimension of the array A.
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input/output, double B[N].  On input, the right hand side.
        //    On output, the solution.
        //
    {
        int k;
        double t;
        //
        //  Solve R' * Y = B.
        //
        for (k = 1; k <= n; k++)
        {
            t = BLAS1D.ddot(k - 1, a, 1, b, 1, aIndex + 0 + (k - 1) * lda, bIndex);
            b[bIndex + k - 1] = (b[bIndex + k - 1] - t) / a[aIndex + k - 1 + (k - 1) * lda];
        }

        //
        //  Solve R * X = Y.
        //
        for (k = n; 1 <= k; k--)
        {
            b[bIndex + k - 1] /= a[aIndex + k - 1 + (k - 1) * lda];
            t = -b[bIndex + k - 1];
            BLAS1D.daxpy(k - 1, t, a, 1, ref b, 1, aIndex + 0 + (k - 1) * lda, bIndex);
        }
    }
        
    public static double[] r8po_to_r8ge ( int n, double[] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PO_TO_R8GE copies an R8PO matrix to an R8GE matrix.
        //
        //  Discussion:
        //
        //    The R8PO storage format is appropriate for a symmetric positive definite 
        //    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
        //    upper triangular matrix, so it will be in R8GE storage format.)
        //
        //    Only the diagonal and upper triangle of the square array are used.
        //    This same storage format is used when the matrix is factored by
        //    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
        //    is set to zero.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 October 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[N*N], the R8PO matrix.
        //
        //    Output, double R8PO_TO_R8GE[N*N], the R8GE matrix.
        //
    {
        double[] b;
        int i;
        int j;

        b = new double[n*n];

        for ( i = 0; i < n; i++ )
        {
            for ( j = 0; j < n; j++ )
            {
                if ( i <= j )
                {
                    b[i+j*n] = a[i+j*n];
                }
                else
                {
                    b[i+j*n] = a[j+i*n];
                }
            }
        }

        return b;
    }

    public static double[] r8po_zeros ( int n )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8PO_ZEROS zeros an R8PO matrix.
        //
        //  Discussion:
        //
        //    The R8PO storage format is appropriate for a symmetric positive definite 
        //    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
        //    upper triangular matrix, so it will be in R8GE storage format.)
        //
        //    Only the diagonal and upper triangle of the square array are used.
        //    This same storage format is used when the matrix is factored by
        //    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
        //    is set to zero.
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
        //    Input, int N, the number of rows and columns of the matrix.
        //    N must be positive.
        //
        //    Output, double R8PO_ZERO[N*N], the R8PO matrix.
        //
    {
        double[] a;
        int i;
        int j;

        a = new double[n*n];

        for ( j = 0; j < n; j++ )
        {
            for ( i = 0; i < n; i++ )
            {
                a[i+j*n] = 0.0;
            }
        }

        return a;
    }
}