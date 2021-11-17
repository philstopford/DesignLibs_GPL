using System;
using Burkardt.Uniform;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double r8lt_det(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8LT_DET computes the determinant of an R8LT matrix.
        //
        //  Discussion:
        //
        //    The R8LT storage format is used for an M by N lower triangular 
        //    matrix A, and allocates storage even for the zero entries.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 September 2003
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
        //    Input, double A[N*N], the R8LT matrix.
        //
        //    Output, double R8LT_DET, the determinant of the matrix.
        //
    {
        double det;
        int i;

        det = 1.0;
        for (i = 0; i < n; i++)
        {
            det *= a[i + i * n];
        }

        return det;
    }

    public static double[] r8lt_indicator(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8LT_INDICATOR sets up an R8LT indicator matrix.
        //
        //  Discussion:
        //
        //    The R8LT storage format is used for an M by N lower triangular 
        //    matrix A, and allocates storage even for the zero entries.  
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 February 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns of the matrix.
        //    M and N must be positive.
        //
        //    Output, double R8LT_INDICATOR[M*N], the R8LT matrix.
        //
    {
        double[] a;
        int fac;
        int i;
        int j;

        a = new double[m * n];

        fac = (int) Math.Pow(10, (int) Math.Log10(n) + 1);

        for (i = 1; i <= m; i++)
        {
            for (j = 1; j <= Math.Min(i, n); j++)
            {
                a[i - 1 + (j - 1) * m] = fac * i + j;
            }

            for (j = i + 1; j <= n; j++)
            {
                a[i - 1 + (j - 1) * m] = 0.0;
            }
        }

        return a;
    }

    public static double[] r8lt_inverse(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8LT_INVERSE computes the inverse of an R8LT matrix.
        //
        //  Discussion:
        //
        //    The R8LT storage format is used for an M by N lower triangular 
        //    matrix A, and allocates storage even for the zero entries.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 February 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Albert Nijenhuis, Herbert Wilf,
        //    Combinatorial Algorithms,
        //    Second edition,
        //    Academic Press, 1978,
        //    ISBN 0-12-519260-6
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[N*N], the R8LT matrix.
        //
        //    Output, double R8LT_INVERSE[N*N], the inverse of the matrix.
        //
    {
        double[] b;
        int i;
        int j;
        int k;
        double t;
        //
        //  Check.
        //
        for (i = 0; i < n; i++)
        {
            switch (a[i + i * n])
            {
                case 0.0:
                    Console.WriteLine("");
                    Console.WriteLine("R8LT_INVERSE - Fatal error!");
                    Console.WriteLine("  Zero diagonal element.");
                    return null;
            }
        }

        b = new double[n * n];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                b[i + j * n] = a[i + j * n];
            }
        }

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                if (i < j)
                {
                    b[i + j * n] = 0.0;
                }
                else if (i == j)
                {
                    b[i + j * n] = 1.0 / b[i + j * n];
                }
                else if (j < i)
                {
                    t = 0.0;
                    for (k = j; k <= i - 1; k++)
                    {
                        t -= b[i + k * n] * b[k + j * n];
                    }

                    b[i + j * n] = t / b[i + i * n];
                }
            }
        }

        return b;
    }

    public static double[] r8lt_mm(int n, double[] a, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8LT_MM computes C=A*B for R8LT matrices.
        //
        //  Discussion:
        //
        //    The R8LT storage format is used for an M by N lower triangular 
        //    matrix A, and allocates storage even for the zero entries.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 February 2004
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
        //    Input, double A[N*N], B[N*N], the R8LT factor matrices.
        //
        //    Output, double R8LT_MM[N*N], the R8LT product matrix.
        //
    {
        double[] c;
        int i;
        int j;
        int k;

        c = new double[n * n];

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                c[i + j * n] = 0.0;
                for (k = j; k <= i; k++)
                {
                    c[i + j * n] += a[i + k * n] * b[k + j * n];
                }
            }
        }

        return c;
    }

    public static double[] r8lt_mtm(int n, double[] a, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8LT_MTM computes C=A'*B for R8LT matrices.
        //
        //  Discussion:
        //
        //    The R8LT storage format is used for an M by N lower triangular 
        //    matrix A, and allocates storage even for the zero entries.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 August 2015
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
        //    Input, double A[N*N], B[N*N], the R8LT factor matrices.
        //
        //    Output, double R8LT_MTM[N*N], the R8LT product matrix.
        //
    {
        double[] c;
        int i;
        int j;
        int k;

        c = new double[n * n];

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                c[i + j * n] = 0.0;
                for (k = Math.Min(i, j); k < n; k++)
                {
                    c[i + j * n] += a[k + i * n] * b[k + j * n];
                }
            }
        }

        return c;
    }

    public static double[] r8lt_mtv(int m, int n, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8LT_MTV multiplies a vector times an R8LT matrix.
        //
        //  Discussion:
        //
        //    The R8LT storage format is used for an M by N lower triangular 
        //    matrix A, and allocates storage even for the zero entries.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 September 2003
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
        //    Input, double A[M*N], the R8LT matrix.
        //
        //    Input, double X[M], the vector to be multiplied by A.
        //
        //    Output, double R8LT_MTV[N], the product A * x.
        //
    {
        double[] b;
        int i;
        int j;

        b = new double[n];

        for (j = 0; j < n; j++)
        {
            b[j] = 0.0;
            for (i = j; i < m; i++)
            {
                b[j] += x[i] * a[i + j * m];
            }
        }

        return b;
    }

    public static double[] r8lt_mv(int m, int n, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8LT_MV multiplies an R8LT matrix times a vector.
        //
        //  Discussion:
        //
        //    The R8LT storage format is used for an M by N lower triangular 
        //    matrix A, and allocates storage even for the zero entries.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 September 2003
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
        //    Input, double A[M*N], the R8LT matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double R8LT_MV[M], the product A * x.
        //
    {
        double[] b;
        int i;
        int j;
        int jmax;

        b = new double[m];

        for (i = 0; i < m; i++)
        {
            b[i] = 0.0;
            jmax = Math.Min(i, n - 1);
            for (j = 0; j <= jmax; j++)
            {
                b[i] += a[i + j * m] * x[j];
            }
        }

        return b;
    }

    public static void r8lt_print(int m, int n, double[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8LT_PRINT prints an R8LT matrix.
        //
        //  Discussion:
        //
        //    The R8LT storage format is used for an M by N lower triangular 
        //    matrix A, and allocates storage even for the zero entries.
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
        //    Input, double A[M*N], the R8LT matrix.
        //
        //    Input, string TITLE, a title.
        //
    {
        r8lt_print_some(m, n, a, 1, 1, m, n, title);
    }

    public static void r8lt_print_some(int m, int n, double[] a, int ilo, int jlo, int ihi,
            int jhi, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8LT_PRINT_SOME prints some of an R8LT matrix.
        //
        //  Discussion:
        //
        //    The R8LT storage format is used for an M by N lower triangular 
        //    matrix A, and allocates storage even for the zero entries.
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
        //    Input, double A[M*N], the R8LT matrix.
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
            cout = "  Col: ";

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
            i2lo = Math.Max(i2lo, j2lo);

            i2hi = Math.Min(ihi, m);

            for (i = i2lo; i <= i2hi; i++)
            {
                //
                //  Print out (up to) 5 entries in row I, that lie in the current strip.
                //
                cout = i.ToString().PadLeft(4) + "  ";

                for (j = j2lo; j <= j2hi; j++)
                {
                    if (i < j)
                    {
                        cout += "              ";
                    }
                    else
                    {
                        cout += a[i - 1 + (j - 1) * m].ToString().PadLeft(12) + "  ";
                    }
                }

                Console.WriteLine(cout);
            }
        }

        Console.WriteLine("");
    }

    public static double[] r8lt_random(int m, int n, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8LT_RANDOM randomizes an R8LT matrix.
        //
        //  Discussion:
        //
        //    The R8LT storage format is used for an M by N lower triangular 
        //    matrix A, and allocates storage even for the zero entries.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 August 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns of the matrix.
        //    M and N must be positive.
        //
        //    Input/output, int SEED, a seed for the random number generator.
        //
        //    Output, double R8LT_RANDOM[M*N], the R8LT matrix.
        //
    {
        double[] a;
        int i;
        int j;

        a = new double[m * n];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < j; i++)
            {
                a[i + j * m] = 0.0;
            }

            for (i = j; i < m; i++)
            {
                a[i + j * m] = UniformRNG.r8_uniform_01(ref seed);
            }
        }

        return a;
    }

    public static double[] r8lt_sl(int n, double[] a, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8LT_SL solves A*x=b, where A is an R8LT matrix.
        //
        //  Discussion:
        //
        //    The R8LT storage format is used for an M by N lower triangular 
        //    matrix A, and allocates storage even for the zero entries.
        //
        //    No factorization of the lower triangular matrix is required.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 August 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[N*N], the R8LT matrix.
        //
        //    Input, double B[N], the right hand side.
        //
        //    Output, double R8LT_SL[N], the solution vector.
        //
    {
        int i;
        int j;
        double[] x;

        x = new double[n];

        for (i = 0; i < n; i++)
        {
            x[i] = b[i];
        }

        for (j = 0; j < n; j++)
        {
            x[j] /= a[j + j * n];
            for (i = j + 1; i < n; i++)
            {
                x[i] -= a[i + j * n] * x[j];
            }
        }

        return x;
    }

    public static double[] r8lt_slt(int n, double[] a, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8LT_SLT solves A'*x=b, where A is an R8LT matrix.
        //
        //  Discussion:
        //
        //    The R8LT storage format is used for an M by N lower triangular 
        //    matrix A, and allocates storage even for the zero entries.
        //
        //    No factorization of the lower triangular matrix is required.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 August 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[N*N], the R8LT matrix.
        //
        //    Input, double B[N], the right hand side.
        //
        //    Output, double R8LT_SLT[N], the solution vector.
        //
    {
        int i;
        int j;
        double[] x;

        x = new double[n];

        for (i = 0; i < n; i++)
        {
            x[i] = b[i];
        }

        for (j = n - 1; 0 <= j; j--)
        {
            x[j] /= a[j + j * n];
            for (i = 0; i < j; i++)
            {
                x[i] -= a[j + i * n] * x[j];
            }
        }

        return x;
    }

    public static double[] r8lt_to_r8ge(int m, int n, double[] a_lt)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8LT_TO_R8GE copies an R8LT matrix to an R8GE matrix.
        //
        //  Discussion:
        //
        //    The R8LT storage format is used for an M by N lower triangular matrix,
        //    and allocates space even for the zero entries.
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
        //    21 August 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Input, double A_LT[M,N], the R8LT matrix.
        //
        //    Output, double R8LT_TO_R8GE[M,N], the R8GE matrix.
        //
    {
        double[] a_ge;
        int i;
        int j;

        a_ge = new double[m * n];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                if (j <= i)
                {
                    a_ge[i + j * m] = a_lt[i + j * m];
                }
                else
                {
                    a_ge[i + j * m] = 0.0;
                }
            }
        }

        return a_ge;
    }

    public static double[] r8lt_zeros(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8LT_ZEROS zeros an R8LT matrix.
        //
        //  Discussion:
        //
        //    The R8LT storage format is used for an M by N lower triangular 
        //    matrix A, and allocates storage even for the zero entries.
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
        //    Input, int M, the number of rows of the matrix.
        //    M must be positive.
        //
        //    Input, int N, the number of columns of the matrix.
        //    N must be positive.
        //
        //    Output, double R8LT_ZERO[M*N], the R8LT matrix.
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
                a[i + j * m] = 0.0;
            }
        }

        return a;
    }
}