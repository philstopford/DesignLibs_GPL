using System;
using Burkardt.Uniform;

namespace Burkardt.Types;

public static partial class typeMethods
{

    public static double r8ut_det(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8UT_DET computes the determinant of an R8UT matrix.
        //
        //  Discussion:
        //
        //    The R8UT storage format is used for an M by N upper triangular matrix,
        //    and allocates space even for the zero entries.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 September 2003
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
        //    Input, double A[N*N], the R8UT matrix.
        //
        //    Output, double R8UT_DET, the determinant of the matrix.
        //
    {
        int i;

        double det = 1.0;

        for (i = 0; i < n; i++)
        {
            det *= a[i + i * n];
        }

        return det;
    }

    public static double[] r8ut_indicator(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8UT_INDICATOR sets up an R8UT indicator matrix.
        //
        //  Discussion:
        //
        //    The R8UT storage format is used for an M by N upper triangular matrix,
        //    and allocates space even for the zero entries.
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
        //    Output, double R8UT_INDICATOR[M*N], the R8UT matrix.
        //
    {
        int i;

        double[] a = new double[m * n];

        int fac = (int) Math.Pow(10, (int) Math.Log10(n) + 1);

        for (i = 1; i <= m; i++)
        {
            int j;
            for (j = 1; j <= Math.Min(i - 1, n); j++)
            {
                a[i - 1 + (j - 1) * m] = 0.0;
            }

            for (j = i; j <= n; j++)
            {
                a[i - 1 + (j - 1) * m] = fac * i + j;
            }
        }

        return a;
    }

    public static double[] r8ut_inverse(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8UT_INVERSE computes the inverse of a R8UT matrix.
        //
        //  Discussion:
        //
        //    The R8UT storage format is used for an M by N upper triangular matrix,
        //    and allocates space even for the zero entries.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Albert Nijenhuis, Herbert Wilf,
        //    Combinatorial Algorithms,
        //    Academic Press, 1978, second edition,
        //    ISBN 0-12-519260-6
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[N*N], the R8UT matrix.
        //
        //    Output, double R8UT_INVERSE[N*N], the inverse of the upper
        //    triangular matrix.
        //
    {
        int i;
        int j;
        //
        //  Check.
        //
        for (i = 0; i < n; i++)
        {
            switch (a[i + i * n])
            {
                case 0.0:
                    Console.WriteLine("");
                    Console.WriteLine("R8UT_INVERSE - Fatal error!");
                    Console.WriteLine("  Zero diagonal element.");
                    return null;
            }
        }

        double[] b = new double[n * n];

        for (j = n - 1; 0 <= j; j--)
        {
            for (i = n - 1; 0 <= i; i--)
            {
                if (j < i)
                {
                    b[i + j * n] = 0.0;
                }
                else if (i == j)
                {
                    b[i + j * n] = 1.0 / a[i + j * n];
                }
                else if (i < j)
                {
                    b[i + j * n] = 0.0;

                    int k;
                    for (k = i + 1; k <= j; k++)
                    {
                        b[i + j * n] -= a[i + k * n] * b[k + j * n];
                    }

                    b[i + j * n] /= a[i + i * n];
                }
            }
        }

        return b;
    }

    public static double[] r8ut_mm(int n, double[] a, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8UT_MM computes C = A * B, where A and B are R8UT matrices.
        //
        //  Discussion:
        //
        //    The R8UT storage format is used for an M by N upper triangular matrix,
        //    and allocates space even for the zero entries.
        //
        //    The product C will also be an upper trangular matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 August 2015
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
        //    Input, double A[N*N], B[N*N], the R8UT factor matrices.
        //
        //    Output, double R8UT_MM[N*N], the R8UT product matrix.
        //
    {
        int i;

        double[] c = new double[n * n];

        for (i = 0; i < n; i++)
        {
            int j;
            for (j = 0; j < i; j++)
            {
                c[i + j * n] = 0.0;
            }

            for (j = i; j < n; j++)
            {
                c[i + j * n] = 0.0;
                int k;
                for (k = i; k <= j; k++)
                {
                    c[i + j * n] += a[i + k * n] * b[k + j * n];
                }
            }
        }

        return c;
    }

    public static double[] r8ut_mtm(int n, double[] a, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8UT_MTM computes C = A' * B, where A and B are R8UT matrices.
        //
        //  Discussion:
        //
        //    The R8UT storage format is used for an M by N upper triangular matrix,
        //    and allocates space even for the zero entries.
        //
        //    The product C will NOT be an R8UT matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 August 2015
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
        //    Input, double A[N*N], B[N*N], the factors.
        //
        //    Output, double R8UT_MTM[N*N], the product.
        //
    {
        int i;

        double[] c = new double[n * n];

        for (i = 0; i < n; i++)
        {
            int j;
            for (j = 0; j < n; j++)
            {
                int k_hi = Math.Min(i, j);
                c[i + j * n] = 0.0;
                int k;
                for (k = 0; k <= k_hi; k++)
                {
                    c[i + j * n] += a[k + i * n] * b[k + j * n];
                }
            }
        }

        return c;
    }

    public static double[] r8ut_mtv(int m, int n, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8UT_MTV multiplies a vector times an R8UT matrix.
        //
        //  Discussion:
        //
        //    The R8UT storage format is used for an M by N upper triangular matrix,
        //    and allocates space even for the zero entries.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 September 2003
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
        //    Input, double A[M*N], the R8UT matrix.
        //
        //    Input, double X[M], the vector to be multiplied by A.
        //
        //    Output, double R8UT_MTV[N], the product A' * x.
        //
    {
        int j;

        double[] b = new double[n];

        for (j = 0; j < n; j++)
        {
            b[j] = 0.0;
            int i;
            for (i = 0; i <= j && i < m; i++)
            {
                b[j] += x[i] * a[i + j * m];
            }
        }

        return b;
    }

    public static double[] r8ut_mv(int m, int n, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8UT_MV multiplies an R8UT matrix times a vector.
        //
        //  Discussion:
        //
        //    The R8UT storage format is used for an M by N upper triangular matrix,
        //    and allocates space even for the zero entries.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 September 2003
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
        //    Input, double A[M*N], the R8UT matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double R8UT_MV[M], the product A * x.
        //
    {
        int i;

        double[] b = new double[m];

        for (i = 0; i < m; i++)
        {
            b[i] = 0.0;
            int j;
            for (j = i; j < n; j++)
            {
                b[i] += a[i + j * m] * x[j];
            }
        }

        return b;
    }

    public static void r8ut_print(int m, int n, double[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8UT_PRINT prints an R8UT matrix.
        //
        //  Discussion:
        //
        //    The R8UT storage format is used for an M by N upper triangular matrix,
        //    and allocates space even for the zero entries.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 August 2015
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
        //    Input, double A[M*N], the R8UT matrix.
        //
        //    Input, string TITLE, a title.
        //
    {
        r8ut_print_some(m, n, a, 0, 0, m - 1, n - 1, title);
    }

    public static void r8ut_print_some(int m, int n, double[] a, int ilo, int jlo, int ihi,
            int jhi, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8UT_PRINT_SOME prints some of an R8UT matrix.
        //
        //  Discussion:
        //
        //    The R8UT storage format is used for an M by N upper triangular matrix,
        //    and allocates space even for the zero entries.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 August 2015
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
        //    Input, double A[M*N], the R8UT matrix.
        //
        //    Input, int ILO, JLO, IHI, JHI, designate the first row and
        //    column, and the last row and column to be printed.
        //    0 <= ILO <= IHI < M.
        //    0 <= JLO <= JHI < N.
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

            int i2hi = Math.Min(ihi, m - 1);

            int i;
            for (i = i2lo; i <= i2hi; i++)
            {
                //
                //  Print out (up to) 5 entries in row I, that lie in the current strip.
                //
                cout = i.ToString().PadLeft(4) + "  ";

                for (j = j2lo; j <= j2hi; j++)
                {
                    if (j < i)
                    {
                        cout += "              ";
                    }
                    else
                    {
                        cout += a[i + j * m].ToString().PadLeft(12) + "  ";
                    }
                }

                Console.WriteLine(cout);
            }
        }

    }

    public static double[] r8ut_random(int m, int n, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8UT_RANDOM randomizes an R8UT matrix.
        //
        //  Discussion:
        //
        //    The R8UT storage format is used for an M by N upper triangular matrix,
        //    and allocates space even for the zero entries.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 August 2015
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
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double R8UT_RANDOM[M*N], the R8UT matrix.
        //
    {
        int j;

        double[] a = new double[m * n];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i <= Math.Min(j, m - 1); i++)
            {
                a[i + j * m] = UniformRNG.r8_uniform_01(ref seed);
            }

            for (i = j + 1; i < m; i++)
            {
                a[i + j * m] = 0.0;
            }
        }

        return a;
    }

    public static double[] r8ut_sl(int n, double[] a, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8UT_SL solves a linear system A*x=b with R8UT matrix.
        //
        //  Discussion:
        //
        //    The R8UT storage format is used for an M by N upper triangular matrix,
        //    and allocates space even for the zero entries.
        //
        //    No factorization of the upper triangular matrix is required.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 August 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[N*N], the R8UT matrix.
        //
        //    Input, double B[N], the right hand side.
        //
        //    Output, double R8UT_SL[N], the solution vector.
        //
    {
        int i;
        int j;

        double[] x = new double[n];

        for (i = 0; i < n; i++)
        {
            x[i] = b[i];
        }

        for (j = n - 1; 0 <= j; j--)
        {
            x[j] /= a[j + j * n];
            for (i = 0; i < j; i++)
            {
                x[i] -= a[i + j * n] * x[j];
            }
        }

        return x;
    }

    public static double[] r8ut_slt(int n, double[] a, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8UT_SLT solves a linear system A'*x=b with R8UT matrix.
        //
        //  Discussion:
        //
        //    The R8UT storage format is used for an M by N upper triangular matrix,
        //    and allocates space even for the zero entries.
        //
        //    No factorization of the upper triangular matrix is required.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 August 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[N*N], the R8UT matrix.
        //
        //    Input, double B[N], the right hand side.
        //
        //    Output, double R8UT_SLT[N], the solution vector.
        //
    {
        int i;
        int j;

        double[] x = new double[n];

        for (i = 0; i < n; i++)
        {
            x[i] = b[i];
        }

        for (j = 0; j < n; j++)
        {
            x[j] /= a[j + j * n];
            for (i = j + 1; i < n; i++)
            {
                x[i] -= a[j + i * n] * x[j];
            }
        }

        return x;
    }

    public static double[] r8ut_to_r8ge(int m, int n, double[] a_ut)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8UT_TO_R8GE copies an R8UT matrix to an R8GE matrix.
        //
        //  Discussion:
        //
        //    The R8UT storage format is used for an M by N upper triangular matrix,
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
        //    19 August 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Input, double A_UT[M,N], the R8UT matrix.
        //
        //    Output, double R8UT_TO_R8GE[M,N], the R8GE matrix.
        //
    {
        int j;

        double[] a_ge = new double[m * n];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                if (i <= j)
                {
                    a_ge[i + j * m] = a_ut[i + j * m];
                }
                else
                {
                    a_ge[i + j * m] = 0.0;
                }
            }
        }

        return a_ge;
    }

    public static double[] r8ut_zeros(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8UT_ZEROS zeros an R8UT matrix.
        //
        //  Discussion:
        //
        //    The R8UT storage format is used for an M by N upper triangular matrix,
        //    and allocates space even for the zero entries.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 September 2003
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
        //    Output, double R8UT_ZERO[M*N], the R8UT matrix.
        //
    {
        int j;

        double[] a = new double[m * n];

        for (j = 0; j < n; j++)
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