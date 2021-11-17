using System;
using Burkardt.Uniform;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double r8vm_det(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VM_DET computes the determinant of an R8VM matrix.
        //
        //  Discussion:
        //
        //    The R8VM storage format is used for an M by N Vandermonde matrix.
        //    An M by N Vandermonde matrix is defined by the values in its second
        //    row, which will be written here as X(1:N).  The matrix has a first 
        //    row of 1's, a second row equal to X(1:N), a third row whose entries
        //    are the squares of the X values, up to the M-th row whose entries
        //    are the (M-1)th powers of the X values.  The matrix can be stored
        //    compactly by listing just the values X(1:N).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of rows and columns of the matrix.
        //
        //    Input, double A[N], the R8VM matrix.
        //
        //    Output, double R8VM_DET, the determinant of the matrix.
        //
    {
        int j;

        double det = 1.0;

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = j + 1; i < n; i++)
            {
                det *= (a[i] - a[j]);
            }
        }

        return det;
    }

    public static double[] r8vm_indicator(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VM_INDICATOR returns an R8VM indicator matrix.
        //
        //  Discussion:
        //
        //    The R8VM storage format is used for an M by N Vandermonde matrix.
        //    An M by N Vandermonde matrix is defined by the values in its second
        //    row, which will be written here as X(1:N).  The matrix has a first 
        //    row of 1's, a second row equal to X(1:N), a third row whose entries
        //    are the squares of the X values, up to the M-th row whose entries
        //    are the (M-1)th powers of the X values.  The matrix can be stored
        //    compactly by listing just the values X(1:N).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 August 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns of the matrix.
        //
        //    Output, double R8VM_INDICATOR[N], the R8VM matrix.
        //
    {
        int j;

        double[] a = new double[n];

        for (j = 0; j < n; j++)
        {
            a[j] = j + 1;
        }

        return a;
    }

    public static double r8vm_indicator_det(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VM_INDICATOR_DET returns the determinant of an R8VM indicator matrix.
        //
        //  Discussion:
        //
        //    The R8VM storage format is used for an M by N Vandermonde matrix.
        //    An M by N Vandermonde matrix is defined by the values in its second
        //    row, which will be written here as X(1:N).  The matrix has a first 
        //    row of 1's, a second row equal to X(1:N), a third row whose entries
        //    are the squares of the X values, up to the M-th row whose entries
        //    are the (M-1)th powers of the X values.  The matrix can be stored
        //    compactly by listing just the values X(1:N).
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
        //    Input, int N, the number of rows and columns of 
        //    the matrix.
        //
        //    Output, double R8VM_INDICATOR_DET, the determinant.
        //
    {
        int i;

        double value = 1.0;
        for (i = 0; i < n; i++)
        {
            value *= r8_factorial(i);
        }

        return value;
    }

    public static double[] r8vm_mtv(int m, int n, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VM_MTV multiplies a vector times an R8VM matrix.
        //
        //  Discussion:
        //
        //    The R8VM storage format is used for an M by N Vandermonde matrix.
        //    An M by N Vandermonde matrix is defined by the values in its second
        //    row, which will be written here as X(1:N).  The matrix has a first 
        //    row of 1's, a second row equal to X(1:N), a third row whose entries
        //    are the squares of the X values, up to the M-th row whose entries
        //    are the (M-1)th powers of the X values.  The matrix can be stored
        //    compactly by listing just the values X(1:N).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns of the matrix.
        //
        //    Input, double A[N], the R8VM matrix.
        //
        //    Input, double X[M], the vector to be multiplied by A.
        //
        //    Output, double R8VM_MTV[N], the product A' * x.
        //
    {
        int j;

        double[] b = new double[n];

        for (j = 0; j < n; j++)
        {
            b[j] = 0.0;
            int i;
            for (i = 0; i < m; i++)
            {
                b[j] += i switch
                {
                    0 => x[i],
                    _ => Math.Pow(a[j], i) * x[i]
                };
            }
        }

        return b;
    }

    public static double[] r8vm_mv(int m, int n, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VM_MV multiplies an R8VM matrix times a vector.
        //
        //  Discussion:
        //
        //    The R8VM storage format is used for an M by N Vandermonde matrix.
        //    An M by N Vandermonde matrix is defined by the values in its second
        //    row, which will be written here as X(1:N).  The matrix has a first 
        //    row of 1's, a second row equal to X(1:N), a third row whose entries
        //    are the squares of the X values, up to the M-th row whose entries
        //    are the (M-1)th powers of the X values.  The matrix can be stored
        //    compactly by listing just the values X(1:N).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns of the matrix.
        //
        //    Input, double A[N], the R8VM matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double R8VM_MV[M], the product A * x.
        //
    {
        int i;

        double[] b = new double[m];

        for (i = 0; i < m; i++)
        {
            b[i] = 0.0;
            int j;
            for (j = 0; j < n; j++)
            {
                b[i] += i switch
                {
                    0 => x[j],
                    _ => Math.Pow(a[j], i) * x[j]
                };
            }
        }

        return b;
    }

    public static void r8vm_print(int m, int n, double[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VM_PRINT prints an R8VM matrix.
        //
        //  Discussion:
        //
        //    The R8VM storage format is used for an M by N Vandermonde matrix.
        //    An M by N Vandermonde matrix is defined by the values in its second
        //    row, which will be written here as X(1:N).  The matrix has a first 
        //    row of 1's, a second row equal to X(1:N), a third row whose entries
        //    are the squares of the X values, up to the M-th row whose entries
        //    are the (M-1)th powers of the X values.  The matrix can be stored
        //    compactly by listing just the values X(1:N).
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
        //    Input, int M, N, the number of rows and columns of the matrix.
        //
        //    Input, double A[N], the R8VM matrix.
        //
        //    Input, string TITLE, a title.
        //
    {
        r8vm_print_some(m, n, a, 1, 1, m, n, title);

    }

    public static void r8vm_print_some(int m, int n, double[] a, int ilo, int jlo, int ihi,
            int jhi, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VM_PRINT_SOME prints some of an R8VM matrix.
        //
        //  Discussion:
        //
        //    The R8VM storage format is used for an M by N Vandermonde matrix.
        //    An M by N Vandermonde matrix is defined by the values in its second
        //    row, which will be written here as X(1:N).  The matrix has a first 
        //    row of 1's, a second row equal to X(1:N), a third row whose entries
        //    are the squares of the X values, up to the M-th row whose entries
        //    are the (M-1)th powers of the X values.  The matrix can be stored
        //    compactly by listing just the values X(1:N).
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
        //    Input, int M, N, the number of rows and columns of the matrix.
        //
        //    Input, double A[N], the R8VM matrix.
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
                cout += j.ToString().PadLeft(7) + "       ";
            }

            Console.WriteLine(cout);
            cout = "";
            Console.WriteLine("  Row");
            Console.WriteLine("  ---");
            //
            //  Determine the range of the rows in this strip.
            //
            int i2lo = Math.Max(ilo, 1);
            int i2hi = Math.Min(ihi, m);

            int i;
            for (i = i2lo; i <= i2hi; i++)
            {
                //
                //  Print out (up to) 5 entries in row I, that lie in the current strip.
                //
                cout = "";
                for (j = j2lo; j <= j2hi; j++)
                {
                    double aij = i switch
                    {
                        1 => 1.0,
                        _ => Math.Pow(a[j - 1], i - 1)
                    };

                    cout += aij.ToString().PadLeft(12) + "  ";
                }

                Console.WriteLine(cout);
            }
        }

        Console.WriteLine("");
    }

    public static double[] r8vm_random(int m, int n, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VM_RANDOM randomizes an R8VM matrix.
        //
        //  Discussion:
        //
        //    The R8VM storage format is used for an M by N Vandermonde matrix.
        //    An M by N Vandermonde matrix is defined by the values in its second
        //    row, which will be written here as X(1:N).  The matrix has a first 
        //    row of 1's, a second row equal to X(1:N), a third row whose entries
        //    are the squares of the X values, up to the M-th row whose entries
        //    are the (M-1)th powers of the X values.  The matrix can be stored
        //    compactly by listing just the values X(1:N).
        //
        //    The parameter M is not actually needed by this routine.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 February 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns of the matrix.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double R8VM_RANDOM[N], the R8VM matrix.
        //
    {
        double[] a = UniformRNG.r8vec_uniform_01_new(n, ref seed);

        return a;
    }

    public static void r8vm_sl(int n, double[] a, double[] b, int job, double[] x, ref int info)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VM_SL solves a R8VM system.
        //
        //  Discussion:
        //
        //    The R8VM storage format is used for an M by N Vandermonde matrix.
        //    An M by N Vandermonde matrix is defined by the values in its second
        //    row, which will be written here as X(1:N).  The matrix has a first 
        //    row of 1's, a second row equal to X(1:N), a third row whose entries
        //    are the squares of the X values, up to the M-th row whose entries
        //    are the (M-1)th powers of the X values.  The matrix can be stored
        //    compactly by listing just the values X(1:N).
        //
        //    Vandermonde systems are very close to singularity.  The singularity
        //    gets worse as N increases, and as any pair of values defining
        //    the matrix get close.  Even a system as small as N = 10 will
        //    involve the 9th power of the defining values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 September 2003
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Golub, VanLoan.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Gene Golub, Charles Van Loan,
        //    Matrix Computations,
        //    Third Edition,
        //    Johns Hopkins, 1996.
        //
        //  Parameters:
        //
        //    Input, int N, the number of rows and columns of the matrix.
        //
        //    Input, double A[N], the R8VM matrix.
        //
        //    Input, double B[N], the right hand side.
        //
        //    Input, int JOB, specifies the system to solve.
        //    0, solve A * x = b.
        //    nonzero, solve A' * x = b.
        //
        //    Output, int &INFO.
        //    0, no error.
        //    nonzero, at least two of the values in A are equal.
        //
        //    Output, double X[N], the solution of the linear system.
        //
    {
        int i;
        int j;
        //
        //  Check for explicit singularity.
        //
        info = 0;

        for (j = 0; j < n; j++)
        {
            for (i = j + 1; i < n; i++)
            {
                if (a[i] == a[j])
                {
                    info = 1;
                    return;
                }
            }
        }

        for (i = 0; i < n; i++)
        {
            x[i] = b[i];
        }

        switch (job)
        {
            case 0:
            {
                for (j = 1; j <= n - 1; j++)
                {
                    for (i = n; j + 1 <= i; i--)
                    {
                        x[i - 1] -= a[j - 1] * x[i - 2];
                    }
                }

                for (j = n - 1; 1 <= j; j--)
                {
                    for (i = j + 1; i <= n; i++)
                    {
                        x[i - 1] /= (a[i - 1] - a[i - j - 1]);
                    }

                    for (i = j; i <= n - 1; i++)
                    {
                        x[i - 1] -= x[i];
                    }
                }

                break;
            }
            default:
            {
                for (j = 1; j <= n - 1; j++)
                {
                    for (i = n; j + 1 <= i; i--)
                    {
                        x[i - 1] = (x[i - 1] - x[i - 2]) / (a[i - 1] - a[i - j - 1]);
                    }
                }

                for (j = n - 1; 1 <= j; j--)
                {
                    for (i = j; i <= n - 1; i++)
                    {
                        x[i - 1] -= x[i] * a[j - 1];
                    }
                }

                break;
            }
        }
    }

    public static double[] r8vm_sl_new(int n, double[] a, double[] b, int job, ref int info)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VM_SL_NEW solves a R8VM system.
        //
        //  Discussion:
        //
        //    The R8VM storage format is used for an M by N Vandermonde matrix.
        //    An M by N Vandermonde matrix is defined by the values in its second
        //    row, which will be written here as X(1:N).  The matrix has a first 
        //    row of 1's, a second row equal to X(1:N), a third row whose entries
        //    are the squares of the X values, up to the M-th row whose entries
        //    are the (M-1)th powers of the X values.  The matrix can be stored
        //    compactly by listing just the values X(1:N).
        //
        //    Vandermonde systems are very close to singularity.  The singularity
        //    gets worse as N increases, and as any pair of values defining
        //    the matrix get close.  Even a system as small as N = 10 will
        //    involve the 9th power of the defining values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 September 2003
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Golub, VanLoan.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Gene Golub, Charles Van Loan,
        //    Matrix Computations,
        //    Third Edition,
        //    Johns Hopkins, 1996.
        //
        //  Parameters:
        //
        //    Input, int N, the number of rows and columns of the matrix.
        //
        //    Input, double A[N], the R8VM matrix.
        //
        //    Input, double B[N], the right hand side.
        //
        //    Input, int JOB, specifies the system to solve.
        //    0, solve A * x = b.
        //    nonzero, solve A' * x = b.
        //
        //    Output, int &INFO.
        //    0, no error.
        //    nonzero, at least two of the values in A are equal.
        //
        //    Output, double R8VM_SL_NEW[N], the solution of the linear system.
        //
    {
        int i;
        int j;
        //
        //  Check for explicit singularity.
        //
        info = 0;

        for (j = 0; j < n; j++)
        {
            for (i = j + 1; i < n; i++)
            {
                if (a[i] == a[j])
                {
                    info = 1;
                    return null;
                }
            }
        }

        double[] x = new double[n];

        for (i = 0; i < n; i++)
        {
            x[i] = b[i];
        }

        switch (job)
        {
            case 0:
            {
                for (j = 1; j <= n - 1; j++)
                {
                    for (i = n; j + 1 <= i; i--)
                    {
                        x[i - 1] -= a[j - 1] * x[i - 2];
                    }
                }

                for (j = n - 1; 1 <= j; j--)
                {
                    for (i = j + 1; i <= n; i++)
                    {
                        x[i - 1] /= (a[i - 1] - a[i - j - 1]);
                    }

                    for (i = j; i <= n - 1; i++)
                    {
                        x[i - 1] -= x[i];
                    }
                }

                break;
            }
            default:
            {
                for (j = 1; j <= n - 1; j++)
                {
                    for (i = n; j + 1 <= i; i--)
                    {
                        x[i - 1] = (x[i - 1] - x[i - 2]) / (a[i - 1] - a[i - j - 1]);
                    }
                }

                for (j = n - 1; 1 <= j; j--)
                {
                    for (i = j; i <= n - 1; i++)
                    {
                        x[i - 1] -= x[i] * a[j - 1];
                    }
                }

                break;
            }
        }

        return x;
    }

    public static void r8vm_slt(int n, double[] a, double[] b, ref double[] x, ref int info)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VM_SLT solves A'*x=b, where A is an R8VM matrix.
        //
        //  Discussion:
        //
        //    The R8VM storage format is used for an M by N Vandermonde matrix.
        //    An M by N Vandermonde matrix is defined by the values in its second
        //    row, which will be written here as X(1:N).  The matrix has a first 
        //    row of 1's, a second row equal to X(1:N), a third row whose entries
        //    are the squares of the X values, up to the M-th row whose entries
        //    are the (M-1)th powers of the X values.  The matrix can be stored
        //    compactly by listing just the values X(1:N).
        //
        //    Vandermonde systems are very close to singularity.  The singularity
        //    gets worse as N increases, and as any pair of values defining
        //    the matrix get close.  Even a system as small as N = 10 will
        //    involve the 9th power of the defining values.
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
        //    John Burkardt.
        //
        //  Reference:
        //
        //    Gene Golub, Charles Van Loan,
        //    Matrix Computations,
        //    Third Edition,
        //    Johns Hopkins, 1996.
        //
        //  Parameters:
        //
        //    Input, int N, the number of rows and columns of the matrix.
        //
        //    Input, double A[N], the R8VM matrix.
        //
        //    Input, double B[N], the right hand side.
        //
        //    Output, double X[N], the solution of the linear system.
        //
        //    Output, int *INFO.
        //    0, no error.
        //    nonzero, at least two of the values in A are equal.
        //
    {
        int i;
        int j;
        //
        //  Check for explicit singularity.
        //
        info = 0;

        for (j = 0; j < n; j++)
        {
            for (i = j + 1; i < n; i++)
            {
                if (a[i] == a[j])
                {
                    info = 1;
                    return;
                }
            }
        }

        for (i = 0; i < n; i++)
        {
            x[i] = b[i];
        }

        for (j = 1; j <= n - 1; j++)
        {
            for (i = n; j + 1 <= i; i--)
            {
                x[i - 1] = (x[i - 1] - x[i - 2]) / (a[i - 1] - a[i - j - 1]);
            }
        }

        for (j = n - 1; 1 <= j; j--)
        {
            for (i = j; i <= n - 1; i++)
            {
                x[i - 1] -= x[i] * a[j - 1];
            }
        }

    }

    public static double[] r8vm_slt_new(int n, double[] a, double[] b, ref int info)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VM_SLT_NEW solves A'*x = b, where A is an R8VM matrix.
        //
        //  Discussion:
        //
        //    The R8VM storage format is used for an M by N Vandermonde matrix.
        //    An M by N Vandermonde matrix is defined by the values in its second
        //    row, which will be written here as X(1:N).  The matrix has a first 
        //    row of 1's, a second row equal to X(1:N), a third row whose entries
        //    are the squares of the X values, up to the M-th row whose entries
        //    are the (M-1)th powers of the X values.  The matrix can be stored
        //    compactly by listing just the values X(1:N).
        //
        //    Vandermonde systems are very close to singularity.  The singularity
        //    gets worse as N increases, and as any pair of values defining
        //    the matrix get close.  Even a system as small as N = 10 will
        //    involve the 9th power of the defining values.
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
        //    John Burkardt.
        //
        //  Reference:
        //
        //    Gene Golub, Charles Van Loan,
        //    Matrix Computations,
        //    Third Edition,
        //    Johns Hopkins, 1996.
        //
        //  Parameters:
        //
        //    Input, int N, the number of rows and columns of the matrix.
        //
        //    Input, double A[N], the R8VM matrix.
        //
        //    Input, double B[N], the right hand side.
        //
        //    Output, int *INFO.
        //    0, no error.
        //    nonzero, at least two of the values in A are equal.
        //
        //    Output, double R8VM_SLT_NEW[N], the solution of the linear system.
        //
    {
        int i;
        int j;
        //
        //  Check for explicit singularity.
        //
        info = 0;

        for (j = 0; j < n; j++)
        {
            for (i = j + 1; i < n; i++)
            {
                if (a[i] == a[j])
                {
                    info = 1;
                    return null;
                }
            }
        }

        double[] x = new double[n];

        for (i = 0; i < n; i++)
        {
            x[i] = b[i];
        }

        for (j = 1; j <= n - 1; j++)
        {
            for (i = n; j + 1 <= i; i--)
            {
                x[i - 1] = (x[i - 1] - x[i - 2]) / (a[i - 1] - a[i - j - 1]);
            }
        }

        for (j = n - 1; 1 <= j; j--)
        {
            for (i = j; i <= n - 1; i++)
            {
                x[i - 1] -= x[i] * a[j - 1];
            }
        }

        return x;
    }

    public static double[] r8vm_to_r8ge(int m, int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VM_TO_R8GE copies an R8VM matrix to an R8GE matrix.
        //
        //  Discussion:
        //
        //    The R8VM storage format is used for an M by N Vandermonde matrix.
        //    An M by N Vandermonde matrix is defined by the values in its second
        //    row, which will be written here as X(1:N).  The matrix has a first 
        //    row of 1's, a second row equal to X(1:N), a third row whose entries
        //    are the squares of the X values, up to the M-th row whose entries
        //    are the (M-1)th powers of the X values.  The matrix can be stored
        //    compactly by listing just the values X(1:N).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 January 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns of the matrix.
        //
        //    Input, double A[N], the R8VM matrix.
        //
        //    Output, double R8VM_TO_R8GE[M*N], the R8GE matrix.
        //
    {
        int i;

        double[] b = new double[m * n];

        for (i = 0; i < m; i++)
        {
            int j;
            for (j = 0; j < n; j++)
            {
                b[i + j * m] = i switch
                {
                    0 => 1.0,
                    _ => b[i - 1 + j * m] * a[j]
                };
            }
        }

        return b;
    }

    public static double[] r8vm_zeros(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VM_ZEROS zeros an R8VM matrix.
        //
        //  Discussion:
        //
        //    The R8VM storage format is used for an M by N Vandermonde matrix.
        //    An M by N Vandermonde matrix is defined by the values in its second
        //    row, which will be written here as X(1:N).  The matrix has a first 
        //    row of 1's, a second row equal to X(1:N), a third row whose entries
        //    are the squares of the X values, up to the M-th row whose entries
        //    are the (M-1)th powers of the X values.  The matrix can be stored
        //    compactly by listing just the values X(1:N).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 August 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns of the matrix.
        //
        //    Output, double R8VM_ZERO[N], the zero R8VM matrix.
        //
    {
        int j;

        double[] a = new double[n];

        for (j = 0; j < n; j++)
        {
            a[j] = 0.0;
        }

        return a;
    }

}