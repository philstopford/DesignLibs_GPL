using System;
using System.Globalization;
using Burkardt.Uniform;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double r8blt_det(int n, int ml, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BLT_DET computes the determinant of an R8BLT matrix.
        //
        //  Discussion:
        //
        //    The R8BLT storage format is appropriate for a banded lower triangular matrix.
        //    The matrix is assumed to be zero below the ML-th subdiagonal.
        //    The matrix is stored in an ML+1 by N array, in which the diagonal
        //    appears in the first row, followed by successive subdiagonals.
        //    Columns are preserved.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 October 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int ML, the lower bandwidth.
        //
        //    Input, double A[(ML+1)*N], the R8BLT matrix.
        //
        //    Output, double R8BLT_DET, the determinant of A.
        //
    {
        int j;

        double det = 1.0;
        for (j = 0; j < n; j++)
        {
            det *= a[0 + j * (ml + 1)];
        }

        return det;
    }

    public static double[] r8blt_indicator(int n, int ml)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BLT_INDICATOR sets up an R8BLT indicator matrix.
        //
        //  Discussion:
        //
        //    The R8BLT storage format is appropriate for a banded lower triangular matrix.
        //    The matrix is assumed to be zero below the ML-th subdiagonal.
        //    The matrix is stored in an ML+1 by N array, in which the diagonal
        //    appears in the first row, followed by successive subdiagonals.
        //    Columns are preserved.
        //
        //  Example:
        //
        //    N = 5, ML = 2
        //
        //    A11   0   0   0   0
        //    A21 A22   0   0   0
        //    A31 A32 A33   0   0
        //      0 A42 A43 A44   0
        //      0   0 A53 A54 A55
        //                --- ---
        //                    ---
        //
        //    The indicator matrix is stored as:
        //
        //      11  22  33  44  55
        //      21  32  43  54   0
        //      31  42  53   0   0
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
        //    Input, int N, the number of columns of the matrix.
        //
        //    Input, int ML, the lower bandwidth.
        //
        //    Output, double R8BLT_INDICATOR[(ML+1)*N], the R8BLT matrix.
        //
    {
        int i;
        int j;

        double[] a = new double[(ml + 1) * n];

        int fac = (int) Math.Pow(10, Math.Log10(n) + 1);

        for (i = 1; i <= n; i++)
        {
            int jlo = Math.Max(1, i - ml);
            for (j = jlo; j <= i; j++)
            {
                a[i - j + (j - 1) * (ml + 1)] = fac * i + j;
            }
        }

        //
        //  The junk entries can be thought of as corresponding to
        //  elements of a phantom portion of the matrix.
        //
        for (i = n; i < n + ml; i++)
        {
            for (j = i - ml; j < n; j++)
            {
                a[i - j + j * (ml + 1)] = 0.0;
            }
        }

        return a;
    }

    public static double[] r8blt_mtv(int n, int ml, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BLT_MTV multiplies a vector by an R8BLT matrix.
        //
        //  Discussion:
        //
        //    The R8BLT storage format is appropriate for a banded lower triangular matrix.
        //    The matrix is assumed to be zero below the ML-th subdiagonal.
        //    The matrix is stored in an ML+1 by N array, in which the diagonal
        //    appears in the first row, followed by successive subdiagonals.
        //    Columns are preserved.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 October 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int ML, the lower bandwidth.
        //
        //    Input, double A[(ML+1)*N], the R8BLT matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double R8BLT_MTV[N], the product X*A.
        //
    {
        int i;

        double[] b = r8vec_zeros_new(n);

        for (i = 0; i < n; i++)
        {
            int jlo = Math.Max(0, i - ml);
            int j;
            for (j = jlo; j <= i; j++)
            {
                b[j] += x[i] * a[i - j + j * (ml + 1)];
            }
        }

        return b;
    }

    public static double[] r8blt_mv(int n, int ml, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BLT_MV multiplies an R8BLT matrix times a vector.
        //
        //  Discussion:
        //
        //    The R8BLT storage format is appropriate for a banded lower triangular matrix.
        //    The matrix is assumed to be zero below the ML-th subdiagonal.
        //    The matrix is stored in an ML+1 by N array, in which the diagonal
        //    appears in the first row, followed by successive subdiagonals.
        //    Columns are preserved.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 January 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int ML, the lower bandwidth.
        //
        //    Input, double A[(ML+1)*N], the R8BLT matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double R8BLT_MV[N], the product A * x.
        //
    {
        int i;

        double[] b = r8vec_zeros_new(n);

        for (i = 0; i < n; i++)
        {
            int jlo = Math.Max(0, i - ml);
            int j;
            for (j = jlo; j <= i; j++)
            {
                b[i] += a[i - j + j * (ml + 1)] * x[j];
            }
        }

        return b;
    }

    public static void r8blt_print(int n, int ml, double[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BLT_PRINT prints an R8BLT matrix.
        //
        //  Discussion:
        //
        //    The R8BLT storage format is appropriate for a banded lower triangular matrix.
        //    The matrix is assumed to be zero below the ML-th subdiagonal.
        //    The matrix is stored in an ML+1 by N array, in which the diagonal
        //    appears in the first row, followed by successive subdiagonals.
        //    Columns are preserved.
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
        //    Input, int ML, the lower bandwidth.
        //
        //    Input, double A[(ML+1)*N], the R8BLT matrix.
        //
        //    Input, string TITLE, a title.
        //
    {
        r8blt_print_some(n, ml, a, 1, 1, n, n, title);
    }

    public static void r8blt_print_some(int n, int ml, double[] a, int ilo, int jlo, int ihi,
            int jhi, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BLT_PRINT_SOME prints some of an R8BLT matrix.
        //
        //  Discussion:
        //
        //    The R8BLT storage format is appropriate for a banded lower triangular matrix.
        //    The matrix is assumed to be zero below the ML-th subdiagonal.
        //    The matrix is stored in an ML+1 by N array, in which the diagonal
        //    appears in the first row, followed by successive subdiagonals.
        //    Columns are preserved.
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
        //    Input, int ML, the lower bandwidth.
        //
        //    Input, double A[(ML+1)*N], the R8BLT matrix.
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
            i2lo = Math.Max(i2lo, j2lo);
            int i2hi = Math.Min(ihi, n);
            i2hi = Math.Min(i2hi, j2hi + ml);

            int i;
            for (i = i2lo; i <= i2hi; i++)
            {
                cout = i.ToString(CultureInfo.InvariantCulture).PadLeft(5) + "  ";
                //
                //  Print out (up to) 5 entries in row I, that lie in the current strip.
                //
                for (j = j2lo; j <= j2hi; j++)
                {
                    if (ml < i - j || 0 < j - i)
                    {
                        cout += "              ";
                    }
                    else
                    {
                        cout += a[i - j + (j - 1) * (ml + 1)].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "  ";
                    }
                }

                Console.WriteLine(cout);
            }
        }

        Console.WriteLine("");

    }

    public static double[] r8blt_random(int n, int ml, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BLT_RANDOM randomizes an R8BLT matrix.
        //
        //  Discussion:
        //
        //    The R8BLT storage format is appropriate for a banded lower triangular matrix.
        //    The matrix is assumed to be zero below the ML-th subdiagonal.
        //    The matrix is stored in an ML+1 by N array, in which the diagonal
        //    appears in the first row, followed by successive subdiagonals.
        //    Columns are preserved.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 October 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of columns of the matrix.
        //
        //    Input, int ML, the lower bandwidth.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double R8BLT_RANDOM[(ML+1)*N], the R8BLT matrix.
        //
    {
        int i;
        int j;

        double[] a = new double[(ml + 1) * n];

        for (i = 0; i < n; i++)
        {
            int jlo = Math.Max(0, i - ml);
            for (j = jlo; j <= i; j++)
            {
                a[i - j + j * (ml + 1)] = UniformRNG.r8_uniform_01(ref seed);
            }
        }

        //
        //  The junk entries can be thought of as corresponding to
        //  elements of a phantom portion of the matrix.
        //
        for (i = n; i < n + ml; i++)
        {
            for (j = i - ml; j < n; j++)
            {
                a[i - j + j * (ml + 1)] = 0.0;
            }
        }

        return a;
    }

    public static double[] r8blt_sl(int n, int ml, double[] a, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BLT_SL solves A*x=b, where A is an R8BLT matrix.
        //
        //  Discussion:
        //
        //    The R8BLT storage format is appropriate for a banded lower triangular matrix.
        //    The matrix is assumed to be zero below the ML-th subdiagonal.
        //    The matrix is stored in an ML+1 by N array, in which the diagonal
        //    appears in the first row, followed by successive subdiagonals.
        //    Columns are preserved.
        //
        //    No factorization of the lower triangular matrix is required.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int ML, the lower bandwidth.
        //
        //    Input, double A[(ML+1)*N], the R8BLT matrix.
        //
        //    Input, double B(N), the right hand side.
        //
        //    Output, double R8BLT_SL[N], the solution vector.
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
            x[j] /= a[0 + j * (ml + 1)];
            int ihi = Math.Min(j + ml, n - 1);
            for (i = j + 1; i <= ihi; i++)
            {
                x[i] -= a[i - j + j * (ml + 1)] * x[j];
            }
        }

        return x;
    }

    public static double[] r8blt_slt(int n, int ml, double[] a, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BLT_SLT solves A'*x=b, where A is an R8BLT matrix.
        //
        //  Discussion:
        //
        //    The R8BLT storage format is appropriate for a banded lower triangular matrix.
        //    The matrix is assumed to be zero below the ML-th subdiagonal.
        //    The matrix is stored in an ML+1 by N array, in which the diagonal
        //    appears in the first row, followed by successive subdiagonals.
        //    Columns are preserved.
        //
        //    No factorization of the lower triangular matrix is required.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int ML, the lower bandwidth.
        //
        //    Input, double A[(ML+1)*N], the R8BLT matrix.
        //
        //    Input, double B(N), the right hand side.
        //
        //    Output, double R8BLT_SLT[N], the solution vector.
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
            x[j] /= a[0 + j * (ml + 1)];
            int ilo = Math.Max(j - ml, 0);
            for (i = ilo; i <= j - 1; i++)
            {
                x[i] -= a[j - i + i * (ml + 1)] * x[j];
            }
        }

        return x;
    }

    public static double[] r8blt_to_r8ge(int n, int ml, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BLT_TO_R8GE copies an R8BLT matrix to an R8GE matrix.
        //
        //  Discussion:
        //
        //    The R8BLT storage format is used for a banded lower triangular matrix.
        //    The matrix is assumed to be zero below the ML-th subdiagonal.
        //    The matrix is stored in an ML+1 by N array, in which the diagonal
        //    appears in the first row, followed by successive subdiagonals.
        //    Columns are preserved.
        //
        //  Example:
        //
        //    N = 5, ML = 2
        //
        //    A11   0   0   0   0
        //    A21 A22   0   0   0
        //    A31 A32 A33   0   0
        //      0 A42 A43 A44   0
        //      0   0 A53 A54 A55
        //                --- ---
        //                    ---
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 March 2004
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
        //    Input, int ML, the lower bandwidth of A.
        //    ML must be nonnegative, and no greater than N-1.
        //
        //    Input, double A[(ML+1)*N], the R8BLT matrix.
        //
        //    Output, double R8BLT_TO_R8GE[N*N], the R8GE matrix.
        //
    {
        int i;

        double[] b = new double[n * n];

        for (i = 1; i <= n; i++)
        {
            int j;
            for (j = 1; j <= n; j++)
            {
                if (j <= i && i <= j + ml)
                {
                    b[i - 1 + (j - 1) * n] = a[i - j + (j - 1) * (ml + 1)];
                }
                else
                {
                    b[i - 1 + (j - 1) * n] = 0.0;
                }
            }
        }

        return b;
    }

    public static double[] r8blt_zeros(int n, int ml)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BLT_ZEROS zeros an R8BLT matrix.
        //
        //  Discussion:
        //
        //    The R8BLT storage format is appropriate for a banded lower triangular matrix.
        //    The matrix is assumed to be zero below the ML-th subdiagonal.
        //    The matrix is stored in an ML+1 by N array, in which the diagonal
        //    appears in the first row, followed by successive subdiagonals.
        //    Columns are preserved.
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
        //    Input, int N, the number of columns of the matrix.
        //
        //    Input, int ML, the lower bandwidth.
        //
        //    Output, double R8BLT_ZERO[(ML+1)*N], the R8BLT matrix.
        //
    {
        int j;

        double[] a = new double[(ml + 1) * n];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < ml + 1; i++)
            {
                a[i + j * (ml + 1)] = 0.0;
            }
        }

        return a;
    }

}