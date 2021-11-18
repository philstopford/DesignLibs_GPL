using System;
using System.Globalization;
using Burkardt.Uniform;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double[] r85_dif2(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R85_DIF2 sets up an R85 second difference matrix.
        //
        //  Discussion:
        //
        //    The R85 storage format represents a pentadiagonal matrix as a 5
        //    by N array, in which each row corresponds to a diagonal, and
        //    column locations are preserved.  Thus, the original matrix is
        //    "collapsed" vertically into the array.
        //
        //  Example:
        //
        //    Here is how an R85 matrix of order 6 would be stored:
        //
        //       *   *  A13 A24 A35 A46
        //       *  A12 A23 A34 A45 A56
        //      A11 A22 A33 A44 A55 A66
        //      A21 A32 A43 A54 A65  *
        //      A31 A42 A53 A64  *   *
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 July 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Output, double R85_DIF2[5*N], the R85 matrix.
        //
    {
        int j;

        double[] a = r8vec_zeros_new(5 * n);

        for (j = 2; j <= n; j++)
        {
            a[1 + (j - 1) * 5] = -1.0;
        }

        for (j = 1; j <= n; j++)
        {
            a[2 + (j - 1) * 5] = 2.0;
        }

        for (j = 1; j <= n - 1; j++)
        {
            a[3 + (j - 1) * 5] = -1.0;
        }

        return a;
    }

    public static double[] r85_indicator(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R85_INDICATOR sets up an R85 indicator matrix.
        //
        //  Discussion:
        //
        //    The R85 storage format represents a pentadiagonal matrix as a 5
        //    by N array, in which each row corresponds to a diagonal, and
        //    column locations are preserved.  Thus, the original matrix is
        //    "collapsed" vertically into the array.
        //
        //  Example:
        //
        //    Here is how an R85 matrix of order 6 would be stored:
        //
        //       *   *  A13 A24 A35 A46
        //       *  A12 A23 A34 A45 A56
        //      A11 A22 A33 A44 A55 A66
        //      A21 A32 A43 A54 A65  *
        //      A31 A42 A53 A64  *   *
        //
        //    Here are the values as stored in an indicator matrix:
        //
        //      00 00 13 24 35 46
        //      00 12 23 34 45 56
        //      11 22 33 44 55 66
        //      21 32 43 54 65 00
        //      31 42 53 64 00 00
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 January 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be at least 2.
        //
        //    Output, double R85_INDICATOR[3*N], the R85 indicator matrix.
        //
    {
        int i;
        int j;

        double[] a = r8vec_zeros_new(5 * n);

        int fac = (int) Math.Pow(10, (int) Math.Log10(n) + 1);

        for (j = 3; j <= n; j++)
        {
            i = j - 2;
            a[0 + (j - 1) * 5] = fac * i + j;
        }

        for (j = 2; j <= n; j++)
        {
            i = j - 1;
            a[1 + (j - 1) * 5] = fac * i + j;
        }

        for (j = 1; j <= n; j++)
        {
            i = j;
            a[2 + (j - 1) * 5] = fac * i + j;
        }

        for (j = 1; j <= n - 1; j++)
        {
            i = j + 1;
            a[3 + (j - 1) * 5] = fac * i + j;
        }

        for (j = 1; j <= n - 2; j++)
        {
            i = j + 2;
            a[4 + (j - 1) * 5] = fac * i + j;
        }

        return a;
    }

    public static double[] r85_np_fs(int n, double[] a, ref double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R85_NP_FS factors and solves an R85 system.
        //
        //  Discussion:
        //
        //    The R85 storage format represents a pentadiagonal matrix as a 5
        //    by N array, in which each row corresponds to a diagonal, and
        //    column locations are preserved.  Thus, the original matrix is
        //    "collapsed" vertically into the array.
        //
        //    This algorithm requires that each diagonal entry be nonzero.
        //
        //    No pivoting is performed, and therefore the algorithm may fail
        //    in simple cases where the matrix is not singular.
        //
        //  Example:
        //
        //    Here is how an R85 matrix of order 6 would be stored:
        //
        //       *   *  A13 A24 A35 A46
        //       *  A12 A23 A34 A45 A56
        //      A11 A22 A33 A44 A55 A66
        //      A21 A32 A43 A54 A65  *
        //      A31 A42 A53 A64  *   *
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 September 2003
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Cheney, Kincaid.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Ward Cheney, David Kincaid,
        //    Numerical Mathematics and Computing,
        //    1985, pages 233-236.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the linear system.
        //
        //    Input/output, double A[5*N],
        //    On input, the pentadiagonal matrix.
        //    On output, the data in these vectors has been overwritten
        //    by factorization information.
        //
        //    Input/output, double B[N].
        //    On input, B contains the right hand side of the linear system.
        //    On output, B has been overwritten by factorization information.
        //
        //    Output, double R85_NP_FS[N], the solution of the linear system.
        //
    {
        int i;
        double xmult;

        for (i = 0; i < n; i++)
        {
            switch (a[2 + i * 5])
            {
                case 0.0:
                    return null;
            }
        }

        double[] x = r8vec_zeros_new(n);

        for (i = 2; i <= n - 1; i++)
        {
            xmult = a[1 + (i - 1) * 5] / a[2 + (i - 2) * 5];
            a[2 + (i - 1) * 5] -= xmult * a[3 + (i - 2) * 5];
            a[3 + (i - 1) * 5] -= xmult * a[4 + (i - 2) * 5];

            b[i - 1] -= xmult * b[i - 2];

            xmult = a[0 + i * 5] / a[2 + (i - 2) * 5];
            a[1 + i * 5] -= xmult * a[3 + (i - 2) * 5];
            a[2 + i * 5] -= xmult * a[4 + (i - 2) * 5];

            b[i] -= xmult * b[i - 2];
        }

        xmult = a[1 + (n - 1) * 5] / a[2 + (n - 2) * 5];
        a[2 + (n - 1) * 5] -= xmult * a[3 + (n - 2) * 5];

        x[n - 1] = (b[n - 1] - xmult * b[n - 2]) / a[2 + (n - 1) * 5];
        x[n - 2] = (b[n - 2] - a[3 + (n - 2) * 5] * x[n - 1]) / a[2 + (n - 2) * 5];

        for (i = n - 2; 1 <= i; i--)
        {
            x[i - 1] = (b[i - 1] - a[3 + (i - 1) * 5] * x[i] - a[4 + (i - 1) * 5] * x[i + 1])
                       / a[2 + (i - 1) * 5];
        }

        return x;
    }

    public static double[] r85_mtv(int n, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R85_MTV multiplies a vector times an R85 matrix.
        //
        //  Discussion:
        //
        //    The R85 storage format represents a pentadiagonal matrix as a 5
        //    by N array, in which each row corresponds to a diagonal, and
        //    column locations are preserved.  Thus, the original matrix is
        //    "collapsed" vertically into the array.
        //
        //  Example:
        //
        //    Here is how an R85 matrix of order 6 would be stored:
        //
        //       *   *  A13 A24 A35 A46
        //       *  A12 A23 A34 A45 A56
        //      A11 A22 A33 A44 A55 A66
        //      A21 A32 A43 A54 A65  *
        //      A31 A42 A53 A64  *   *
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the linear system.
        //
        //    Input, double A[5*N], the pentadiagonal matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A'.
        //
        //    Output, double R85_MTV[N], the product A' * x.
        //
    {
        int j;

        double[] b = r8vec_zeros_new(n);

        for (j = 0; j < n; j++)
        {
            b[j] = a[2 + j * 5] * x[j];
        }

        for (j = 1; j < n; j++)
        {
            b[j] += a[3 + j * 5] * x[j - 1];
        }

        for (j = 2; j < n; j++)
        {
            b[j] += a[4 + j * 5] * x[j - 2];
        }

        for (j = 0; j < n - 1; j++)
        {
            b[j] += a[1 + j * 5] * x[j + 1];
        }

        for (j = 0; j < n - 2; j++)
        {
            b[j] += a[0 + j * 5] * x[j + 2];
        }

        return b;
    }

    public static double[] r85_mv(int n, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R85_MV multiplies an R85 matrix times a vector.
        //
        //  Discussion:
        //
        //    The R85 storage format represents a pentadiagonal matrix as a 5
        //    by N array, in which each row corresponds to a diagonal, and
        //    column locations are preserved.  Thus, the original matrix is
        //    "collapsed" vertically into the array.
        //
        //  Example:
        //
        //    Here is how an R85 matrix of order 6 would be stored:
        //
        //       *   *  A13 A24 A35 A46
        //       *  A12 A23 A34 A45 A56
        //      A11 A22 A33 A44 A55 A66
        //      A21 A32 A43 A54 A65  *
        //      A31 A42 A53 A64  *   *
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the linear system.
        //
        //    Input, double A[5*N], the pentadiagonal matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double R85_MV[N], the product A * x.
        //
    {
        int i;

        double[] b = r8vec_zeros_new(n);

        for (i = 0; i < n; i++)
        {
            b[i] = a[2 + i * 5] * x[i];
        }

        for (i = 2; i < n; i++)
        {
            b[i] += a[0 + i * 5] * x[i - 2];
        }

        for (i = 1; i < n; i++)
        {
            b[i] += a[1 + i * 5] * x[i - 1];
        }

        for (i = 0; i < n - 1; i++)
        {
            b[i] += a[3 + i * 5] * x[i + 1];
        }

        for (i = 0; i < n - 2; i++)
        {
            b[i] += a[4 + i * 5] * x[i + 2];
        }

        return b;
    }

    public static void r85_print(int n, double[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R85_PRINT prints an R85 matrix.
        //
        //  Discussion:
        //
        //    The R85 storage format represents a pentadiagonal matrix as a 5
        //    by N array, in which each row corresponds to a diagonal, and
        //    column locations are preserved.  Thus, the original matrix is
        //    "collapsed" vertically into the array.
        //
        //  Example:
        //
        //    Here is how an R85 matrix of order 6 would be stored:
        //
        //       *   *  A13 A24 A35 A46
        //       *  A12 A23 A34 A45 A56
        //      A11 A22 A33 A44 A55 A66
        //      A21 A32 A43 A54 A65  *
        //      A31 A42 A53 A64  *   *
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
        //    Input, double A[5*N], the pentadiagonal matrix.
        //
        //    Input, string TITLE, a title.
        //
    {
        r85_print_some(n, a, 1, 1, n, n, title);
    }

    public static void r85_print_some(int n, double[] a, int ilo, int jlo, int ihi, int jhi,
            string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R85_PRINT_SOME prints some of an R85 matrix.
        //
        //  Discussion:
        //
        //    The R85 storage format represents a pentadiagonal matrix as a 5
        //    by N array, in which each row corresponds to a diagonal, and
        //    column locations are preserved.  Thus, the original matrix is
        //    "collapsed" vertically into the array.
        //
        //  Example:
        //
        //    Here is how an R85 matrix of order 6 would be stored:
        //
        //       *   *  A13 A24 A35 A46
        //       *  A12 A23 A34 A45 A56
        //      A11 A22 A33 A44 A55 A66
        //      A21 A32 A43 A54 A65  *
        //      A31 A42 A53 A64  *   *
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
        //    Input, double A[5*N], the pentadiagonal matrix.
        //
        //    Input, int ILO, JLO, IHI, JHI, designate the first row and
        //    column, and the last row and column, to be printed.
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

            int inc = j2hi + 1 - j2lo;

            Console.WriteLine("");
            string cout = "  Col:  ";
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
            i2lo = Math.Max(i2lo, j2lo - 2);

            int i2hi = Math.Min(ihi, n);
            i2hi = Math.Min(i2hi, j2hi + 2);

            int i;
            for (i = i2lo; i <= i2hi; i++)
            {
                //
                //  Print out (up to) 5 entries in row I, that lie in the current strip.
                //
                cout = i.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  ";

                int j2;
                for (j2 = 1; j2 <= inc; j2++)
                {
                    j = j2lo - 1 + j2;

                    if (2 < i - j || 2 < j - i)
                    {
                        cout += "            ";
                    }
                    else if (j == i + 2)
                    {
                        cout += a[0 + (j - 1) * 5].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  ";
                    }
                    else if (j == i + 1)
                    {
                        cout += a[1 + (j - 1) * 5].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  ";
                    }
                    else if (j == i)
                    {
                        cout += a[2 + (j - 1) * 5].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  ";
                    }
                    else if (j == i - 1)
                    {
                        cout += a[3 + (j - 1) * 5].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  ";
                    }
                    else if (j == i - 2)
                    {
                        cout += a[4 + (j - 1) * 5].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  ";
                    }
                }

                Console.WriteLine(cout);
            }

            Console.WriteLine("");
        }

    }

    public static double[] r85_random(int n, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R85_RANDOM randomizes an R85 matrix.
        //
        //  Discussion:
        //
        //    The R85 storage format represents a pentadiagonal matrix as a 5
        //    by N array, in which each row corresponds to a diagonal, and
        //    column locations are preserved.  Thus, the original matrix is
        //    "collapsed" vertically into the array.
        //
        //  Example:
        //
        //    Here is how an R85 matrix of order 6 would be stored:
        //
        //       *   *  A13 A24 A35 A46
        //       *  A12 A23 A34 A45 A56
        //      A11 A22 A33 A44 A55 A66
        //      A21 A32 A43 A54 A65  *
        //      A31 A42 A53 A64  *   *
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 March 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the linear system.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double R85_RANDOM[5*N], the pentadiagonal matrix.
        //
    {
        int j;

        double[] a = r8vec_zeros_new(5 * n);

        int i = 0;
        for (j = 2; j < n; j++)
        {
            a[i + j * 5] = UniformRNG.r8_uniform_01(ref seed);
        }

        i = 1;
        for (j = 1; j < n; j++)
        {
            a[i + j * 5] = UniformRNG.r8_uniform_01(ref seed);
        }

        i = 2;
        for (j = 0; j < n; j++)
        {
            a[i + j * 5] = UniformRNG.r8_uniform_01(ref seed);
        }

        i = 3;
        for (j = 0; j < n - 1; j++)
        {
            a[i + j * 5] = UniformRNG.r8_uniform_01(ref seed);
        }

        i = 4;
        for (j = 0; j < n - 2; j++)
        {
            a[i + j * 5] = UniformRNG.r8_uniform_01(ref seed);
        }

        return a;
    }

    public static double[] r85_to_r8ge(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R85_TO_R8GE copies an R85 matrix to an R8GE matrix.
        //
        //  Discussion:
        //
        //    The R85 storage format represents a pentadiagonal matrix as a 5
        //    by N array, in which each row corresponds to a diagonal, and
        //    column locations are preserved.  Thus, the original matrix is
        //    "collapsed" vertically into the array.
        //
        //  Example:
        //
        //    Here is how an R85 matrix of order 6 would be stored:
        //
        //       *   *  A13 A24 A35 A46
        //       *  A12 A23 A34 A45 A56
        //      A11 A22 A33 A44 A55 A66
        //      A21 A32 A43 A54 A65  *
        //      A31 A42 A53 A64  *   *
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be at least 3.
        //
        //    Input, double A[5*N], the nonzero diagonals of the matrix.
        //
        //    Output, double R85_TO_R8GE[N*N], the pentadiagonal matrix.
        //
    {
        int j;

        double[] b = r8vec_zeros_new(n * n);

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < n; i++)
            {
                if (j == i - 2)
                {
                    b[i + j * 5] = a[0 + i * 5];
                }
                else if (j == i - 1)
                {
                    b[i + j * 5] = a[1 + i * 5];
                }
                else if (i == j)
                {
                    b[i + j * 5] = a[2 + i * 5];
                }
                else if (j == i + 1)
                {
                    b[i + j * 5] = a[3 + i * 5];
                }
                else if (j == i + 2)
                {
                    b[i + j * 5] = a[4 + i * 5];
                }
                else
                {
                    b[i + j * 5] = 0.0;
                }
            }
        }

        return b;
    }

    public static double[] r85_zeros(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R85_ZEROS zeros an R85 matrix.
        //
        //  Discussion:
        //
        //    The R85 storage format represents a pentadiagonal matrix as a 5
        //    by N array, in which each row corresponds to a diagonal, and
        //    column locations are preserved.  Thus, the original matrix is
        //    "collapsed" vertically into the array.
        //
        //  Example:
        //
        //    Here is how an R85 matrix of order 6 would be stored:
        //
        //       *   *  A13 A24 A35 A46
        //       *  A12 A23 A34 A45 A56
        //      A11 A22 A33 A44 A55 A66
        //      A21 A32 A43 A54 A65  *
        //      A31 A42 A53 A64  *   *
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
        //    Input, int N, the order of the linear system.
        //
        //    Output, double R85_ZERO[5*N], the R85 matrix.
        //
    {
        double[] a = r8vec_zeros_new(5 * n);

        return a;
    }

}