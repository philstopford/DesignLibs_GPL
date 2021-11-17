using System;
using Burkardt.Uniform;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double[] r8bto_dif2(int m, int l)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BTO_DIF2 sets up an R8BTO second difference matrix.
        //
        //  Discussion:
        //
        //    To get the second difference matrix, it is assumed that M will be 1!
        //
        //    The R8BTO storage format is for a block Toeplitz matrix. The matrix
        //    can be regarded as an L by L array of blocks, each of size M by M.
        //    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
        //    that is, along its diagonal, the blocks repeat.
        //
        //    Storage for the matrix consists of the L blocks of the first row,
        //    followed by the L-1 blocks of the first column (skipping the first row).
        //    These items are stored in the natural way in an (M,M,2*L-1) array.
        //
        //  Example:
        //
        //    M = 2, L = 3
        //
        //    1 2 | 3 4 | 5 6
        //    5 5 | 6 6 | 7 7
        //    ----+-----+-----
        //    7 8 | 1 2 | 3 4
        //    8 8 | 5 5 | 6 6
        //    ----+-----+-----
        //    9 0 | 7 8 | 1 2
        //    9 9 | 8 8 | 5 5
        //
        //    X = (/ 1, 2, 3, 4, 5, 6 /)
        //
        //    B = (/ 91, 134, 73, 125, 97, 129 /)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 July 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the order of the blocks of the matrix A.
        //
        //    Input, int L, the number of blocks in a row or column of A.
        //
        //    Output, double R8BTO_INDICATOR[M*M*(2*L-1)], the R8BTO matrix.
        //
    {
        double[] a;
        int i;
        int i2;
        int j;
        int j2;
        int k;
        double value = 0;

        a = r8vec_zeros_new(m * m * (2 * l - 1));
        //
        //  Blocks 1 to L form the first row.
        //
        j = 0;

        for (k = 1; k <= l; k++)
        {
            value = k switch
            {
                1 => 2.0,
                2 => -1.0,
                _ => 0.0
            };

            for (j2 = 1; j2 <= m; j2++)
            {
                j += 1;
                for (i = 1; i <= m; i++)
                {
                    a[i - 1 + (j2 - 1) * m + (k - 1) * m * m] = value;
                }
            }
        }

        //
        //  Blocks L+1 through 2*L-1 form the remainder of the first column.
        //
        i = m;

        for (k = l + 1; k <= 2 * l - 1; k++)
        {
            if (k == l + 1)
            {
                value = -1.0;
            }
            else
            {
                value = 0.0;
            }

            for (i2 = 1; i2 <= m; i2++)
            {
                i += 1;
                for (j = 1; j <= m; j++)
                {
                    a[i2 - 1 + (j - 1) * m + (k - 1) * m * m] = value;
                }
            }
        }

        return a;
    }

    public static double[] r8bto_indicator(int m, int l)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BTO_INDICATOR sets up an R8BTO indicator matrix.
        //
        //  Discussion:
        //
        //    The R8BTO storage format is for a block Toeplitz matrix. The matrix
        //    can be regarded as an L by L array of blocks, each of size M by M.
        //    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
        //    that is, along its diagonal, the blocks repeat.
        //
        //    Storage for the matrix consists of the L blocks of the first row,
        //    followed by the L-1 blocks of the first column (skipping the first row).
        //    These items are stored in the natural way in an (M,M,2*L-1) array.
        //
        //  Example:
        //
        //    M = 2, L = 3
        //
        //    1 2 | 3 4 | 5 6
        //    5 5 | 6 6 | 7 7
        //    ----+-----+-----
        //    7 8 | 1 2 | 3 4
        //    8 8 | 5 5 | 6 6
        //    ----+-----+-----
        //    9 0 | 7 8 | 1 2
        //    9 9 | 8 8 | 5 5
        //
        //    X = (/ 1, 2, 3, 4, 5, 6 /)
        //
        //    B = (/ 91, 134, 73, 125, 97, 129 /)
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
        //    Input, int M, the order of the blocks of the matrix A.
        //
        //    Input, int L, the number of blocks in a row or column of A.
        //
        //    Output, double R8BTO_INDICATOR[M*M*(2*L-1)], the R8BTO matrix.
        //
    {
        double[] a;
        int fac;
        int i;
        int i2;
        int j;
        int j2;
        int k;

        a = r8vec_zeros_new(m * m * (2 * l - 1));

        fac = (int) Math.Pow(10, (int) Math.Log10(m * l) + 1);
        //
        //  Blocks 1 to L form the first row.
        //
        j = 0;

        for (k = 1; k <= l; k++)
        {
            for (j2 = 1; j2 <= m; j2++)
            {
                j += 1;
                for (i = 1; i <= m; i++)
                {
                    a[i - 1 + (j2 - 1) * m + (k - 1) * m * m] = fac * i + j;
                }
            }
        }

        //
        //  Blocks L+1 through 2*L-1 form the remainder of the first column.
        //
        i = m;

        for (k = l + 1; k <= 2 * l - 1; k++)
        {
            for (i2 = 1; i2 <= m; i2++)
            {
                i += 1;
                for (j = 1; j <= m; j++)
                {
                    a[i2 - 1 + (j - 1) * m + (k - 1) * m * m] = fac * i + j;
                }
            }
        }

        return a;
    }

    public static double[] r8bto_mtv(int m, int l, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BTO_MTV multiplies a vector times an R8BTO matrix.
        //
        //  Discussion:
        //
        //    The R8BTO storage format is for a block Toeplitz matrix. The matrix
        //    can be regarded as an L by L array of blocks, each of size M by M.
        //    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
        //    that is, along its diagonal, the blocks repeat.
        //
        //    Storage for the matrix consists of the L blocks of the first row,
        //    followed by the L-1 blocks of the first column (skipping the first row).
        //    These items are stored in the natural way in an (M,M,2*L-1) array.
        //
        //  Example:
        //
        //    M = 2, L = 3
        //
        //    1 2 | 3 4 | 5 6
        //    5 5 | 6 6 | 7 7
        //    ----+-----+-----
        //    7 8 | 1 2 | 3 4
        //    8 8 | 5 5 | 6 6
        //    ----+-----+-----
        //    9 0 | 7 8 | 1 2
        //    9 9 | 8 8 | 5 5
        //
        //    X = (/ 1, 2, 3, 4, 5, 6 /)
        //
        //    B = (/ 163, 122, 121, 130, 87, 96 /)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 January 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the order of the blocks of the matrix A.
        //
        //    Input, int L, the number of blocks in a row or column of A.
        //
        //    Input, double A[M*M*(2*L-1)], the R8BTO matrix.
        //
        //    Input, double X[M*L], the vector to be multiplied.
        //
        //    Output, double R8BTO_MTV[M*L], the product X * A.
        //
    {
        double[] b;
        int i;
        int i2;
        int j;
        int k;

        b = r8vec_zeros_new(m * l);
        //
        //  Construct the right hand side by blocks.
        //
        for (j = 1; j <= l; j++)
        {
            for (k = 1; k <= j; k++)
            {
                for (i = 1; i <= m; i++)
                {
                    for (i2 = 1; i2 <= m; i2++)
                    {
                        b[i - 1 + (j - 1) * m] += a[i2 - 1 + (i - 1) * m + (j - k) * m * m] *
                                                  x[i2 - 1 + (k - 1) * m];
                    }
                }
            }

            for (k = j + 1; k <= l; k++)
            {
                for (i = 1; i <= m; i++)
                {
                    for (i2 = 1; i2 <= m; i2++)
                    {
                        b[i - 1 + (j - 1) * m] += a[i2 - 1 + (i - 1) * m + (l + k - j - 1) * m * m] *
                                                  x[i2 - 1 + (k - 1) * m];
                    }
                }
            }
        }

        return b;
    }

    public static double[] r8bto_mv(int m, int l, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BTO_MV multiplies an R8BTO matrix times a vector.
        //
        //  Discussion:
        //
        //    The R8BTO storage format is for a block Toeplitz matrix. The matrix
        //    can be regarded as an L by L array of blocks, each of size M by M.
        //    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
        //    that is, along its diagonal, the blocks repeat.
        //
        //    Storage for the matrix consists of the L blocks of the first row,
        //    followed by the L-1 blocks of the first column (skipping the first row).
        //    These items are stored in the natural way in an (M,M,2*L-1) array.
        //
        //  Example:
        //
        //    M = 2, L = 3
        //
        //    1 2 | 3 4 | 5 6
        //    5 5 | 6 6 | 7 7
        //    ----+-----+-----
        //    7 8 | 1 2 | 3 4
        //    8 8 | 5 5 | 6 6
        //    ----+-----+-----
        //    9 0 | 7 8 | 1 2
        //    9 9 | 8 8 | 5 5
        //
        //    X = (/ 1, 2, 3, 4, 5, 6 /)
        //
        //    B = (/ 91, 134, 73, 125, 79, 138 /)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 January 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the order of the blocks of the matrix A.
        //
        //    Input, int L, the number of blocks in a row or column of A.
        //
        //    Input, double A[M*M*(2*L-1)], the R8BTO matrix.
        //
        //    Input, double X[M*L], the vector to be multiplied.
        //
        //    Output, double R8BTO_MV[M*L], the product A * X.
        //
    {
        double[] b;
        int i;
        int i2;
        int j;
        int k;

        b = r8vec_zeros_new(m * l);
        //
        //  Construct the right hand side by blocks.
        //
        for (j = 0; j < l; j++)
        {
            for (k = 0; k <= j - 1; k++)
            {
                for (i = 0; i < m; i++)
                {
                    for (i2 = 0; i2 < m; i2++)
                    {
                        b[i + j * m] += a[i + i2 * m + (l + j - k - 1) * m * m] * x[i2 + k * m];
                    }
                }
            }

            for (k = j; k < l; k++)
            {
                for (i = 0; i < m; i++)
                {
                    for (i2 = 0; i2 < m; i2++)
                    {
                        b[i + j * m] += a[i + i2 * m + (k - j) * m * m] * x[i2 + k * m];
                    }
                }
            }
        }

        return b;
    }

    public static void r8bto_print(int m, int l, double[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BTO_PRINT prints an R8BTO matrix.
        //
        //  Discussion:
        //
        //    The R8BTO storage format is for a block Toeplitz matrix. The matrix
        //    can be regarded as an L by L array of blocks, each of size M by M.
        //    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
        //    that is, along its diagonal, the blocks repeat.
        //
        //    Storage for the matrix consists of the L blocks of the first row,
        //    followed by the L-1 blocks of the first column (skipping the first row).
        //    These items are stored in the natural way in an (M,M,2*L-1) array.
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
        //    Input, int M, the order of the blocks of the matrix A.
        //
        //    Input, int L, the number of blocks in a row or column of A.
        //
        //    Input, double A[M*M*(2*L-1)], the R8BTO matrix.
        //
        //    Input, string TITLE, a title.
        //
    {
        r8bto_print_some(m, l, a, 1, 1, m * l, m * l, title);
    }

    public static void r8bto_print_some(int m, int l, double[] a, int ilo, int jlo, int ihi,
            int jhi, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BTO_PRINT_SOME prints some of an R8BTO matrix.
        //
        //  Discussion:
        //
        //    The R8BTO storage format is for a block Toeplitz matrix. The matrix
        //    can be regarded as an L by L array of blocks, each of size M by M.
        //    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
        //    that is, along its diagonal, the blocks repeat.
        //
        //    Storage for the matrix consists of the L blocks of the first row,
        //    followed by the L-1 blocks of the first column (skipping the first row).
        //    These items are stored in the natural way in an (M,M,2*L-1) array.
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
        //    Input, int M, the order of the blocks of the matrix A.
        //
        //    Input, int L, the number of blocks in a row or column of A.
        //
        //    Input, double A[M*M*(2*L-1)], the R8BTO matrix.
        //
        //    Input, int ILO, JLO, IHI, JHI, designate the first row and
        //    column, and the last row and column to be printed.
        //
        //    Input, string TITLE, a title.
        //
    {
        int INCX = 5;

        int i;
        int i1;
        int i2;
        int i3hi;
        int i3lo;
        int inc;
        int j;
        int j1;
        int j2;
        int j3hi;
        int j3lo;
        int n;
        string cout = "";

        n = m * l;

        Console.WriteLine("");
        Console.WriteLine(title + "");
        //
        //  Print the columns of the matrix, in strips of 5.
        //
        for (j3lo = jlo; j3lo <= jhi; j3lo += INCX)
        {
            j3hi = j3lo + INCX - 1;
            j3hi = Math.Min(j3hi, n);
            j3hi = Math.Min(j3hi, jhi);

            inc = j3hi + 1 - j3lo;

            Console.WriteLine("");
            cout = "  Col: ";
            for (j = j3lo; j <= j3hi; j++)
            {
                cout += j.ToString().PadLeft(7) + "       ";
            }

            Console.WriteLine(cout);
            Console.WriteLine("  Row");
            Console.WriteLine("  ---");
            //
            //  Determine the range of the rows in this strip.
            //
            i3lo = Math.Max(ilo, 1);
            i3hi = Math.Min(ihi, n);

            for (i = i3lo; i <= i3hi; i++)
            {
                cout = i.ToString().PadLeft(4) + "  ";
                //
                //  Print out (up to) 5 entries in row I, that lie in the current strip.
                //
                for (j = j3lo; j <= j3lo + inc - 1; j++)
                {
                    //
                    //  i = M * ( i1 - 1 ) + i2
                    //  j = M * ( j1 - 1 ) + j2
                    //
                    i1 = (i - 1) / m + 1;
                    i2 = i - m * (i1 - 1);
                    j1 = (j - 1) / m + 1;
                    j2 = j - m * (j1 - 1);

                    if (i1 <= j1)
                    {
                        cout += a[i2 - 1 + (j2 - 1) * m + (j1 - i1) * m * m].ToString().PadLeft(12) + "  ";
                    }
                    else
                    {
                        cout += a[i2 - 1 + (j2 - 1) * m + (l - 1 + i1 - j1) * m * m].ToString().PadLeft(12) + "  ";
                    }
                }

                Console.WriteLine(cout);
            }
        }

    }

    public static double[] r8bto_random(int m, int l, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BTO_RANDOM randomizes an R8BTO matrix.
        //
        //  Discussion:
        //
        //    The R8BTO storage format is for a block Toeplitz matrix. The matrix
        //    can be regarded as an L by L array of blocks, each of size M by M.
        //    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
        //    that is, along its diagonal, the blocks repeat.
        //
        //    Storage for the matrix consists of the L blocks of the first row,
        //    followed by the L-1 blocks of the first column (skipping the first row).
        //    These items are stored in the natural way in an (M,M,2*L-1) array.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 January 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the order of the blocks of the matrix A.
        //
        //    Input, int L, the number of blocks in a row or column of A.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double R8BTO_RANDOM[M*M*(2*L-1)], the R8BTO matrix.
        //
    {
        double[] a;
        int i;
        int j;
        int k;

        a = r8vec_zeros_new(m * m * (2 * l - 1));

        for (i = 0; i < m; i++)
        {
            for (j = 0; j < m; j++)
            {
                for (k = 0; k < 2 * l - 1; k++)
                {
                    a[i + j * m + k * m * m] = UniformRNG.r8_uniform_01(ref seed);
                }
            }
        }

        return a;
    }

    public static double[] r8bto_to_r8ge(int m, int l, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BTO_TO_R8GE copies an R8BTO matrix to an R8GE matrix.
        //
        //  Discussion:
        //
        //    The R8BTO storage format is for a block Toeplitz matrix. The matrix
        //    can be regarded as an L by L array of blocks, each of size M by M.
        //    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
        //    that is, along its diagonal, the blocks repeat.
        //
        //    Storage for the matrix consists of the L blocks of the first row,
        //    followed by the L-1 blocks of the first column (skipping the first row).
        //    These items are stored in the natural way in an (M,M,2*L-1) array.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 November 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the order of the blocks of the R8BTO matrix.
        //
        //    Input, int L, the number of blocks in a row or column of the
        //    R8BTO matrix.
        //
        //    Input, double A[M*M*(2*L-1)], the R8BTO matrix.
        //
        //    Output, double R8BTO_TO_R8GE[(M*L)*(M*L)], the R8GE matrix.
        //
    {
        double[] b;
        int i;
        int i1;
        int i2;
        int j;
        int j1;
        int j2;
        int n;

        n = m * l;

        b = r8vec_zeros_new(n * n);

        for (i = 1; i <= n; i++)
        {
            i1 = (i - 1) / m + 1;
            i2 = i - m * (i1 - 1);

            for (j = 1; j <= n; j++)
            {
                j1 = (j - 1) / m + 1;
                j2 = j - m * (j1 - 1);

                if (i1 <= j1)
                {
                    b[i - 1 + (j - 1) * n] = a[i2 - 1 + (j2 - 1) * m + (j1 - i1) * m * m];
                }
                else
                {
                    b[i - 1 + (j - 1) * n] = a[i2 - 1 + (j2 - 1) * m + (l + i1 - j1 - 1) * m * m];
                }
            }
        }

        return b;
    }

    public static double[] r8bto_zeros(int m, int l)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8BTO_ZEROS zeros an R8BTO matrix.
        //
        //  Discussion:
        //
        //    The R8BTO storage format is for a block Toeplitz matrix. The matrix
        //    can be regarded as an L by L array of blocks, each of size M by M.
        //    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
        //    that is, along its diagonal, the blocks repeat.
        //
        //    Storage for the matrix consists of the L blocks of the first row,
        //    followed by the L-1 blocks of the first column (skipping the first row).
        //    These items are stored in the natural way in an (M,M,2*L-1) array.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 July 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the order of the blocks of the matrix A.
        //
        //    Input, int L, the number of blocks in a row or column of A.
        //
        //    Output, double R8BTO_ZERO[M*M*(2*L-1)], the R8BTO matrix.
        //
    {
        double[] a;

        a = r8vec_zeros_new(m * m * (2 * l - 1));

        return a;
    }

}